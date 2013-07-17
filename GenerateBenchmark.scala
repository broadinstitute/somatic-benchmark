package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.Logging

import scala.collection.immutable.Map
import org.broadinstitute.sting.queue.extensions.gatk.{CommandLineGATK, PrintReads}
import org.broadinstitute.sting.queue.extensions.picard.SortSam

class GenerateBenchmark extends QScript with Logging {
  qscript =>
  //TODO implement these as cmdline paramaters isntead of hard coding them
  val indelFile : File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf")
  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bam : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
  val spikeContributorBAM : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")

  val libDir : File = new File(".")
  
  val intervalFile = new File(libDir,"chr20.interval_list" )
  val spikeSitesVCF = new File(libDir,"vcf_data/na12878_ref_NA12891_het_chr1_high_conf.vcf" )

  val gatk : File = new File("/xchip/cga2/louisb/gatk-protected/dist/GenomeAnalysisTK.jar")
  val prefix = "chr1"
  val validateIndelsFile = "indels.vcf"

  val PICARD_PATH = "/seq/software/picard/current/bin"
  val sortSamPath = new File(PICARD_PATH, "SortSam.jar")
  val mergeSamPath = new File(PICARD_PATH, "MergeSamFiles.jar")
  val tmpdir = "/broad/hptmp/louisb/sim"
   
  val outDir = new File(libDir)

  //TODOAn ugly hardcoded hack.  Must eventually be replaced when the number of divisions while fracturing is allowed to be changed from 6.
  val bamMapFile = new File(libDir,"louis_bam_1g_info.txt" )
  var bamDigitToNameMap :Map[Char,String]= null 
  var bamNameToFileMap : Map[String,File]= null 
  
  
  //Used by make_fn_data 
  val alleleFractions = Set( 0.04, .1 , .2, .4, .8)
  val maxDepth = "123456789ABC"
  val depths = for(i <- 1 to maxDepth.length) yield maxDepth.substring(0, i) 
 
  val PIECES = 6;
  val LIBRARIES = List("Solexa-18483","Solexa-18484","Solexa-23661")



  def script()= {
    qscript.bamDigitToNameMap = loadBamMap(bamMapFile)
   
    //make vcfs 
    val makeVcfs = new MakeVcfs
    val vcfOut = new File(libDir)  
    makeVcfs.indelFile = qscript.indelFile  
    makeVcfs.vcfOutFile = spikeSitesVCF 
    add(makeVcfs)
    
    //fracture bams 
    val fractureOutDir = new File(outDir,"data_1g_wgs" )
    val (splitBams, fractureCmds) = FractureBams.makeFractureJobs(bam, referenceFile, LIBRARIES, intervalFile, PIECES, fractureOutDir)
    fractureCmds.foreach(add(_)) 
    
    qscript.bamNameToFileMap = splitBams.map( (bam: File) => (bam.getName, bam)).toMap
    
    //use SomaticSpike to create false negative test data 
    val makeFnCommands = new FalseNegativeSim(spikeSitesVCF,spikeContributorBAM)
    val falseNegativeCmds = makeFnCommands.makeFnSimCmds( alleleFractions, depths)
    falseNegativeCmds.foreach(add(_))

    //merge bams
    val mergers = MergeBams.makeMergeBamsJobs(fractureOutDir)
    mergers.foreach(add(_)) 

  }

  trait JobQueueArguments extends CommandLineFunction {
    this.jobQueue = "week"

  }

  trait GeneratorArguments extends CommandLineGATK with JobQueueArguments{
    this.reference_sequence = referenceFile
  }

  class FilterByLibrary extends PrintReads with GeneratorArguments{
    @Argument(doc="library name")
    var library :String = _

    this.javaMemoryLimit = 2
    this.read_filter ++= List("DuplicateRead", "FailsVendorQualityCheck","UnmappedRead")

    override def commandLine = super.commandLine + required("-rf","LibraryRead") + required("--library",library)

  }

  object FractureBams {
      val FILE_NAME_PREFIX = "NA12878.WGS" 
      class NameSortBamByLibrary extends CommandLineFunction with JobQueueArguments{
        @Input(doc="bam file to sort")
        var inputBam : File = _
        
        @Input(doc="interval file")
        var interval : File = qscript.intervalFile 
        
        
        @Argument(doc="library name")
        var library :String = _

        @Output(doc="sorted bam of only the reads from one library")
        var nameSortedBam : File = _
        this.memoryLimit = 18

        lazy val filter = new FilterByLibrary {
          this.library = library
          this.input_file :+= inputBam
          this.out = nameSortedBam
        }

        lazy val sort = new SortSam with JobQueueArguments {
          this.javaMemoryLimit = 16
          this.maxRecordsInRam = 4000000
          this.input :+= "/dev/stdin"
          this.output = nameSortedBam
          this.sortOrder = net.sf.samtools.SAMFileHeader.SortOrder.queryname
          this.compressionLevel = 1
        }


        def commandLine = filter.commandLine + required("|",escape=false) + sort.commandLine
      } 
  
       class SplitBam extends CommandLineFunction{
          @Input(doc="bam file to copy header from")
          var headerBam: File = _

          @Input(doc="bam to be split, the one filtered by library")
          var nameSortedBam: File = _

          @Output(doc="list of output files")
          var outFiles: List[File] = _

          def commandLine = required("%s/splitBam.pl".format(libDir))+
                            required(headerBam)+
                            required(nameSortedBam)+
                            repeat(outFiles)
      }
      
      class CoordinateSortAndConvertToBAM extends CommandLineFunction {
        @Input(doc="input Sam files")
        var inputSam: File = _

        @Output(doc="output Bam files")
        var outputBam: File = _
        
        def commandLine = required("java")+
                          required("-Xmx2g")+
                          required("-jar",qscript.sortSamPath)+
                          required("TMP_DIR=",qscript.tmpdir, spaceSeparated=false)+
                          required("I=",inputSam, spaceSeparated=false)+
                          required("O=",outputBam, spaceSeparated=false)+
                          required("SO=","coordinate", spaceSeparated=false)+
                          required("CREATE_INDEX=","true", spaceSeparated=false)+
                          required("QUIET=","true", spaceSeparated=false)    
      }
      def makeFractureJobs(bam: File, reference: File,  libraries: Traversable[String], interval : File, pieces: Int, outDir : File) = {
        
        def makeSingleFractureJob(library: String):(List[File],List[CommandLineFunction])={
          def getSplitBamNames(library:String, pieces:Int):Traversable[String] ={
            val outmask = FILE_NAME_PREFIX+".somatic.simulation.%s.%03d.sam"
            for(i <- 1 to pieces) yield outmask.format(library, i)
          }
           
          val sort = new NameSortBamByLibrary
          sort.reference = reference
          sort.inputBam = bam
          sort.interval = interval 
          sort.library = library
          
          val sortedBam = new File(outDir,FILE_NAME_PREFIX+".original.regional.namesorted.%s.bam".format(library))
          sort.nameSortedBam = sortedBam
          
          
          val split = new SplitBam
          split.headerBam = bam
          split.nameSortedBam = sortedBam

          val splitSams :List[File]= getSplitBamNames(library,pieces).map( new File(outDir, _)).toList
          split.outFiles = splitSams
          
          val splitBams: List[File] = splitSams.map((sam:File) =>( swapExt(outDir, sam ,"sam", "bam" )))
            
          val converters = (splitSams, splitBams).zipped map {(samFile, outputBam) =>
              val convert = new CoordinateSortAndConvertToBAM
              convert.inputSam = samFile
              convert.outputBam = outputBam
              convert
          }
           
          (splitBams, sort::split::converters)
      }

      val (splitBams, cmds)= (for ( library <- libraries) yield ( makeSingleFractureJob(library)) ).unzip
      
      (splitBams.flatten, cmds.flatten)
   }  
  }

  class MergeBams extends CommandLineFunction {
        @Input(doc="list of files to merge")
        var toMerge: List[File] = Nil 

        @Output(doc="output bam file")
        var merged: File = _

        @Output(doc="merged bam index file")
        var mergedBai: File = _
        def commandLine = required("java") +
                          required("-Xmx2g") +
                          required("-jar",mergeSamPath) + 
                          required("CREATE_INDEX=", true, spaceSeparated=false) +
                          required("USE_THREADING=", true, spaceSeparated=false) +
                          required("O=", merged, spaceSeparated=false) + 
                          repeat("I=", toMerge, spaceSeparated=false) +
                          required("&&",escape=false) +
                          required("cp", mergedBai, merged+".bai" )
  } 

  object MergeBams {
    private val outFileNameTemplate = "NA12878.somatic.simulation.merged.%s.bam"
    private val BAMGROUPS = List( 
             "123456789ABC",
             "123456789AB",
             "123456789A",
             "123456789",
             "12345678",
             "1234567",
             "123456",
             "12345",
             "1234",
             "123",
             "12",
             "1",
             "DEFGHI",
             "DEFGH",
             "DEFG", 
             "DEF",
             "DE",
             "FG",
             "HI",
             "D"
           )
       def makeMergeBamsJobs(dir: File) = {
        BAMGROUPS.map{ name =>
            val mergedFile = new File(dir, outFileNameTemplate.format(name) )
            val inputBams = getBams(name)
            val merge = new MergeBams
            val mergedBai = swapExt(dir, mergedFile, "bam", "bai")
            merge.merged = mergedFile
            merge.toMerge = inputBams
            merge.mergedBai = mergedBai
            merge
           }
        }

    }


  class SomaticSpike extends CommandLineFunction{ 
   @Input(doc="spike location intervals file")
   var spikeSitesVCF: File = _

   @Input(doc="tumorBams")
   var tumorBams: List[File]= Nil
   
   @Argument(doc="tumor allele fraction to simulate", required=false)
   var alleleFraction : Double = _

   @Output(doc="interval file of all the spiked in locations")
   var outIntervals: File = _

   @Output(doc="output Bam")
   var outBam : File = _

   def commandLine = required("java")+
                          required("-Xmx4g")+
                          required("-jar", qscript.gatk)+
                          required("-T","SomaticSpike")+
                          required("--reference_sequence", qscript.referenceFile)+
                          repeat("-I", tumorBams)+
                          required("-I:spike", qscript.spikeContributorBAM)+
                          required("--intervals",spikeSitesVCF)+
                          optional("--simulation_fraction", alleleFraction)+
                          optional("--minimum_qscore", 20)+
                          optional("--out", outBam)+
                          optional("--spiked_intervals_out", outIntervals) 
  }
  
 
  class MakeVcfs extends CommandLineFunction{
    @Input(doc="vcf file containing indels to use as true indel sites") var indelFile : File = _
    @Output(doc="dummy output for queue ordering") var vcfOutFile : File = _ 
    def commandLine = "%s/make_vcfs.pl %s".format(libDir, indelFile)
  }

  class FalseNegativeSim(spikeSitesVCF : File, spikeInBam : File) {
    val outputDir = new File(libDir,"fn_data" )
     
 	def makeFnSimCmds(alleleFractions : Traversable[Double], depths:Traversable[String]):Traversable[CommandLineFunction] = {
 	  for {
 	      fraction <- alleleFractions 
 	      depth <- depths
 	  } yield (makeMixedBam(fraction, depth))

 	}
 	
 	private def makeMixedBam(alleleFraction: Double, depth: String): CommandLineFunction = {
        val tumorBams = getBams(depth) 
        val outBam = new File(outputDir, deriveBamName(alleleFraction, depth))
        val outIntervals  =  swapExt(outputDir, outBam, "bam", "interval_list")
        val spike = new SomaticSpike
        spike.alleleFraction = alleleFraction
        spike.outBam = outBam
        spike.tumorBams = tumorBams
        spike.outIntervals = outIntervals 
        spike.spikeSitesVCF = spikeSitesVCF
        spike
     }

    private def deriveBamName(alleleFraction: Double, depth : String) :String = {
      val bamNameTemplate = "NA12878_%s_NA12891_%s_spikein.bam"
      bamNameTemplate.format(depth, alleleFraction) 
    }
  } 


  def getBams(hexDigitString : String):List[File] = {
     try { 
      hexDigitString.map( digit => bamNameToFileMap( bamDigitToNameMap(digit)) ).toList
     } catch { 
         case e: Exception => 
             println(bamNameToFileMap)
             println(bamDigitToNameMap)
             throw e
     }   
  }
  
  def loadBamMap(bamMapFile : File):Map[Char, String] = {
        println("loading file")
        def splitLine(line: String): Option[(Char,String)] = {
          try{ val segments = line.split("\\s+")
           val char = segments(0).charAt(0)
           val name = segments(1).trim()
           Some(char, name)
          }catch { case e =>
           None
          } 
        }
        val fileLines = io.Source.fromFile(bamMapFile).getLines();
        val map : Map[Char, String] =  fileLines.map(splitLine).flatten.toMap
       
        map
  }
}


