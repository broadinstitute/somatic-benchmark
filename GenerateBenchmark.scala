package org.broadinstitute.org.cga

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.TaggedFile
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source
import java.io.{File, IOException}

import scala.collection.immutable.{HashMap, Map}

class GenerateBenchmark extends QScript {
  qscript =>
  //TODO implement these as cmdline paramaters isntead of hard coding them
  val indelFile : File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf")
  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bam : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
  val spikeContributorBAM : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")

  val libDir : File = new File("/xchip/cga2/louisb/indelocator-benchmark/sim-updated")
  
  val intervalFile = new File(libDir,"chr20.interval_list" )
  val spikeSitesVCF = new File(libDir,"vcf_data/na12878_ref_NA12891_het_chr1_high_conf.vcf" )

  val gatk : File = new File("/xchip/cga2/louisb/gatk-protected/dist/GenomeAnalysisTK.jar")
  val prefix = "chr1"
  val validateIndelsFile = "indels.vcf"

  val sortSamPath = "/seq/software/picard/current/bin/SortSam.jar"
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
    
    qscript.bamNameToFileMap = splitBams.map( bam => (fractureOutDir +"/" + bam.toString(), bam)).toMap 
    
    //use SomaticSpike to create false negative test data 
    val makeFnCommands = new FalseNegativeSim(spikeSitesVCF,spikeContributorBAM)
    val falseNegativeCmds = makeFnCommands.makeFnSimCmds( alleleFractions, depths)
    falseNegativeCmds.foreach(add(_))
  }

  object FractureBams {
      val FILE_NAME_PREFIX = "NA12878.WGS" 
      class NameSortBamByLibrary extends CommandLineFunction {
        @Input(doc="reference fasta file")
        var reference :File = qscript.referenceFile
         
        @Input(doc="bam file to sort")
        var inputBam : File = _
        
        @Input(doc="interval file")
        var interval : File = qscript.intervalFile 
        
        
        @Argument(doc="library name")
        var library :String = _

        @Output(doc="sorted bam of only the reads from one library")
        var nameSortedBam : File = _
       
        val printReadsCmd = required("java") +
                            required("-Xmx2g") +
                            required("-jar", gatk)+
                            required("-T","PrintReads")+
                            required("-l","Error")+
                            required("-log",nameSortedBam + "printreads.log")+
                            required("-rf","DuplicateRead")+
                            required("-rf","FailsVendorQualityCheck")+
                            required("-rf","UnmappedRead")+
                            required("-rf","LibraryRead")+
                            required("--library", library)+
                            required("-R",reference)+
                            required("-I",inputBam)+
                            required("-L", interval)

        val sortSamCmd = required("java")+
                         required("-Xmx16g")+
                         required("-jar", qscript.sortSamPath)+
                         required("VALIDATION_STRINGENCY=","SILENT",spaceSeparated=false)+
                         required("MAX_RECORDS_IN_RAM=","4000000", spaceSeparated=false)+
                         required("TMP_DIR=",qscript.tmpdir, spaceSeparated=false)+
                         required("I=","/dev/stdin",spaceSeparated=false)+
                         required("O=",nameSortedBam, spaceSeparated=false)+
                         required("SO=","queryname",spaceSeparated=false)+
                         required("COMPRESSION_LEVEL",1,spaceSeparated=false)
        def commandLine = printReadsCmd + required("|",escape=false) + sortSamCmd
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
          
          val splitBams: List[File] = splitSams.map(swapExt(_, "sam", "bam")) 
            
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

  class SomaticSpike extends CommandLineFunction with Logging{ 
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

  class create_1g_data extends CommandLineFunction{
    def commandLine = "%s/create_1g_sim_data.pl".format(libDir)
  }


  class FalseNegativeSim(spikeSitesVCF : File, spikeInBam : File) {
    val outputDir = new File(libDir,"fn_data" )
    val bamNameTemplate = outputDir+"/NA12878_%s_NA12891_%s_spikein.bam"
     
 	def makeFnSimCmds(alleleFractions : Traversable[Double], depths:Traversable[String]):Traversable[CommandLineFunction] = {
 	  for {
 	      fraction <- alleleFractions 
 	      depth <- depths
 	  } yield (makeMixedBam(fraction, depth))

 	}
 	
 	private def makeMixedBam(alleleFraction: Double, depth: String): CommandLineFunction = {
        val tumorBams = getBams(depth) 
        val outBam = deriveBamName(alleleFraction, depth)
        val outIntervals  = swapExt(outBam, "bam", "interval_list")
        val spike = new SomaticSpike
        spike.alleleFraction = alleleFraction
        spike.outBam = outBam
        spike.tumorBams = tumorBams
        spike.outIntervals = outIntervals 
        spike.spikeSitesVCF = spikeSitesVCF
        spike
     }

    private def deriveBamName(alleleFraction: Double, depth : String) :String = {
         bamNameTemplate.format(depth, alleleFraction) 
    }
  } 


  def getBams(hexDigitString : String):List[File] = {
      hexDigitString.map( digit => bamNameToFileMap( bamDigitToNameMap(digit)) ).toList
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


