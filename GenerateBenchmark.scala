package org.broadinstitute.org.cga

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.TaggedFile
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source
import java.io.{FileReader,IOException}

import scala.collection.immutable.{HashMap, Map}

class GenerateBenchmark extends QScript {
  qscript =>

  val indelFile : File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf")
  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bam1 : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
  val spikeContributorBAM : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")

  val intervalFile = new File("%s/chr20.interval_list".format(libDir) )
  val spikeSitesVCF = new File("%s/vcf_data/na12878_ref_NA12891_het_chr1_high_conf.vcf".format(libDir) )


  val libDir : File = new File("/xchip/cga2/louisb/indelocator-benchmark/sim-updated")
  val gatkDir : File = new File("/xchip/cga2/louisb/gatk-protected/dist/GenomeAnalysisTK.jar")
  val prefix = "chr1"
  val validateIndelsFile = "indels.vcf"

  //TODOAn ugly hardcoded hack.  Must eventually be replaced when the number of divisions while fracturing is allowed to be changed from 6.
  val bamMapFile = new File("%s/louis_bam_1g_info.txt".format(libDir) )
  var bamDigitToNameMap :Map[Char,File]= null 
 
  val alleleFractions = Set( 0.04, .1 , .2, .4, .8)

  val maxdepth = "123456789ABC"
  val depths = for(i <- 1 to maxdepth.length) yield maxdepth.substring(0, i) 
 

  def script()= {
    qscript.bamDigitToNameMap = loadBamMap(bamMapFile)
   
   
    val makeVcfs = new MakeVcfs
    val vcfOut = new File(libDir)  
    makeVcfs.indelFile = qscript.indelFile  
    makeVcfs.vcfOutFile = vcfOut
      
    add(makeVcfs)

    

    val makeFnCommands = new FalseNegativeSim(spikeSitesVCF,spikeContributorBAM)
    val cmds = makeFnCommands.makeFnSimCmds( alleleFractions, depths)
  
    cmds.foreach(cmd => add(cmd)) 
  }
  
  
  class SomaticSpike extends CommandLineFunction with Logging{ 

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
                          required("-jar", qscript.gatkDir)+
                          required("-T SomaticSpike")+
                          required("--reference_sequence", qscript.referenceFile)+
                          repeat("-I", tumorBams)+
                          required("-I:spike", qscript.spikeContributorBAM)+
                          required("--intervals",qscript.spikeSitesVCF)+
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


  class  FalseNegativeSim(spikeSitesVCF : File, spikeInBam : File) {
    val outputDir = new File("fn_data" )
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
        spike
     }

    private def deriveBamName(alleleFraction: Double, depth : String) :String = {
         bamNameTemplate.format(depth, alleleFraction) 
    }
  } 


  def getBams(hexDigitString : String):List[File] = {
      hexDigitString.map(  bamDigitToNameMap ).toList
  }
  
  def loadBamMap(bamMapFile : File):Map[Char, File] = {
        println("loading file")
        def splitLine(line: String): Option[(Char, File)] = {
          try{ val segments = line.split("\\s+")
           val char = segments(0).charAt(0)
           val file = new File (segments(1));
           Some(char, file)
          }catch { case e =>
           None
          } 
        }
        val fileLines = io.Source.fromFile(bamMapFile).getLines();
        val map : Map[Char, File] =  fileLines.map(splitLine).flatten.toMap
       
        map
  }
}


