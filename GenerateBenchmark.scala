package org.broadinstitute.sting.gatk.queue.qscripts.examples


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

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

  val prefix = "chr1"
  val validateIndelsFile = "indels.vcf"

  //TODOAn ugly hardcoded hack.  Must eventually be replaced when the number of divisions while fracturing is allowed to be changed from 6.
  val bamMapFile = new File("%s/louis_bam_1g_info.txt".format(libDir) )
  val bamDigitToNameMap :Map[Char,File]= loadBamMap(bamMapFile)
 
 




  def script()= {

    val mv = new make_vcfs 
    mv.indelFile = indelFile
    val mvOut = new File("dummyFile")
    mv.vcfOutFile = mvOut
    
    val sim = new make_fn_sim_dataset
    sim.inputFile = mvOut 

    val ss = new SomaticSpike
   
    add(mv)
    add( new create_1g_data)
    }

  class make_vcfs extends CommandLineFunction{
    @Input(doc="vcf file containing indels to use as true indel sites") var indelFile : File = _
    @Output(doc="dummy output for queue ordering") var vcfOutFile : File = _ 
    def commandLine = "%s/make_vcfs.pl %s".format(libDir, indelFile)
  }
  
  class make_fn_sim_dataset extends CommandLineFunction{
    @Input(doc="dummy input for queue ordering") var inputFile : File = _
    def commandLine = "%s/make_fn_sim_dataset.pl".format(libDir)
  }

  class create_1g_data extends CommandLineFunction{
    def commandLine = "%s/create_1g_sim_data.pl".format(libDir)
  }

  class SomaticSpike extends CommandLineGATK {
    this.analysis_type = "SomaticSpike"
    this.memoryLimit = 4 
  }

  class makeFalseNegatives(spikeSitesVCF : File, spikeInBam : File ,alleleFractions: Set[Double], depths : Set[String]) {
    val outputDir = new File("fn_data" )
     
 	def makeFalseNegativeDataSet(alleleFractions : Set[Double], depths:Set[String]):Set[SomaticSpike] = {
 	  for {
 	      fraction <- alleleFractions 
 	      depth <- depths
 	  } yield (makeMixedBam(fraction, depth))

 	}
 	
 	def makeMixedBam(alleleFraction: Double, depth: String): SomaticSpike = {
        val tumorBams = getBams(depth) 
        val spike = new SomaticSpike
        spike.reference_sequence = qscript.referenceFile
        spike.input_file ++= tumorBams
        spike.input_file :+= taggedFile( qscript.spikeContributorBAM, "spike")
        spike.intervals :+= qscript.spikeSitesVCF 
        spike.simulation_fraction = alleleFraction
        spike.minimum_qscore = 20 
        spike.out = outBAM
        spike.spiked_intervals_out = outIntervals
    }

  } 


  def getBams(hexDigitString : String):Set[File] = {
      hexDigitString.map(  bamDigitToNameMap ).toSet
  }
  def loadBamMap(bamMapFile : File):Map[Char, File] = {
      val fileLines = io.Source.fromFile(bamMapFile).getLines.toList;
      val map =  fileLines.map{ line:String => 
                                            val segments = line.split(" ")
                                            val char = segments(0).charAt(0)
                                            val file = new File (segments(1));
                                            (char, file)
                                          }.toMap
      map  
  }
}
        
