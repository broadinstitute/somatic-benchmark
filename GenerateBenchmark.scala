package org.broadinstitute.sting.gatk.queue.qscripts.examples


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class GenerateBenchmark extends QScript {
  qscript =>

  val indelFile : File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf")
  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bam1 : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam")
  val bam2 : File = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam")
  val libDir : File = new File("/xchip/cga2/louisb/indelocator-benchmark/sim-updated")

  var prefix = "chr1"
  var validateIndelsFile = "indels.vcf"


  def script() {

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
  }

  class makeFalseNegatives(spikeSitesVCF : File, spikeInBam : File ,alleleFractions: Set[Double], depths : Set[String])  {
    val outputDir = new File("fn_data" )
     
 	def makeFalseNegativeDataSet(alleleFractions : Set[Double], depths:Set[String]):List[SomaticSpike]{
 	  for {
 	      fraction <- alleleFractions 
 	      depth <- depths
 	  } yield (makeMixedBam(fraction, bamMap(depth)))

 	}
 	
 	def makeMixedBam(fraction: Double, tumorBams: List[File]): SomaticSpike{
        spike = new SomaticSpike
        spike.reference = referenceFile
    
    }
  } 
}
        
