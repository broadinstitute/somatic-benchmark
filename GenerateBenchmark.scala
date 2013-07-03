package org.broadinstitute.sting.gatk.queue.qscripts.examples

import org.broadinstitute.sting.gatk.CommandLineGATK
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper
import org.broadinstitute.sting.gatk.walkers.variantutils.SelectVariants
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
    add(mv)
  }

  class make_vcfs extends CommandLineFunction{
    @Input(doc="vcf file containing indels to use as true indel sites") var indelFile : File = _
    def commandLine = "%s/make_vcfs %s".format(libDir, indelFile)
  }


}
        
