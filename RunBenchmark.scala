package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript

class RunBenchmark extends QScript {
  qscript =>

  val ALLELE_FRACTIONS = List(0.04, 0.1, 0.2, 0.4, 0.8)
  val TUMOR_DEPTHS = List("123456789ABC",
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
                          "1")

  val NORMAL_DEPTHS = List("D", "DE", "DEF", "DEFG", "DEFGHI", "FG", "HI")



  val FN_DIR = new File("fn_data")

  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  def script() {
    
  }



         //invokes <tool> with parameters <libDir><normal><tumor><reference><outputDir>
  class ToolInvocation extends  CommandLineFunction{
            @Input(doc="The tool to run")
            var tool: File = _

            @Input(doc="normal sample bam")
            var normal: File = _

            @Input(doc="tumor sample bam")
            var tumor: File = _

            @Input(doc="reference fasta")
            var reference: File = _

            @Output(doc="output directory")
            var outputDir: File =_
            
            def commandLine = required(tool)+
                              required(tool.getParent)+
                              required(normal)+
                              required(tumor)+
                              required(qscript.referenceFile)+
                              required(outputDir)
  }     

  object ToolInvocation {
    def apply(tool:File, normal:File, tumor:File, reference:File, outputDir:File) = {
        val ti = new ToolInvocation
        ti.tool = tool
        ti.normal = normal
        ti.tumor = tumor
        ti.reference = reference
        ti.outputDir = outputDir
        ti
    }
  } 

  def getNormalTumorPairs: List[(File, File)] = {
    def getPureFalsePositivePairs(normals: List[File], tumors: List[File]): List[(File, File)] ={
        for{
            normal <- normals
            tumor <- tumors
        } yield (normal, tumor)
    }

    def getSomaticSpikedPairs(normals: List[File], alleleFraction: List[Double], depths: List[String]):List[(File,File)] = {
        for{
            normal <- normals
            fraction <- alleleFraction
            depth <- depths
        } yield { 
            val tumorName = "NA12878_%s_NA12891_%s_spikein.bam".format(depth, fraction)
            val tumor = new File(FN_DIR, tumorName)
            (normal, tumor) 
        }
    }
      
    List.empty
  }


  def generateFalseNegativeCmds(toolsToTest: List[File], normal: File, tumors:List[File]):List[CommandLineFunction] = {
    def generateFalseNegativeCmd(tool :File, normal:File, tumor:File):CommandLineFunction ={
        val outputDir = new File("fn/"+tool.getName+"-N"+normal.getName+"-T"+tumor.getName)
        ToolInvocation(tool=tool, normal=normal, tumor=tumor, reference=referenceFile, outputDir=outputDir) 
    }
    
    for{
     tool <- toolsToTest
     tumor <- tumors   
    } yield generateFalseNegativeCmd(tool, normal, tumor)
  }
  
}

