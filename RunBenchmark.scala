package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript

class RunBenchmark extends QScript {
      qscript =>

      val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

      def script() {
      }
      


      def generateFalseNegativeCmds(toolsToTest: List[File], normal: File, tumors:List[File]):List[CommandLineFunction] = {
        def generateFalseNegativeCmd(tool :File, normal:File, tumor:File):CommandLineFunction ={
            val outputDir = tool.getName+"-N"+normal.getName+"-T"+tumor.getName 
                        
            new CommandLineFunction{
                def commandLine = required(tool)+
                                  required(tool.getParent)+
                                  required(normal)+
                                  required(tumor)+
                                  required(qscript.referenceFile)+
                                  required(outputDir)
            }
        }
        
        for{
         tool <- toolsToTest
         tumor <- tumors   
        } yield generateFalseNegativeCmd(tool, normal, tumor)
      }
}

