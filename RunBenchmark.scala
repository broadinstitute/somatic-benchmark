package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.FileExtension

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


  val GERMLINE_NAME_TEMPLATE = "NA12878.somatic.simulation.merged.%s.bam"  
  val GERMLINE_MIX_DIR = new File("data_1g_wgs")

  def germlineMixFile(abrv :String) = AbrvFile.fromTemplate(GERMLINE_MIX_DIR, GERMLINE_NAME_TEMPLATE, abrv) 
  
  val FP_NORMAL_DEPTHS = List("D", "DE", "DEF", "DEFG", "DEFGHI", "FG", "HI").map(germlineMixFile(_) )
  val SPIKE_NORMAL_DEPTHS = List("DEFGHI").map(germlineMixFile(_) ) 
 

  val SPIKE_DIR = new File("fn_data")

  val referenceFile : File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  def script() {
    val tools = List(new AbrvFile("runIndelocator.sh", "Indelocator"))
    val cmds = getCommands(tools)
    cmds.foreach(add(_)) 
  }

  class AbrvFile(file: File , val abrv: String) extends File(file) with FileExtension {
    def withPath(path: String) = new AbrvFile(path, abrv) 
  }
 
  object AbrvFile{
    def fromTemplate(dir: File, nameTemplate: String, abrv: String) = {
        val name = nameTemplate.format(abrv)
        val file = new File(dir,name)
        new AbrvFile(file, abrv)
    }
  }

         //invokes <tool> with parameters <normal><tumor><reference><outputDir>
  class ToolInvocation extends  CommandLineFunction{
            @Input(doc="The script to run")
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

  def getCommands(tools: List[AbrvFile]) = {
    def getPureFalsePositivePairs(normals: List[AbrvFile], tumors: List[AbrvFile]) = {
        for{
            normal <- normals
            tumor <- tumors
        } yield (normal, tumor)
    }

    def getSomaticSpikedPairs(normals: List[AbrvFile], alleleFraction: List[Double], depths: List[String]) = {
        def tumorFile(fraction: Double, depth: String): AbrvFile = {
            val tumorName = "NA12878_%s_NA12891_%s_spikein.bam".format(depth, fraction)
            val tumorFile = new File(SPIKE_DIR, tumorName)
            val abrv = "%s_%s".format(depth,fraction)
            new AbrvFile(tumorFile, abrv)
        }
        
        for{
            normal <- normals
            fraction <- alleleFraction
            depth <- depths
        } yield { 
            val tumor = tumorFile(fraction, depth)          
            (normal, tumor) 
        }
    }
    
    val pureGermline = getPureFalsePositivePairs(FP_NORMAL_DEPTHS, TUMOR_DEPTHS.map(germlineMixFile(_)) )
    val spiked = getSomaticSpikedPairs(SPIKE_NORMAL_DEPTHS, ALLELE_FRACTIONS, TUMOR_DEPTHS)
    val pureGermlineCmds = generateCmds(tools, pureGermline, "germline_mix")
    val spikedCmds = generateCmds(tools, spiked, "spiked")  

    pureGermlineCmds ++ spikedCmds
  }


  def generateCmds(toolsToTest: List[AbrvFile], normalTumorPairs: List[(AbrvFile, AbrvFile)], outputDir: File):List[CommandLineFunction] = {
    def generateCmd(tool: AbrvFile, normal:AbrvFile, tumor: AbrvFile, outputDir: File): CommandLineFunction ={
        val individualOutputDir = new File(outputDir, "%s_N%S_T%S".format(tool.abrv, normal.abrv, tumor.abrv))
        ToolInvocation(tool=tool, normal=normal, tumor=tumor, reference=referenceFile, outputDir=individualOutputDir) 
    }
    
    for{
     tool <- toolsToTest
     (normal, tumor) <- normalTumorPairs   
    } yield generateCmd(tool, normal, tumor, outputDir)
  }
  
}

