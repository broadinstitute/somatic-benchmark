package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.FileExtension
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile

class RunBenchmark extends QScript {
  qscript =>

  @Argument(doc="If this is set than only 1 set of files will be run instead of the complete array.", required=false)
  var is_test = false


  @Argument(fullName="tool", shortName="t", doc="The name of a tool to run.  A matching script named run<Tool>.sh must be placed in the tool-scripts directory.")
  var tool_names: List[String]  = Nil

  @Argument(fullName="no_false_positives", shortName="nofp", doc="Run false positive analysis.", required=false)
  var no_false_positives: Boolean = false

  @Argument(fullName="no_false_negatives", shortName="nofn", doc="Run false negative analysis.", required=false)
  var no_false_negatives: Boolean = false

  lazy val ALLELE_FRACTIONS = if(is_test) List(0.8) else List(0.04, 0.1, 0.2, 0.4, 0.8)
  lazy val TUMOR_DEPTHS = if(is_test) List("123456789ABC") else List("123456789ABC",
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

  lazy val FP_NORMAL_DEPTHS = (if (is_test) List("DEFGHI") else List("D", "DE", "DEF", "DEFG", "DEFGH","DEFGHI", "FG", "HI")).map(germlineMixFile )
  val SPIKE_NORMAL_DEPTHS = List("DEFGHI").map(germlineMixFile )


  val SPIKE_DIR = new File("fn_data")

  val referenceFile : File = new File("/home/unix/louisb/cga_home/reference/human_g1k_v37_decoy.fasta")

  def script() {
    val tools = getTools(tool_names)

    if (!no_false_positives) {
        val (outFPDirs, fpCmds) = getFalsePositiveCommands(tools).unzip
        fpCmds.foreach(add(_))
    }

    if (!no_false_negatives) {
        val (outSpikedDirs, spikeCmds ) = getSpikedCommands(tools).unzip
        spikeCmds.foreach(add(_))
    }
  }

  def getTools(names: List[String]):List[AbrvFile] = {
      val toolDir = "tool-scripts"
      names.map{name =>
          val toolFile = new File(toolDir, "run%s.sh".format(name))
          if( ! toolFile.exists() ) throw new CouldNotReadInputFile(toolFile, " file does not exist.")
          new AbrvFile(toolFile, name)
      }
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
  class ToolInvocation extends  CommandLineFunction with RetryMemoryLimit{
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
        ti.memoryLimit = 4
        ti
    }
  }

  def getFalsePositiveCommands(tools: List[AbrvFile]) = {
    def getPureFalsePositivePairs(normals: List[AbrvFile], tumors: List[AbrvFile]) = {
        for{
            normal <- normals
            tumor <- tumors
        } yield (normal, tumor)
    }

    val pureGermline = getPureFalsePositivePairs(FP_NORMAL_DEPTHS, TUMOR_DEPTHS.map(germlineMixFile) )
    generateCmds(tools, pureGermline, "germline_mix")
 }

 def getSpikedCommands(tools:List[AbrvFile]) = {
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

    val spiked = getSomaticSpikedPairs(SPIKE_NORMAL_DEPTHS, ALLELE_FRACTIONS, TUMOR_DEPTHS)
    generateCmds(tools, spiked, "spiked")
 }


  def generateCmds(toolsToTest: List[AbrvFile], normalTumorPairs: List[(AbrvFile, AbrvFile)], outputDir: File):List[(File, CommandLineFunction)] = {
    def generateCmd(tool: AbrvFile, normal:AbrvFile, tumor: AbrvFile, outputDir: File): (File, CommandLineFunction) ={
        val individualOutputDir = new File(outputDir, "%s_N%S_T%S".format(tool.abrv, normal.abrv, tumor.abrv))
        (individualOutputDir, ToolInvocation(tool=tool, normal=normal, tumor=tumor, reference=referenceFile, outputDir=individualOutputDir) )
    }

    for{
     tool <- toolsToTest
     (normal, tumor) <- normalTumorPairs
    } yield generateCmd(tool, normal, tumor, outputDir)
  }

}

