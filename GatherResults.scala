package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import java.io.{BufferedWriter, FileWriter}
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source

/**
Traverse the output directories of RunBenchmark and gather results.

It assumes that each output is in it's own directory and the final file is called sample.final.indels.vcf"

*/

class GatherResults extends QScript with Logging{
    qscript =>

    @Input(doc = "False positive test root directories.", required = false)
    var false_positive: Seq[File] = List(new File("germline_mix") )

    @Input(doc = "False negative test root directories.", required = false)
    var false_negative: Seq[File] = List(new File("spiked") )


    def analyzePositives(files: Seq[File]) = {
    }


    def analyzeNegatives(files: Seq[File]) = {
        val (jobs, diffOuts) = files.map{ file =>
            val metaData = new DirectoryMetaData(file.getParentFile)
            val vcfdiff = new VcfDiff
            vcfdiff.outputPrefix = file.getParent + "/vcfout"
            vcfdiff.comparisonVcf = swapExt(metaData.tumor.getParentFile, metaData.tumor, "bam","vcf")
            vcfdiff.vcf = file
            val diffOut = new File( file.getParent, "vcfout.diff.sites_in_files")
            vcfdiff.differenceFile = diffOut
            (vcfdiff, diffOut)
        }.unzip

        jobs.foreach(add(_))

        val stats = new ComputeIndelStats{
            this.sitesFiles = diffOuts.toList
            this.results = new File("diffResults.txt")
        }

        add(stats)

    }


    /**
     * Container for directory level meta data.
     */
    class DirectoryMetaData(val file: File) {
        private def tumorFileFromName(name: String, fraction: Double)={
            new File("fn_data","NA12878_%s_NA12891_%s_spikein.bam".format(name, fraction))
        }

        private def normalFileFromName(name: String)={
            new File("data_1g_wgs", "NA12878.somatic.simulation.merged.%s.bam".format(name))
        }

        def hasSpikeIn: Boolean =( splits.length==4)

        //expecting filename in the format of Indelocator_NDEFGHI_T1234_0.4
        private val splits = file.getName.split("_")
        val tool: String = splits(0)

        val normalName: String = splits(1).drop(1)
        val normal: File = normalFileFromName(normalName)

        val tumorName :String = splits(2).drop(1)
        val (tumor, fraction) = if (hasSpikeIn) {
            val fraction = splits(3).toDouble
            ( tumorFileFromName(tumorName, fraction), Some(fraction) )
        } else {
            (normalFileFromName(tumorName), None)
        }

    }



    class ComputeIndelStats extends InProcessFunction(){
        @Input(doc="all files from vcfdiff")
        var sitesFiles: List[File] = _

        @Output(doc="file to print results in")
        var results: File = _

        def run() {
            val indels = sitesFiles.map(extractResultsFromVcftoolsDiff)
            val fw = new FileWriter(results.getAbsoluteFile)
            val bw = new BufferedWriter(fw)

            bw.write("Tool\tNormal\tTumor\tFraction\tFP\tFN\tMatched\n")

            indels.zip(sitesFiles).foreach{ pair =>
                val ((first,second), file ) = pair
                val onlyFirst = first.diff(second).size
                val onlySecond = second.diff(first).size
                val matches = first.intersect(second).size
                val metaData = new DirectoryMetaData(file.getParentFile)


                bw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n".format(metaData.tool, metaData.normalName, metaData.tumorName,
                                                           metaData.fraction, onlyFirst, onlySecond,matches))
            }

            bw.close()
        }
    }

    class VcfDiff extends CommandLineFunction {
        @Input(doc="vcf file to be compared")
        var vcf: File = _

        @Input(doc="vcf file to compare against")
        var comparisonVcf: File =_

        @Argument(doc="prefix string for output files")
        var outputPrefix: String =_

        @Output(doc="output file")
        var differenceFile: File = _

        def commandLine: String = required("vcftools") +
                                  required("--vcf", vcf) +
                                  required("--diff", comparisonVcf) +
                                  required("--out", outputPrefix)
    }



    

    def script() {
        val fpResults  = searchForOutputFiles( false_positive )
        val fnResults  = searchForOutputFiles( false_negative )
        logger.debug("Fp results:"+ fpResults)
        logger.debug("Fn results:"+ fnResults)

        analyzePositives(fpResults)
        analyzeNegatives(fnResults)

    }

 
    def searchForOutputFiles(roots: Seq[File]) = {
       val finalVcfs = roots.flatMap( searchRootDirectory )
       finalVcfs
    }
    
    def searchRootDirectory(dir: File) = {
        val resultDirs = dir.listFiles()
        val results = resultDirs.flatMap(checkForResultFile)
        results
    }

    def checkForResultFile(dir: File):Option[File] = {
        val files = dir.listFiles()
        if (files != null) {
            files.find( _.getName == "final.indels.vcf" )
        } else {
            None
        }
    }


    def extractResultsFromVcftoolsDiff(vcftoolsOut : File): (List[String], List[String]) = {
        def assignIndels(line: String):(Option[String], Option[String]) = {
                val tokens = line.split('\t')
                val indel = tokens(1)

                tokens(2) match{
                    case "1" => (Some(indel), None)
                    case "2" => (None, Some(indel))
                    case "B" => (Some(indel), Some(indel))
                    case _ => (None, None)
                }
        }

        val indels: List[(Option[String], Option[String])] = try {
            Source.fromFile(vcftoolsOut).getLines().toList.map(assignIndels)
        } catch {
            case e =>
            logger.error(e.getMessage)
            List((None, None))
        }

        val (first, second) = indels.unzip
        (first.flatten, second.flatten)

    }
}
