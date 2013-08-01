package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.FileExtension
import java.nio.file.{Files, Paths}
import java.io.{BufferedWriter, FileWriter, FilenameFilter}
import java.io
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.commandline

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


    def analyzePositives(files: Seq[File]) = ()
    def analyzeNegatives(files: Seq[File]) = {
        val (jobs, outputs) = files.map{ file =>
            val (tool, normalName, tumorName, fraction) = parseFilePath(file)
            val vcfdiff = new VcfDiff
            vcfdiff.prefix = file.getParent + "/vcfout"
            vcfdiff.diff = getTumorFromName(tumorName, fraction)
            vcfdiff.vcf = file
            val diffOut = new File( file.getParent, "vcfout.diff.sites_in_files")
            vcfdiff.diffOut = diffOut
            (vcfdiff, diffOut)
        }.unzip

        jobs.foreach(add(_))

        val stats = new ComputeIndelStats{
            this.sitesFiles = outputs.toList
        }

    }

    def getTumorFromName(name: String, fraction: Double)={
        new File("fn_data","NA12878_%s_NA12891_%s".format(name, fraction))
    }

    class ComputeIndelStats extends InProcessFunction(){
        @Input(doc="all files from vcfdiff")
        var sitesFiles: List[File] = _

        @Output(doc="file to print results in")
        var results: File = _

        def run() {
            val indels = sitesFiles.map(extractResultsFromVcftoolsDiff(_))
            val fw = new FileWriter(results.getAbsoluteFile());
            val bw = new BufferedWriter(fw);

            indels.zip(sitesFiles).foreach{ pair =>
                val ((first,second), file ) = pair
                val onlyFirst = first.diff(second).size
                val onlySecond =
                bw.write("")
            }

            bw.close();
        }
    }

    class VcfDiff extends CommandLineFunction {
        @Input(doc="vcf file to be compared")
        var vcf: File = _

        @Input(doc="vcf file to compare against")
        var diff: File =_

        @Argument(doc="prefix string for output files")
        var prefix: String =_

        @Output(doc="output file")
        var diffOut: File = _

        def commandLine: String = required("vcftools") +
                                  required("--vcf", vcf) +
                                  required("--diff", diff) +
                                  required("--out", prefix)
    }

    def parseFilePath(file: File){
        val parent = file.getParentFile.getName
        val (tool:String, normalWithN:String, tumorWithT:String, fraction: Double) = parent.split('_')
        (tool, normalWithN.drop(1), tumorWithT.drop(1), fraction)
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
            files.find( _.getName() == "final.indels.vcf" )
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

        val indels: List[(Option[String], Option[String])] = Source.fromFile(vcftoolsOut).getLines().toList.map(assignIndels)

        val (first, second) = indels.unzip
        (first.flatten, second.flatten)

    }
}
