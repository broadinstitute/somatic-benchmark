package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.FileExtension
import java.nio.file.{Files, Paths}
import java.io.FilenameFilter
import java.io
import org.broadinstitute.sting.queue.util.Logging

/**
Traverse the output directories of RunBenchmark and gather results.
*/

class GatherResults extends QScript with Logging{
    qscript =>

    @Input(doc = "False positive test root directories.", required = false)
    var false_positive: Seq[File] = List(new File("germline_mix") )

    @Input(doc = "False negative test root directories.", required = false)
    var false_negative: Seq[File] = List(new File("spiked") )

      
    def script() {
        val fpResults  = searchForOutputFiles( false_positive )
        val fnResults  = searchForOutputFiles( false_negative )
        logger.info("FP:"+ fpResults)
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
            files.find( _.getName() == "sample.final.indels.vcf" )
        } else {
            None
        }
    }
}

