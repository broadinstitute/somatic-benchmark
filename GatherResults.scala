package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.FileExtension
/**
Traverse the output directories of RunBenchmark and gather results.
*/

class RunBenchmark extends QScript {
    qscript =>

    @Input(doc = "False positive test root directories.", required = false)
    var false_positive: Seq[File] = List(new File("germline") )

    @Input(doc = "False negative test root directories.", required = false)
    var false_negative: Seq[File] = List(new File("spiked") )

      
    def script() {
        val (fpResults, fnResults) = (false_positive, false_negative).map( searchForOutputFiles ) 

        
    
    }

 
    def searchForOutputFiles(roots: Seq[File]) = {
       finalVcfs = roots.flatMap( searchDirectory ) 
       finalVcfs
    }
    
    def searchDirectory(dir: File){
                   
    }
}

