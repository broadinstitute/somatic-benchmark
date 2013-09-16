This directory is where the benchmarker looks for wrappers to invoke when given a --tool parameter.


Instructions for adding another tool to the benchmarking suite

-Copy runExampleScript.sh and name it runYourTool.sh
-Complete the wrapper script so that it runs your tool on the given inputs.
--Ensure that the variant call output files are vcfs labeled final.indels.vcf or final.snps.vcf

Now you can run it on the generated bams by using RunBenchmark.scala --tool YourTool

