Tools for comparing somatic variant calling tools.

7/25/2013:  This is a work in progress.  Tools are incomplete, untested, undocumented, and subject to change.

This is the new home for the SomaticSpike walker.  

INSTALLATION
Assuming you have perl and python 2 installed already and are using bash as your shell.
Also requires java 7.

To run with indelocator(assuming access to the broad network):
1. clone necessary repos

git clone git@github.com:broadinstitute/somatic-benchmark
git clone git@github.com:broadinstitute/indelocator
git clone git@github.com:broadgsa/gatk-protected

2. build gatk-protected including SomaticSpike:

cd gatk-protected
ant -Dexternal.dir=../somatic-benchmark

3. build gatk-protected again with indelocator
ant -Dexternal.dir=../indelocator

4. be sure that you built Queue
ant queue

5. move to the somatic-benchmark directory
cd ../somatic-benchmark

6. edit make_vcfs.pl and set GATK_BIN to your gatk jar file which is gatk-protected/dist/GenomeAnalysisTK.jar

7. edit tool-scripts/runIndelocator.sh
    -set LIBDIR to point to your indelocator directory
    -set GATKJAR to point to your gatk jar (or set the GATK environmental variable)

8. install this particular unstable build of oncotator
    a. deactivate any existing python virtual environment (will do nothing if you don't have one active)
        deactivate

    b. activate a new virtual environment for oncotator
        source /xchip/cga1/lichtens/test_oncotator/test_env/bin/activate
        ( it helps to add this to your .my.bashrc so you don't forget) 
    c. check that you have the right version
        which oncotator #should display /xchip/cga1/lichtens/test_oncotator/test_env/bin/oncotator

    (d. to return to your regular python environment later, use deactivate)

9. run GenerateBenchmark.scala 
    java -jar gatk-protected/dist/Queue.jar -S GenerateBenchmark.scala -run -bsub
    
10. run RunBenchmark.scala with indelocator 
    java -jar gatk-protected/dist/Queue.jar -S RunBenchmark.scala -t Indelocator -run -bsub

11. analysis module is incomplete...

