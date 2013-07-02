from subprocess import call


bamBasePath= "/xchip/cga/home/kcibul/analysis/mutect/paper/sim_3.0/data_1g_wgs/"

def main():
    normals="""NA12878.somatic.simulation.merged.D.bam
NA12878.somatic.simulation.merged.DE.bam
NA12878.somatic.simulation.merged.DEF.bam
NA12878.somatic.simulation.merged.DEFG.bam
NA12878.somatic.simulation.merged.DEFGH.bam
NA12878.somatic.simulation.merged.DEFGHI.bam
NA12878.somatic.simulation.merged.FG.bam
NA12878.somatic.simulation.merged.HI.bam""".split("\n")

    tumors ="""NA12878.somatic.simulation.merged.1.bam
NA12878.somatic.simulation.merged.12.bam
NA12878.somatic.simulation.merged.123.bam
NA12878.somatic.simulation.merged.1234.bam
NA12878.somatic.simulation.merged.12345.bam
NA12878.somatic.simulation.merged.123456.bam
NA12878.somatic.simulation.merged.1234567.bam
NA12878.somatic.simulation.merged.12345678.bam
NA12878.somatic.simulation.merged.123456789.bam
NA12878.somatic.simulation.merged.123456789A.bam
NA12878.somatic.simulation.merged.123456789AB.bam
NA12878.somatic.simulation.merged.123456789ABC.bam""".split("\n")

    INDIV = "NA12878"
    RUN_STRELKA_PATH ="/home/unix/louisb/xchip/indelocator-benchmark/falsepositives/runStrelka.sh"
    STRELKA_LIBDIR="/home/unix/louisb/xchip/indelocator-benchmark/falsepositives"
    GATK_PATH="/xchip/cga2/louisb/gatk-protected/dist/GenomeAnalysisTK.jar"
    RUN_PIPELINE_PATH="/home/unix/louisb/xchip/indelocator-benchmark/falsepositives/runPipeline.sh"
    RUN_PIPELINE_MODULE_PATH="/xchip/cga2/louisb/indelocator"
    REFERENCE = "/xchip/cga/home/kcibul/analysis/mutect/paper/sim_3.0/refs/human_g1k_v37_decoy.fasta"
    RUN_FILTERS_PATH="/home/unix/louisb/xchip/indelocator-benchmark/falsepositives/runFilters.sh"

    normalBam = bamBasePath+"NA12878.somatic.simulation.merged.D.bam"
    tumorBam = bamBasePath+"NA12878.somatic.simulation.merged.12345689ABC.bam"

    bsub(*runFilters(RUN_FILTERS_PATH,normalBam,tumorBam,RUN_PIPELINE_MODULE_PATH))
    normalBam = bamBasePath+"NA12878.somatic.simulation.merged.DEFGH.bam"
    tumorBam = bamBasePath+"NA12878.somatic.simulation.merged.123.bam"

    bsub(*runFilters(RUN_FILTERS_PATH,normalBam,tumorBam,RUN_PIPELINE_MODULE_PATH))

   # normals = ["NA12878.somatic.simulation.merged.D.bam"]
   # tumors = ["NA12878.somatic.simulation.merged.1.bam"]
    #for normal in normals:
    #    for tumor in tumors:
    #        tumorBam = bamBasePath+tumor
    #        normalBam = bamBasePath+normal
    #        print getName(normalBam, tumorBam)
    #        bsub(*runFilters(RUN_FILTERS_PATH,normalBam,tumorBam,RUN_PIPELINE_MODULE_PATH))
   #         bsub(*runPipeline(GATK_PATH, RUN_PIPELINE_PATH, normalBam, tumorBam, RUN_PIPELINE_MODULE_PATH))
   #         bsub(*runStrelka(RUN_STRELKA_PATH,STRELKA_LIBDIR,normalBam, tumorBam, REFERENCE, INDIV))

def runFilters(runFiltersPath,normalBam,tumorBam,modulePath):
    jobName="runFilters"
    cmd = [runFiltersPath,
            modulePath ,
            "indl_{name}".format(name=getName(normalBam,tumorBam)) ]
    return (cmd, jobName )

def runStrelka(runStrelkaPath, libdir, normalBam, tumorBam, reference, indiv):
    jobName="strelka_"+getName(normalBam,tumorBam)
    cmd = [runStrelkaPath,
           libdir,
           normalBam,
           tumorBam,
           reference,
           indiv,
       jobName]
    return (cmd, jobName, "-q week")

def getName(normalBam, tumorBam):
    return "N{norm}_T{tum}".format(norm=normalBam.split(".")[-2], tum=tumorBam.split(".")[-2])


def runPipeline(gatkPath, runPipelinePath, normalBam, tumorBam, modulePath ):
    jobName = "indl_"+getName(normalBam,tumorBam)
    cmd = [ runPipelinePath,
            gatkPath,
            normalBam,
            tumorBam,
            modulePath,
            "indl_{name}".format(name=getName(normalBam,tumorBam)) ]
    return (cmd, jobName, "-q week", "-R rusage[mem=3]")


def bsub(cmd, jobName, *args): 
    print "Running:\n {cmd} \n as jobName={jobName}".format(cmd=cmd,jobName=jobName)
    call(["bsub",
        #"-N",
        "-P {job}".format(job=jobName),
        "-o","stdout-{name}-%J.txt".format(name=jobName)]
        + list(args)    
        + cmd )

if __name__ == "__main__":
    main()
