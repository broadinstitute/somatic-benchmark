#!/bin/sh -e

if [ $# != 4 ]
then
    echo "Usage runVarScan2.sh <normal bam> <tumor bam> <reference> <output dir>"
    exit 1 
fi

#must match <normal><tumor><reference><outputDir>
NORMALBAM=$1
TUMORBAM=$2
REFERENCE=$3
OUTPUTDIR=$4

VarJar=~/cga_home/tools/VarScan2/VarScan.v2.3.6.jar

#Create output directory
mkdir -p $OUTPUTDIR

#Fill in the commands to run your caller here
normal_pileup=" samtools mpileup -q 1 -f $REFERENCE $NORMALBAM"
tumor_pileup=" samtools mpileup -q 1 -f $REFERENCE $TUMORBAM"

set +o posix
java -jar $VarJar somatic \
<($normal_pileup) <($tumor_pileup) ${OUTPUTDIR}/output \
--output-vcf \
--strand-filter

java -jar $VarJar processSomatic ${OUTPUTDIR}/output.snp.vcf
java -jar $VarJar processSomatic ${OUTPUTDIR}/output.indel.vcf

#The final output files must be placed in outputDir and called final.indels.vcf and final.snps.vcf 

cp $OUTPUTDIR/output.snp.Somatic.hc.vcf final.snps.vcf
cp $OUTPUTDIR/output.indel.Somatic.hc.vcf final.indels.vcf
