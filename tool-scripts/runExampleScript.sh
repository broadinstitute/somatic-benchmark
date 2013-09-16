#!/bin/sh -e

if [ $# != 4 ]
then
    echo "Usage runExampleScript.sh <normal bam> <tumor bam> <reference> <output dir>"
    exit 1 
fi

#must match <normal><tumor><reference><outputDir>
NORMALBAM=$1
TUMORBAM=$2
REFERENCE=$3
OUTPUTDIR=$4

#Create output directory
mkdir -p $OUTPUTDIR

#Fill in the commands to run your caller here
echo "Pretending to run a caller to run a caller"
echo "normal=$NORMALBAM"
echo "tumor=$TUMORBAM"
echo "reference=$REFERENCE"
echo "Done pretending"

#The final output files must be placed in outputDir and called final.indels.vcf and final.snps.vcf 

touch ${OUTPUTDIR}/final.indels.vcf

