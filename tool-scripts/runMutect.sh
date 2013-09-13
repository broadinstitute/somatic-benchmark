#!/bin/sh -e

if [ $# != 4 ]
then
    echo "Usage runMutect.sh <normal bam> <tumor bam> <reference> <output dir>"
    echo "Requires a working installation of oncotator."
    echo "Please edit this file to set the gatk path."
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

MutectJar=$GATK
IntervalFile='benchmark.interval_list'

java -jar $MutectJar \
'--analysis_type' 'MuTect' \
'--normal_sample_name' "NORMAL" \
'-I:normal' $NORMALBAM \
'--tumor_sample_name' "TUMOR" \
'-I:tumor' $TUMORBAM \
'--reference_sequence' $REFERENCE \
'--dbsnp' '/xchip/cga/reference/hg19/dbsnp_134_b37.leftAligned.vcf' \
'--cosmic' '/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf' \
'--normal_panel' '/xchip/cga/reference/hg19/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf' \
'--out' "${OUTPUTDIR}/call_stats.txt" \
'--vcf' "${OUTPUTDIR}/final.snps.vcf" \
'--only_passing_calls' \
'--coverage_file' "${OUTPUTDIR}/coverage.wig.txt" \
'--power_file' "${OUTPUTDIR}/power.wig.txt" \
'--downsample_to_coverage' '1000' \
'--enable_extended_output' \
'--intervals' $IntervalFile



