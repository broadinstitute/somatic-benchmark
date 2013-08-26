#!/bin/sh -e

if [ $# != 4 ]
then
    echo "Usage runIndelocator.sh <normal bam> <tumor bam> <reference> <output dir>"
    echo "Requires a working installation of oncotator."
    echo "Please edit this file to set the gatk path."
    exit 1 
fi

#must match <normal><tumor><reference><outputDir>
NORMALBAM=$1
TUMORBAM=$2
REFERENCE=$3
OUTPUTDIR=$4
#GATKJAR=/xchip/cga2/louisb/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/modules/IndelGenotyper/GenomeAnalysisTK-53.5759.jar

#individual id
INDIVIDUAL="sample"
#directory to place outputs and intermediates


LIBDIR="/xchip/cga_home/louisb/indelocator"
GATKJAR=$GATK

#Bam files
#NORMALBAM=${TESTDIR}"/data/HCC1143_BL.cghub.ccle.small.bam"
#TUMORBAM=${TESTDIR}"/data/HCC1143.cghub.ccle.small.bam"

PANEL1="/cga/tcga-gsc/benchmark/Indels/NormalDB/normal_panel_indels_1KG_db.sorted.txt"
PANEL2="/cga/tcga-gsc/benchmark/Indels/NormalDB/normal_panel_indels_THCA_db.sorted.txt"


mkdir -p $OUTPUTDIR

#getblackList.sh   --access broad svn to load a lane blacklist
#echo "Getting lane blacklist."
#./LaneBlackList/getBlackList.sh 

#2: IndelGenotyper
echo "Running IndelGenotper"
java -Xmx3g -jar $GATKJAR \
	-T SomaticIndelDetector \
	-R $REFERENCE \
	-I:normal $NORMALBAM \
	-I:tumor $TUMORBAM \
	-verbose ${OUTPUTDIR}/${INDIVIDUAL}.indels.txt \
	-o ${OUTPUTDIR}/${INDIVIDUAL}.indels.vcf \
	--window_size 400 \
	--maxNumberOfReads 8000 \
	-filter 'T_COV<6||N_COV<4||(T_INDEL_F<=0.3&&T_CONS_CNT<7)||T_INDEL_CF<=0.7' \

#3.FilterIndelCalls -- N-FilterIndelCalls
echo "N-Filtering"
${LIBDIR}/FilterIndelCalls/filterSingleSampleCalls.pl \
	--calls ${OUTPUTDIR}/${INDIVIDUAL}.indels.txt \
	--prefix "N_" \
	--max_cons_av_mm 1000 \
	--max_ref_av_mm 1000 \
	--max_cons_nqs_av_mm 1000 \
	--min_ref_nqs_av_qual 0 \
	--min_cons_nqs_av_qual 0 \
	--min_cons_count 0 \
	--min_readpos_median 10 \
	--min_readpos_mad 3 \
	--mode ANNOTATE \
	--output ${OUTPUTDIR}/${INDIVIDUAL}.n_filtered.indels.txt 

#4 FilterIndelCalls -- T_FilterIndelCallsls
echo "T-Filtering"
${LIBDIR}/FilterIndelCalls/filterSingleSampleCalls.pl \
	--calls ${OUTPUTDIR}/${INDIVIDUAL}.n_filtered.indels.txt \
	--prefix "T_" \
	--max_cons_av_mm 4 \
	--max_ref_av_mm 4 \
	--max_cons_nqs_av_mm 0.3 \
	--min_ref_nqs_av_qual 15 \
	--min_cons_nqs_av_qual 15 \
	--min_cons_count 0 \
	--min_readpos_median 10 \
	--min_readpos_mad 3 \
	--mode ANNOTATE \
	--output ${OUTPUTDIR}/${INDIVIDUAL}.t_filtered.indels.txt 

#5. FilterIndelCallsByPanelDB
echo "Filtering by panel of normals"
${LIBDIR}/FilterIndelCallsByGermline/filterIndelCallsByPanelDB.pl \
	--library ${LIBDIR}/FilterIndelCallsByGermline \
	--ref $REFERENCE \
	--calls ${OUTPUTDIR}/${INDIVIDUAL}.t_filtered.indels.txt\
	--panel $PANEL1 \
	--panel $PANEL2 \
	--window 10 \
	--somatic \
	--tag "1KG,THCA" \
	--output $OUTPUTDIR/${INDIVIDUAL}.filtered.panel_marked.indels.txt \
	--filter '$EVT>=10'
echo "Done filtering"

#6 Convert t_filtered  txt file into vcf format
echo "Converting n_filtered indel file to vcf"
${LIBDIR}/indelsToVcf.sh ${LIBDIR} $OUTPUTDIR/${INDIVIDUAL}.n_filtered.indels.txt  ${OUTPUTDIR}/${INDIVIDUAL}.n_filtered.indels.vcf
echo "Done conversion."

#7 Convert t_filtered  txt file into vcf format
echo "Converting t_filtered indel file to vcf"
${LIBDIR}/indelsToVcf.sh ${LIBDIR} $OUTPUTDIR/${INDIVIDUAL}.t_filtered.indels.txt  ${OUTPUTDIR}/${INDIVIDUAL}.t_filtered.indels.vcf
echo "Done conversion."

#8 Convert final indel.txt file into vcf format
echo "Converting filtered indel file to vcf"
${LIBDIR}/indelsToVcf.sh ${LIBDIR} $OUTPUTDIR/${INDIVIDUAL}.filtered.panel_marked.indels.txt  ${OUTPUTDIR}/${INDIVIDUAL}.panel_filtered.indels.vcf
echo "Done conversion."

#9 Copy out final results
echo "Copying indel results final.indels.vcf"
cp ${OUTPUTDIR}/${INDIVIDUAL}.t_filtered.indels.vcf ${OUTPUTDIR}/final.indels.vcf
echo "Copy complete"

