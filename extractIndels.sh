GATK=~/xchip/gatk-protected/dist/GenomeAnalysisTK.jar

java -Xmx2g -jar $GATK \
 -R "/humgen/1kg/reference/human_g1k_v37.fasta" \
 -T SelectVariants \
 --variant $1 \
 -o "indels.vcf" \
 -L "chr20.interval_list" \
 -selectType INDEL
