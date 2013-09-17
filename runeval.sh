java -jar $GATK -T VariantEval -R ~/reference/human_g1k_v37.fasta \
 -eval:.8 fn_data/NA12878_123456789AB_NA12891_0.8_spikein.ug.vcf \
 -eval:.4 fn_data/NA12878_123456789AB_NA12891_0.4_spikein.ug.vcf \
 -eval:.2 fn_data/NA12878_123456789AB_NA12891_0.2_spikein.ug.vcf \
 -eval:.1 fn_data/NA12878_123456789AB_NA12891_0.1_spikein.ug.vcf \
 -eval:.04 fn_data/NA12878_123456789AB_NA12891_0.04_spikein.ug.vcf \
 -noST -ST AltReadFraction -o combo.gsareport

