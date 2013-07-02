#!/usr/bin/perl -w
use strict;
my $dir = "vcf_data";
mkdir($dir);

my $VCF_1KG = "indels.vcf";
my $NA12878_BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam";
my $NA12891_BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam";
my $REGIONS = "chr20.interval_list";
my $REF = "/humgen/1kg/reference/human_g1k_v37_decoy.fasta";

my $GATK_BIN = "~/xchip/gatk-protected/dist/GenomeAnalysisTK.jar";

my $prefix = "chr1";
my $gtSites = "$dir/$prefix.trio.with.genotypes.at.1000G.sites.vcf";

my $cmd = "java -Xmx32g -jar $GATK_BIN " .
       "-T UnifiedGenotyper " .
       "-nt 16 " .
       "-R $REF "  .
       "-I $NA12878_BAM -I $NA12891_BAM " .
       "--genotyping_mode GENOTYPE_GIVEN_ALLELES " .
       "--alleles:VCF $VCF_1KG " .
       "-o $gtSites " .
       "-L $REGIONS";

unless(-e $gtSites) {
    print $cmd . "\n";
    system($cmd) == 0 or die();
}

my $out = "$dir/na12878_ref_NA12891_het_${prefix}_high_conf.vcf";
$cmd =    "java -jar $GATK_BIN " .
          "-T SelectVariants " .
          "-R $REF " .
          "-selectType INDEL " .
          "--variant:VCF $gtSites " .
          "-select 'vc.getGenotype(\"NA12878\").isHomRef() && vc.getGenotype(\"NA12891\").isHet() && vc.getGenotype(\"NA12891\").getPhredScaledQual() > 50 && QUAL > 50' " .
          "-L $REGIONS " .
          "-o $out";

unless(-e $out) {
    print $cmd . "\n";
    system($cmd) == 0 or die();
}

$out = "$dir/na12878_het_or_hom_nonref_${prefix}_high_conf.vcf";

$cmd = "java -jar $GATK_BIN " .
       "-T SelectVariants " .
       "-R $REF " .
       "--variant:VCF $gtSites " .
       "-selectType INDEL " .
       "-select '!vc.getGenotype(\"NA12878\").isHomRef() && vc.getGenotype(\"NA12878\").getPhredScaledQual() > 50 && QUAL > 50' " .
       "-L $REGIONS " .
       "-o $out";

unless(-e $out) {
    print $cmd . "\n";
    system($cmd) == 0 or die();
}

my $intervals = "$dir/na12878_het_or_hom_nonref_${prefix}_high_conf.intervals";
unless (-e $intervals) {
    `cat $out | grep -v "#" | awk '{ print \$1 ":" \$2 }' > $intervals`;
}

