#!/usr/bin/perl -w
use strict;

my $VCF_1KG = $ARGV[0];
my $NA12878_BAM = $ARGV[1];
my $NA12891_BAM = $ARGV[2];
my $REGIONS = $ARGV[3];
my $REF = $ARGV[4];
my $dir = $ARGV[5];

mkdir($dir);

my $GATK_BIN = "~/cga_home/gatk-protected/dist/GenomeAnalysisTK.jar";

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

