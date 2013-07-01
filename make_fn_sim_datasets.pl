#!/usr/bin/perl -w
use strict;
require "./SimulationCommon.pl";

my $outputDir = "fn_data";
my $SOMATIC_SPIKE = "~/xchip/gatk-protected/dist/GenomeAnalysisTK.jar";
my $REF = "/humgen/1kg/reference/human_g1k_v37.fasta";

my $NA12891_BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.bam";
my $SPIKE_SITES_VCF = "vcf_data/na12878_ref_NA12891_het_chr20_high_conf.vcf";
my $QUEUE = "hour";
my $bamMap = "bam_1g_info.txt";

my @mix_af = (0.04, 0.1, 0.2, 0.4, 0.8);
my @depths = ("123456789ABC"); #,
          #   "123456789AB",
          #   "123456789A",
          #   "123456789",
          #   "12345678",
          #   "1234567",
          #   "123456",
          #   "12345",
          #   "1234",
          #   "123",
          #   "12",
          #   "1");

mkdir($outputDir);
foreach my $mf (@mix_af) {
    foreach my $tumorName (@depths) {

        my $bam = sprintf("$outputDir/NA12878_%s_NA12891_%s_spikein.bam", $tumorName, $mf);
        my $int_out = $bam;
        $int_out =~ s/bam/intervals/g;

        my @tumorBams = getBams($tumorName, $bamMap);


      my $cmd = "bsub -R \"rusage[mem=4096]\" -o $bam.lsfout \" java " .
                 "-Xmx4g -jar $SOMATIC_SPIKE -T SomaticSpike " .
                 "-R $REF " .
                 "-I " . join(" -I ", @tumorBams) . " " . 
                 "-I:spike $NA12891_BAM " .
                 "-L $SPIKE_SITES_VCF " .
                 "--simulation_fraction $mf " .
                 "--minimum_qscore 20 " .
                 "-o $bam " .
                 "--spiked_intervals_out $int_out \"";

        unless (-e $bam) {
            print "Generating $bam\n";
            print $cmd . "\n";
            system($cmd) == 0 or die();
        }


    }
}
