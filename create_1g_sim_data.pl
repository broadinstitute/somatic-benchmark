#!/usr/bin/perl -w
use strict;
use IO::File;

my $bam = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.jaffe.b37_decoy.NA12878.bam";
my $ref = "/humgen/1kg/reference/human_g1k_v37_decoy.fasta";
my @libraries = ("Solexa-18483","Solexa-18484","Solexa-23661");
my $interval = "chr20.interval_list";
my $pieces = 6;
my $outdir = "data_1g_wgs";

my $QUEUE = "week";

unless(-e $outdir or mkdir $outdir){
    die "Can't create $outdir\n";
}

foreach my $library (@libraries) {
    my $cmd = "bsub -q week -P benchmark -o stdout-%J.txt -R rusage[mem=20] ./fracture_bam.pl $bam $ref $library $interval $pieces $outdir"; 
    print $cmd . "\n";
    system($cmd) == 0 or die();
}
