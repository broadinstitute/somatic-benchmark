#!/usr/bin/perl -w
use strict;
use IO::File;

if (scalar(@ARGV) < 6) {
	die("usage: <bam> <ref> <interval-list> <pieces> <outdir> <library1> <library2> ... <library n>\n");
}

	
my $bam = @ARGV[0]
my $ref = @ARGV[1]
my $interval = @ARGV[2]
my $pieces = @ARGV[3]
my $outdir = @ARGV[4]
my @libraries = @ARGV[5..scalar(@ARGV)]

foreach my $library (@libraries) {
    my $cmd = "./fracture_bam.pl $bam $ref $library $interval $pieces $outdir"; 
    print $cmd . "\n";
#    system($cmd) == 0 or die();
}
