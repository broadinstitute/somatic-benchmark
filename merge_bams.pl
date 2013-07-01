#!/usr/bin/perl -w
use strict;
require "./SimulationCommon.pl";

my ($out, $bamMap) = ("data_1g_wgs", "bam_1g_info.txt");
my $picardBin="/seq/software/picard/current/bin"

#define the sets to run tumor -> normal, in hex
my @sets = ( 
             "123456789ABC",
             "123456789AB",
             "123456789A",
             "123456789",
             "12345678",
             "1234567",
             "123456",
             "12345",
             "1234",
             "123",
             "12",
             "1",
             "DEFGHI",
             "DEFGH",
             "DEFG", 
             "DEF",
             "DE",
             "FG",
             "HI",
             "D"
           );

for my $name (@sets) {

        my $out = sprintf("$out/NA12878.somatic.simulation.merged.%s.bam",$name);

        my @inputBams = getBams($name, $bamMap);
        my $outBai = $out;
        $outBai =~ s/bam/bai/g;

        my $cmd = "bsub -P benchmark -q week -o $out.lsfout \" java " .
                "-Xmx2g -jar $picardBin/MergeSamFiles.jar ".
                "CREATE_INDEX=True TMP_DIR=/broad/hptmp/louisb USE_THREADING=True " .
                "O=$out " .
                "I=" . join(" I=", @inputBams) . " && cp $outBai $out.bai \" "; 

        unless (-e "$out") {
                print $cmd . "\n";
                system($cmd) == 0 or die();
        }

}
