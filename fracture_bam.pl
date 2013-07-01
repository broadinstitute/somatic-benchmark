#!/usr/bin/perl -w
use strict;
use IO::File;


if (scalar(@ARGV) != 6) {
    die("usage: fracture_bam.pl <bam> <ref> <library> <interval> <pieces> <out-dir>\n");
}

my ($bam, $ref, $library, $interval, $pieces, $outdir) = @ARGV;

my $ns_bam = "$outdir/NA12878.WGS.original.regional.namesorted.$library.bam";
my $QUEUE = "week";
my $outmask = "$outdir/NA12878.WGS.somatic.simulation.$library.%03d";

my $PICARD_SORT_SAM_BIN = "/seq/software/picard/current/bin/SortSam.jar";
my $SAMTOOLS = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/samtools/samtools_0.1.16/bin/samtools";


my $tmpdir = "/broad/hptmp/louisb/sim";
my $gatk = "/xchip/cga2/louisb/gatk-protected/dist/GenomeAnalysisTK.jar";

#initialize the random number generator
srand 28482;

#
# Step 1: create regional name-sorted BAMs by library
# 
unless (-e $ns_bam) {
    my $library_rf = "";
    if ($library ne "all") { $library_rf = " -rf LibraryRead --library $library "; }

    print "Namesorting BAM...\n";

    my $print_reads_cmd = "java -Xmx2g -jar $gatk -T PrintReads -l ERROR -log $ns_bam.printreads.log --simplifyBAM -rf DuplicateRead -rf FailsVendorQualityCheck -rf UnmappedRead $library_rf -R $ref -I $bam -L $interval ";

    my $sort_sam_cmd = "java -Xmx16g -jar $PICARD_SORT_SAM_BIN VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$tmpdir I=/dev/stdin O=$ns_bam SO=queryname COMPRESSION_LEVEL=1";
    my $cmd = "$print_reads_cmd | $sort_sam_cmd";
    print "$cmd\n";
    system($cmd) == 0 or die();
}

#check if the last file exist
unless (-e sprintf("$outmask.sam",$pieces)) {

    # initialize with the the header
    my @handles = ();
    for (my $i=1; $i<=$pieces; $i++) {
        my $name = sprintf("$outmask.sam",$i);
        print "Initializing output file $name\n";
        system("$SAMTOOLS view -H $bam > $name") == 0 or die();
        my $handle = IO::File->new($name, 'a');
        push(@handles, $handle);
    }

    my $cmd = "$SAMTOOLS view $ns_bam ";
    print "Running $cmd\n";
    open(my $fh, "$cmd | ") or die $!;
    my $lastpair = "";
    my $lastpair_handle_id = -1;
 
    my $read_count = 0;
    while (my $line = <$fh> ) {
        my ($readname) = split("\t", $line);    

        # if it's part of a pair, write to the same filehandle
        if ($readname eq $lastpair) {
            $handles[$lastpair_handle_id]->print($line);
        } else {
            # pick a random number between 1 <= x <= $PIECES and write to that handle
            my $r = int(rand($pieces)); 

            $lastpair = $readname;
            $lastpair_handle_id = $r;
            $handles[$lastpair_handle_id]->print($line);
        }

        if ($read_count++ % 1000000 == 1) { print "Processed $read_count reads...\n"; }
    }
    close($fh);

    foreach my $handle (@handles) {
        $handle->close();
    }
}

# now convert them to BAM
for (my $i=1; $i<=$pieces; $i++) {
    my $name = sprintf("$outmask.sam",$i);
    my $out = $name;
    $out =~ s/sam/bam/g;

    unless(-e $out) {
        print "Converting $name to BAM\n";

        my $cmd = "bsub -P benchmark -q $QUEUE -o $out.lsfout java -Xmx2g -jar $PICARD_SORT_SAM_BIN TMP_DIR=$tmpdir I=$name O=$out SO=coordinate CREATE_INDEX=true QUIET=True ";
        print "$cmd\n";
        system($cmd) == 0 or die();
    }
}
