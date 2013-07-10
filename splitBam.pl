#!/usr/bin/perl -w
use strict;
use IO::File;

if (scalar(@ARGV) <= 3) {
    die("usage: fracture_bam.pl <bam> <name sorted bam> <outfile1..> <outfile2..> <outfile_n...> \n");
}

my ($bam, $ns_bam) = @ARGV[0,1];
my @files = @ARGV[2..(scalar(@ARGV)-1)];
#my $outmask = "$outdir/NA12878.WGS.somatic.simulation.$library.%03d";

my $pieces = scalar(@files);
#initialize the random number generator
srand 28482;


    # initialize with the the header
    my @handles = ();
    for (my $i=1; $i<=$pieces; $i++) {
        #my $name = sprintf("$outmask.sam",$i);
        my $name = $files[$i-1] ;
        print "Initializing output file $name\n";
        system("samtools view -H $bam > $name") == 0 or die();
        my $handle = IO::File->new($name, 'a');
        push(@handles, $handle);
    }

    my $cmd = "samtools view $ns_bam ";
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
