sub getBams {
    my ($name, $bamMap) = @_;
    my %bamInfo = ();
    foreach my $line (`cut -f1-2 $bamMap`) {
        chomp($line);
        my ($part, $file) = split("\t", $line);
        $bamInfo{$part} = $file;
    }

    my @bams = ();

    for( my $i = 0; $i<length($name); $i++) {
        my $s = substr($name,$i,1);
        push(@bams, $bamInfo{$s});
    }
    return @bams;
}

1; # need to end with a true value
