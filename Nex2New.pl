#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
    my $usage = "script.pl INFILE OUTFILE\n";
    my $infile = shift or die $usage;
    my $outfile = shift or die $usage;

my ($filein,$fileout) = @ARGV;
my ($format,$oformat) = qw(nexus newick);
my $in = Bio::TreeIO->new(-file => $infile, -format => $format);
my $out= Bio::TreeIO->new(-format => $oformat, -file => ">$outfile");

while( my $t = $in->next_tree ) {
  $out->write_tree($t);
}

