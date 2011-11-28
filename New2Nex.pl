#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
    my $usage = "script.pl INFILE OUTFILE\n";
    my $infile = shift or die $usage;
    my $outfile = shift or die $usage;

my ($filein,$fileout) = @ARGV;
my ($format,$oformat) = qw(newick nexus);
my $in = Bio::TreeIO->new(-file => $infile, -format => 'newick');

my $out= Bio::TreeIO->new(-format => 'nexus', -file => ">nexus.nxs");

while( my $t = $in->next_tree ) {
  $out->write_tree($t);
}

#This doesn't seem to work particulally well - consider running perl -pi -e '$_ =~ s/(\[\])//g;' NexusTree.file on the resulting
#file to remove all the [] rubbish. Or jsut don't use nexus - it's just as flawed as newick