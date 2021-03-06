#!/usr/bin/env perl

use Bio::TreeIO;
use Modern::Perl;

my ($NEWICKFILE, $NEXUSFILE) = @ARGV;
print "newickfile=$NEWICKFILE, nexusfile=$NEXUSFILE\n";

my $treeio = new Bio::TreeIO(-format => 'newick', -file   => "$NEWICKFILE");
my $treeout = new Bio::TreeIO(-format => 'nexus', -file   => ">$NEXUSFILE");

while(my $tree = $treeio->next_tree) {
        $treeout->write_tree($tree);
}
