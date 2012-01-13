#! /usr/bin/env perl

=head1 NAME


AncestralSQL2Newick<.pl>

=head1 USAGE

 AncestralSQL2Newick.pl [options -v,-d,-h] <ARGS>
 
 Example usage:
 
 AnceSQL2Newick.pl -rl 3 -rr 1364 -o Ancestral
 
 To extract genomes from the root of all eukaryotes.
 Don't forget! Bacterial dollo parsmiony results are not present in SUPERFAMILY, as dollo parsimony is a very poor model for bacteria.

=head1 SYNOPSIS

A script to parse an the ancestral_info table from SUPERFAMILY and output a tree with bracnh lengths proportional to deletions in newick format. Simply specifiy the root of the tree and this script will
sort out the rest.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
use DBI;
use Supfam::Utils;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;



# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $rootleft;  #Root node from the tree (left id in SQL table) to be used.
my $rootright; #Root node from the tree (right id in SQL table) to be used.
my $deletiontreeflag; #Produce a tree with branch lengths equal to the number of domain architecture deltions along that edge

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "rootleft|rl=i" => \$rootleft,
           "rootright|rr=i" => \$rootright,
           "deletions|del:i" => \$deletiontreeflag,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
        
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;
        

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my $dbh = dbConnect();
my $sth;

if(defined($rootleft) && ! defined($rootright)){
	
	$sth = $dbh->prepare("SELECT tree.right_id FROM tree WHERE tree.left_id = ?;");
	$sth->execute($rootleft);
	$rootright = $sth->fetchrow_array();
	
}elsif(defined($rootright) && ! defined($rootleft)){
	
	$sth = $dbh->prepare("SELECT tree.left_id FROM tree WHERE tree.right_id = ?;");
	$sth->execute($rootright);
	$rootleft = $sth->fetchrow_array();
	
}else{
	
	die 'Need to provide either the left_id or right_id of root node, preferably both! ';	
}

print $rootleft." Root Left ".$rootright." Root Right\n";

my $SQLTreeCacheHash = supfamSQL2TreeHash($rootleft,$rootright);

EasyDump('./Dump.dat',$SQLTreeCacheHash);

normailseSQLTreeHashBranchForDeletions($SQLTreeCacheHash);

my $DeletionsNormailsedTree = ExtractNewickSubtree($SQLTreeCacheHash,$rootleft,1,0);

print $DeletionsNormailsedTree."\n";

__END__

