#! /usr/bin/perl -w

=head1 NAME

TreeIntersection<.pl>

=head1 USAGE

 TreeIntersection.pl -t1 <tree1.file> -t2 <tree2.file>

=head1 SYNOPSIS

A script to take two tree files and output two new treefiles, both of which contain the same nodes (i.e. they can only differ on topology and branch lengths, not on leaves).

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
use lib "/home/sardar/bin/perl-libs-custom/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Supfam::Utils;
use Supfam::TreeFuncs;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $treeA;
my $treeB;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree1|t1=s" => \$treeA,
           "tree2|t2=s" => \$treeB,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my ($stringTreeA, $stringTreeB);

open FHA, "<$treeA"; open FHB, "<$treeB";

while (my $line = <FHA>){
	$stringTreeA = $line;
	last unless ($line =~ m/^#/); #Only read in the first non comment line
}

while (my $line = <FHB>){
	$stringTreeB = $line;
	last unless ($line =~ m/^#/); #Only read in the first non comment line
}

close FHA; close FHB;

#This framework of storing trees as as an io stream was chosen so as to allow for future expansion

my $treeAio = IO::String->new($stringTreeA); 
my $treeBio = IO::String->new($stringTreeB);

my ($TreeAObject,$TreeBObject) = TreeIntersection($treeAio,$treeBio,1) ;

my $TreeAOut = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">$treeA.isotre");
$TreeAOut->write_tree($TreeAObject);                     

my $TreeBOut = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">$treeB.isotre");
$TreeBOut->write_tree($TreeBObject); 

__END__
