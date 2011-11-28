#! /usr/bin/perl -w

=head1 NAME

I<.pl>

=head1 USAGE

 .pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to take a list of genomes (strains or whatever) and output some summary data regarding the domain architectures in common.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "/home/sardar/bin/perl-libs-custom";


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
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::Utils;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use IO::String;
use Time::HiRes;

use Bio::Tree::TreeI;
use Supfam::TreeFuncs;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $treefile;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$treefile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Sub definitions
#----------------------------------------------------------------------------------------------------------------
=head1 DESCRIPTION

Detailed info about the script goes here

=head2 Methods
=over 4
=cut

=item * func
Function to do something
=cut

# Main Script Content
#----------------------------------------------------------------------------------------------------------------





my $out = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">ConstraintTrees.nwk");
                           

my ($root,$TreeCacheHash,$tree) = BuildTreeCacheHash($treefile); # A massive limitation of this script is the use of BioPerl TreeIO which is SLOOOOW. This is a lookup hash to speed things up.

my @TreeLeaves = grep{$TreeCacheHash->{$_}{'is_Leaf'} == 1}@{$TreeCacheHash->{$root}{'all_Descendents'}};

foreach my $leaf (@TreeLeaves){
	

	my @LeavesOmmitingOne = grep{$_ ne $leaf}@TreeLeaves;
	print scalar(@LeavesOmmitingOne)."\n";
	
	my @LeavesOmmitingOneIDs = map{$TreeCacheHash->{$_}{'node_id'}}@LeavesOmmitingOne;
	
	my $subtree = ExtractSubtree($tree,\@LeavesOmmitingOneIDs,1);
	$out->write_tree($subtree);	#Thanks bioperl, you make everything easier!
}

__END__