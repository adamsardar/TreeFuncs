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
use DBI;

use Supfam::Utils;
use Supfam::TreeFuncsNonBP;
use Supfam::DolloParsmony;
use Supfam::RAxML_Ancestral_States_Parser;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $LeafSpeciesStatesFile;
my $TreeByCreations;
my $TreeByDeletions;
my $TreeByDeletionsAndCreations;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "speciestraits|st=s" => \$LeafSpeciesStatesFile, #trait file used to build RAxMl tree
           "tree|t=s" => \$TreeFile,
            "treebydeletions|tbc:i" => \$TreeByCreations,
           "treebycreations|tbd:i" => \$TreeByDeletions,
           "treebydelsandcreats|tbcd:i" => \$TreeByDeletionsAndCreations,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------


open FH, "<$TreeFile" or die $?;

my $NewickStringOfTree = <FH>;

close FH;

my ($root,$TreeCacheHash) = BuildTreeCacheHash($NewickStringOfTree);
print STDERR "Built TreeCacheHash\n";

DolloParsimonyAncestralState($TreeCacheHash,$root,$LeafSpeciesStatesFile,10); #Using 10 threads
print STDERR "Parsed DOLLOP Ancestral States\n";

my $NumberOfTraits = scalar(keys(%{$TreeCacheHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}}));

DOLLOP_Ancestral_Trait_Changes_in_Clade($TreeCacheHash,$root,$root);

print STDERR $TreeCacheHash->{$root}{'DOLLOP_Total_Number_Created'}." = Total DOLLOP Created below root\n";
print STDERR $TreeCacheHash->{$root}{'DOLLOP_Total_Number_Deleted'}." = Total DOLLOP Deleted below root\n";

my $Tree = ExtractNewickSubtree($TreeCacheHash,$root,1,0);


foreach my $TreeNode (keys(%{$TreeCacheHash})){
	
	if($TreeByDeletions){
	
		$TreeCacheHash->{$TreeNode}{'branch_length'} = $TreeCacheHash->{$TreeNode}{'DOLLOP_Number_Deleted'}/$NumberOfTraits;
	
	}elsif($TreeByCreations){
		
		$TreeCacheHash->{$TreeNode}{'branch_length'} = $TreeCacheHash->{$TreeNode}{'DOLLOP_Number_Created'}/$NumberOfTraits;
	
		
	}elsif($TreeByDeletionsAndCreations){
		
		$TreeCacheHash->{$TreeNode}{'branch_length'} = $TreeCacheHash->{$TreeNode}{'DOLLOP_Number_Created'}+$TreeCacheHash->{$TreeNode}{'DOLLOP_Number_Deleted'}/$NumberOfTraits;
		
	}else{
		
		die "Need to choose a method by which to normalise tree"
	}
}

$Tree = ExtractNewickSubtree($TreeCacheHash,$root,1,0);

print $Tree."\n";

__END__