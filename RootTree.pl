#! /usr/bin/env perl

=head1 NAME

RootTree<.pl>

=head1 USAGE

Example: RootTree.pl -ld <SF life domain> -o <output.nwk> -t <input tree file (can be a list of trees, script will work sequenctially on all of them)> -og <comma seperated (no spaces) sting of nodes to use as ougroup>

RootTree.pl -t hs.newick -o RootedTree -og hs,gx

RootTree.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to root an unrooted (or previously rooted) tree by use of a superfamily life domain code or by a list of genome identifiers.
The latter option will still work on a tree without SUPERFAMILY taxa as nodes.

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

use Bio::TreeIO;
use IO::String;

use Supfam::SQLFunc;
use Supfam::Utils;
use Supfam::TreeFuncs;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debugNo packages will be installed, upgraded, or removed.

my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $RootingLifeDomain;
my $outputfile;
my $CLIOutgroupGenomes;
my $OutgroupGenomeFile;
my $include = 'y';

#Main Script
#----------------------------------------------------------------------------------------------------------------

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t:s" => \$TreeFile,
           "outgrouplifedomain|ld:s" => \$RootingLifeDomain,
           "output|o:s" => \$outputfile,
           "outgroup|og:s" => \$CLIOutgroupGenomes,
           "include|inc:s" => \$include,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Get a collection of nodeIDs by which to root tree#use XML::Simple qw(:strict);          #Load a config file from the local directory


$outputfile = "$TreeFile.Extract" unless(defined($outputfile));

my $OutgroupNodes = []; #This is the outgroup for the tree

#Populate $SelectedGroupNodes with outgroup;

if($RootingLifeDomain){
		
	my $dbh = dbConnect();
	my $sth;

	$sth = $dbh->prepare("SELECT genome FROM genome WHERE include = 'y' AND domain = ?;") if($include);
	$sth = $dbh->prepare("SELECT genome FROM genome WHERE domain = ?;") unless($include);
	$sth->execute($RootingLifeDomain);
	
	while (my $genome = $sth->fetchrow_array() ) { push(@$OutgroupNodes,$genome);} #Populate @$SelectedGroupNodes with genomes from the chosen superfamily life domain
	dbDisconnect($dbh);

}

if($CLIOutgroupGenomes){
	my @CLIGenomes = split(',',$CLIOutgroupGenomes);
	push(@$OutgroupNodes,@CLIGenomes) ;
}

if($OutgroupGenomeFile){
	
	open OUTGENOMES, "$OutgroupGenomeFile" or die $!;
	
	foreach my $line (<OUTGENOMES>){
		
		chomp($line);
		next if ($line =~ m/^#/); #Weed out comment lines
		
		if($line =~ m/\t/){push(@$OutgroupNodes,split("\t",$line));}else{push(@$OutgroupNodes,$line);} 
		#Input file can be any mixture of newlines and tab seperated entries, which this line catches and pushes onto @$SelectedGroupNodes
	}
	
	close OUTGENOMES;
}
	
die 'Selection must be larger than one gneome (a trvial tree). You must pass in a list of genomes so as to serve as a rooting outgroup!'."\n" unless (scalar(@$OutgroupNodes)); # A little error checking


#Read in a file of newick trees which to root. One tree per line
open TREEFILE, "<$TreeFile";

my $out = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">$outputfile");

while (my $UnrootedTree = <TREEFILE>){

	last if($UnrootedTree =~ m/^#/ || $UnrootedTree !~ m/;$/ ); #Weed out comments or lines that don't have a trailing ";" carachter. God knows why these are in there, but you never know ...
			
	chomp($UnrootedTree); #Remove trailing newline carachter
	my $io = IO::String->new($UnrootedTree);
	my $treeio = Bio::TreeIO->new(-fh => $io,
	                              -format => 'newick');
	                          
	my $tree = $treeio->next_tree;
	my $root = $tree->get_root_node;
	
	my @TreeLeafNodeIDs = grep{$_->is_Leaf}$root->get_all_Descendents;

	#Convert the outgroup node names into ids
	my $OutgroupNodeIDs = [];
	
	my $TreeLeafNodeIDDictionary = {};
	map{$TreeLeafNodeIDDictionary->{$_->id} = $_}@TreeLeafNodeIDs;
	my ($Union,$IntersectionOutgroupTreeLeaves,$Outgroup,$TreeExclusive) = IntUnDiff($OutgroupNodes,[keys(%$TreeLeafNodeIDDictionary)]);
	
	map{my $OutgroupNodeID =  $TreeLeafNodeIDDictionary->{$_}; push(@$OutgroupNodeIDs,$OutgroupNodeID);}@$IntersectionOutgroupTreeLeaves;

	die "Outgroup not in tree!" unless(scalar(@$OutgroupNodeIDs));

	FindTrueRoot($tree,$OutgroupNodeIDs); #Operates on $tree and roots it

	$out->write_tree($tree);
}

close TREEFILE;


__END__

