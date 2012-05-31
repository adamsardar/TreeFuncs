#!/usr/bin/env perl


=head1 NAME

SupfamiliesInGenomes-AllvsAllDistances<.pl>

=head1 USAGE

SupfamiliesInGenomes-AllvsAllDistances<.pl>
 
 Example usage: 
 SupfamiliesInGenomes-AllvsAllDistances.pl (no arguments)

=head1 SYNOPSIS

A script to extract data regarding the inter-genome distances and the distribution of superfamilies within them.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2012 Gough Group, University of Bristol.

=head1 EDIT HISTORY

See 

=cut

#----------------------------------------------------------------------------------------------------------------
# Strict Pragmas
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
B<Math::Random> Used in Monte Carlo simulations steps
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__

use Supfam::TreeFuncsNonBP; #Needed for MRCA calculations and summing of branch lengths
use File::Temp; #Needed to create temporary files in RAM to write results to
use Time::HiRes; #Used to track lenght of script

use Parallel::ForkManager; # Used in create forking jobs out and making thins run a bit speedier

use Math::Combinatorics; # Needed to construct all pairs of genomes

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $TotalTic = Time::HiRes::time; #USed in timing the total runtime of script

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $OutputFilename = 'AllVsAllGenomeTreeDistances';
my $TreeFile;
my $maxProcs = 0;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$OutputFilename,
           "tree|i:s" => \$TreeFile,
           "processor_cores|p:i" => \$maxProcs,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
#---------------------------------------------------------------------------------------------------------------
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

#---------------------------------------

`mkdir /dev/shm/temp` unless (-d '/dev/shm/temp');
my $RAMDISKPATH = '/dev/shm/temp';
#Path to a piece of RAM that we can write to. This could be on hard disk, but on *nix systems /dev/shm is a shared RAM disk folder. We want a temporary folder in her that will be cleaned up on Prgram exit

# Make a temporary directory for output data 
my $RAWPATH = File::Temp->newdir( DIR => $RAMDISKPATH , CLEANUP => 1) or die $!;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

#Produce a tree hash, either from SQL or a provided treefile
my ($root,$TreeCacheHash);

if($TreeFile){
	
	open TREE, "<$TreeFile" or die $!.$?;
	my $TreeString = <TREE>;
	close TREE;

	($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

}else{
	
	die "no tree file provided as tree to calculate genome distances\n";	
}


open TREEINFO, ">AllVsAllGenomesDistances_tree.nwk";
my $NewickTree = ExtractNewickSubtree($TreeCacheHash, $root,1,0);
print TREEINFO $NewickTree."\n";
close TREEINFO;
# Dump tree into an output file

#--------------------------------------------------
# All n choose 2 pairs of genomes into a job hash

my @TreeGenomes = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree
my @TreeGenomesNodeIDs = @{$TreeCacheHash->{$root}{'Clade_Leaves'}}; # All of the genomes (leaves) of the tree


my $NoOfForks = $maxProcs;
$NoOfForks = 1 unless($maxProcs);

my $ForkJobsHash = {};

my $combinat = Math::Combinatorics->new(count => 2,
                                          data => [@TreeGenomesNodeIDs],
                                         );
my $jobcount= 0;

while (my @combo = $combinat->next_combination){
	
	my $ForkIndex = $jobcount%$NoOfForks;
	$ForkJobsHash->{$ForkIndex}=[] unless(exists($ForkJobsHash->{$jobcount%$NoOfForks}));
	my $ForkArray = $ForkJobsHash->{$ForkIndex};
	push(@$ForkArray,\@combo);
	$jobcount++;
}#Create lists of jobs to be done by the relative forks

print STDERR "No Genome Pairs in job batch is approx: ".scalar(@{$ForkJobsHash->{1}})."\n";

#----------------------------------------------------

#Main-loop-------------------------------------------

my $pm = new Parallel::ForkManager($maxProcs) if ($maxProcs);# Initialise

foreach my $fork (0 .. $NoOfForks-1){
	
	my $JobList = $ForkJobsHash->{$fork};
		
	# Forks and returns the pid for the child:
if ($maxProcs){$pm->start and next};
		
	open OUT, ">$RAWPATH/$OutputFilename".$$.".-RawData.colsv" or die "Can't open file $RAWPATH/$OutputFilename".$!;

	foreach my $GenomePairArrayRef (@$JobList){
	
		die "Something wrong with genome combinatorics step - not two genomes in the pair!!\n\n" unless(scalar(@$GenomePairArrayRef) == 2);
		
		my ($GenomeA,$GenomeB) = @$GenomePairArrayRef;
		
		my $MRCA = FindMRCA($TreeCacheHash,$root,[$GenomeA,$GenomeB]);
		
		my $DistanceA = MRCADistanceSum($TreeCacheHash, $GenomeA, $MRCA);
		my $DistanceB = MRCADistanceSum($TreeCacheHash, $GenomeB, $MRCA);
		my $BranchSum = $DistanceA + $DistanceB;
		
		print OUT $TreeCacheHash->{$GenomeA}{'node_id'}."\t";
		print OUT $TreeCacheHash->{$GenomeB}{'node_id'}."\t";
		print OUT $BranchSum."\n";
	}
	
	#Output to data file
	
	close OUT;
	
$pm->finish if ($maxProcs); # Terminates the child process

}

print STDERR "Waiting for Children...\n";
$pm->wait_all_children if ($maxProcs);
print STDERR "Everybody is out of the pool!\n";


`cat $RAWPATH/* > ./$OutputFilename.dists`;

my $TotalToc = Time::HiRes::time;

my $TotalTimeTaken = ($TotalToc-$TotalTic);

print STDERR $TotalTimeTaken." seconds\n";

#-------
__END__


