#! /usr/bin/perl -w

=head1 NAME

TreeSQL2Newick<.pl>

=head1 USAGE

TreeSQL2Newick.pl [options -v,-d,-h] <ARGS> -rl|--rootleft node_left_id -rr|--rootright node_right_id 

=head1 SYNOPSIS

A script to parse an sql table from SUPERFAMILY and output a tree in newick format. Simply specifiy the root of the tree and this script will
sort out the rest. Specify the root by passing in eith the left or the right id of the root node, or better still both.

Outputs tree as a newick file.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "/home/sardar/workspace/Oates/lib/";

# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Supfam::TreeFuncs;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $rootleft;  #Root node from the tree (left id in SQL table) to be used.
my $rootright; #Root node from the tree (right id in SQL table) to be used.
my $deletiontreeflag; #Produce a tree with branch lengths equal to the number of domain architecture deltions along that edge
my @excluded; # A list of genomes not wanted in this tree - not yet implemented
my $OutputFile = 'fileout.dat';

my $ExcludeLifeDomain;
my @CLIExlucdeGenomes;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "rootleft|rl:s" => \$rootleft,
           "rootright|rr:s" => \$rootright,
           "exclude|ex:s" => \@excluded,
           "output|o:s" => \$OutputFile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
        
# Main Script Content
#----------------------------------------------------------------------------------------------------------------


#Populate $ExcludedGenomes with genomes not wanted in the resulting tree;

my $ExcludedGenomes = [];

if($ExcludeLifeDomain){
		
	my $dbh = dbConnect();
	my $sth;

	$sth = $dbh->prepare("SELECT genome FROM genome WHERE include = 'y' AND domain = ?;");
	$sth->execute($ExcludeLifeDomain);
	
	while (my $genome = $sth->fetchrow_array() ) { push(@$ExcludedGenomes,$genome);} #Populate @$ExcludedGenomes with genomes from the chosen superfamily life domain
	dbDisconnect($dbh);

}

push(@$ExcludedGenomes,@CLIExlucdeGenomes) if(scalar(@CLIExlucdeGenomes));

if($ExcludeGenomeFile){
	
	open EXCGENOMES, "$ExcludeGenomeFile" or die $!;
	
	foreach my $line (<EXCGENOMES>){
		
		chomp($line);
		next if ($line =~ m/^#/); #Weed out comment lines
		
		if($line =~ m/\t/){push(@$ExcludedGenomes,split("\t",$line));}else{push(@$ExcludedGenomes,$line);} 
		#Input file can be any mixture of newlines and tab seperated entries, which this line catches and pushes onto @$SelectedGroupNodes
	}
	
	close EXCGENOMES;
}

##Yet to actually be implemented in the SQL2Newick: allow for certiain genomes to be excluded!

my ($RightTreeHash,$LeftTreeHash) = SQL2Newick($rootleft,$rootright);

open OUTPUT, "> " or die $!;
my $FullTree = $LeftTreeHash->{$rootleft}[1][2];

print OUTPUT $FullTree."\n";
close OUTPUT;

__END__

