#! /usr/bin/env perl

=head1 NAME

ExtractSubtree<.pl>

=head1 USAGE

ExtractSubtree.pl -t <input tree> -ld <rooting life domain> -o <Output file>

ExtractSubtree.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

This script cuts out a section of the input newick tree(s) and outputs as a new newick file, as decided by a list of given taxa.
Note that this includes all other members of the tree below the MRCA of the set of taxa given, unless the --remove flag is set to 1.
The new trees are rooted on the MRCA of all members of the specified group (this can be a tab-seperated and/or newline seperated list as a file, genomes passed in through CLI or
simply a SUPERFAMILY database life domain).

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
use DBI;
use Supfam::SQLFunc;
use Supfam::TreeFuncs;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use IO::String;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $LifeDomain;
my $outputfile;
my $remove = 0; #Flag for removal of genomes not belonging to the extract genomes requested. Default is to include all genomes below the MRCA in the tree specified of all genomes asked.
my $OutgroupGenomeFile;
my @CLIOutgroupGenomes;
my $outputformat = 'newick';
my $include = 1; #flag as to whether or not to include the 'include ='y'' genomes from the SUPERFAMILY database. This is only relevent to the life-domain option

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t:s" => \$TreeFile,
           "outgrouplifedomain|ld:s" => \$LifeDomain,
           "includeygenomes|iyg:s" => \$include,
           "outgroupgenome|og:s" => \@CLIOutgroupGenomes,
           "outgroupgenomefile|ogf:s" => \$OutgroupGenomeFile,
           "output|o:s" => \$outputfile,
           "outputformat|of:s" => \$outputformat,
           "remove|r:i" => \$remove,
        ) or die "Fatal Error: Problem parsing command-line ".$!;



$outputfile = "$TreeFile.Extract" unless(defined($outputfile));

my $SelectedGroupNodes = []; #This is the ingroup for the tree

#Populate $SelectedGroupNodes with outgroup;

if($LifeDomain){
		
	my $dbh = dbConnect();
	my $sth;

	$sth = $dbh->prepare("SELECT genome FROM genome WHERE include = 'y' AND domain = ?;") ;#if($include);
	#$sth = $dbh->prepare("SELECT genome FROM genome WHERE domain = ?;") unless($include);
	$sth->execute($LifeDomain);
	
	while (my $genome = $sth->fetchrow_array() ) { push(@$SelectedGroupNodes,$genome);} #Populate @$SelectedGroupNodes with genomes from the chosen superfamily life domain
	dbDisconnect($dbh);

}

push(@$SelectedGroupNodes,@CLIOutgroupGenomes) if(scalar(@CLIOutgroupGenomes));

if($OutgroupGenomeFile){
	
	open OUTGENOMES, "$OutgroupGenomeFile" or die $!;
	
	foreach my $line (<OUTGENOMES>){
		
		chomp($line);
		next if ($line =~ m/^#/); #Weed out comment lines
		
		if($line =~ m/\t/){push(@$SelectedGroupNodes,split("\t",$line));}else{push(@$SelectedGroupNodes,$line);} 
		#Input file can be any mixture of newlines and tab seperated entries, which this line catches and pushes onto @$SelectedGroupNodes
	}
	
	close OUTGENOMES;
}
	
die 'Selection must be larger than one gneome (a trvial tree). You must pass in a list of genomes, whose MRCA will serve as the root of the output tree!'."\n" unless (scalar(@$SelectedGroupNodes) >  1); # A little error checking

#Read in a file of newick trees which to root. One tree per line
open TREEFILE, "<$TreeFile";

my $out = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">$outputfile");
                           

while (my $FullTree = <TREEFILE>){
	
	last if($FullTree =~ m/^#/ || $FullTree !~ m/;$/ ); #Weed out comments or lines that don't have a trailing ";" carachter. God knows why these are in there, but you never know ...
			
	chomp($FullTree); #Remove trailing newline carachter
	my $io = IO::String->new($FullTree);
	my $treeio = Bio::TreeIO->new(-fh => $io,
	                              -format => 'newick');
	                           
	my $tree = $treeio->next_tree;
	
	my $subtree = ExtractSubtree($tree,$SelectedGroupNodes,$remove);

	$out->write_tree($subtree);	#Thanks bioperl, you make everything easier!
}

close TREEFILE;


__END__

