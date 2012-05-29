#! /usr/bin/env perl

=head1 NAME

TreeIntersection<.pl>

=head1 USAGE

 TreeIntersection.pl -t1 --tree1 <tree1.file> (-t2 --tree2 <tree2.file> | -l --list new_line_seperated_list_of_nodes.tx) [-b flag for branch legnths] [-vi verbose intersection - gives more information during the intersection process] 

=head1 SYNOPSIS

A script to take two tree files and output two new treefiles, both of which contain the same nodes (i.e. they can only differ on topology and branch lengths, not on leaves). Another option
is to pass in a new line seperated lsit of genomes to be included in the output file. don't forget to specify the branches flag (-b) if you wish to keep branch lengths.

The speed is ... acceptable. For a 1200 leaf tree intersecting with a list of 20 genomes (so nearing the worst case scenario) it took 5 minutes. If you need it to be faster, give me a shout on the github page and maybe I can speed it up.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

See https://github.com/adamsardar/TreeFuncs/commits/master for the package history. Also, this is the appropriate place to raise bugs.

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

use Supfam::Utils;
use Carp;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $treeA;
my $treeB;
my $branchesflag;
my $internalnodesflag;
my $verboseintersection = 1;
my $listofgenomes;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree1|t1=s" => \$treeA,
           "tree2|t2:s" => \$treeB,
           "list|l:s" => \$listofgenomes,
           "internal|i!" => \$internalnodesflag,
           "branches|b!" => \$branchesflag,
           "verboseintersection|vi!" => \$verboseintersection,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
        
#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my ($stringTreeA, $stringTreeB);

my $TreeAOnly = 0; #Flag to prevent analysis and output of the second tree. Flags as one in the case where a list of genomes were passed in

open FHA, "<$treeA"  or die $?; 

while (my $line = <FHA>){
	$stringTreeA = $line;
	last unless ($line =~ m/^#/); #Only read in the first non comment line
}

close FHA; 

if($treeB){
	
	open FHB, "<$treeB"  or die $?;
	
	while (my $line = <FHB>){
		$stringTreeB = $line;
		last unless ($line =~ m/^#/); #Only read in the first non comment line
	}
	
	close FHB;
	
}elsif($listofgenomes){
	
	my @GenomesList;
	
	open LIST, "<$listofgenomes"  or die $?;
	
	while (my $line = <LIST>){
		chomp($line);
		die "Input contains a semicolon ; - a dissallowed newick carachter on line $. - $line!\n" if ($line =~ m/;{1}/);
		die "Input contains a bracket ( or ) - a dissallowed newick carachter on line $. - $line!\n" if ($line =~ m/[\(\)]{1}/);
		push(@GenomesList,$line);		
	}
	close LIST;
	
	$stringTreeB = join(',',@GenomesList);
	$stringTreeB = '('.$stringTreeB.');'; #Create a flat tree out of the genome list provided
	
	$TreeAOnly = 1;
}


my ($TreeAObject,$Aroot,$TreeBObject,$Broot) = TreeIntersection($stringTreeA,$stringTreeB,$verboseintersection,$TreeAOnly) ;

open OUTA, ">$treeA.isotree" or die $?;
my $StringA = ExtractNewickSubtree($TreeAObject,$Aroot,$branchesflag,$internalnodesflag);
print OUTA $StringA."\n";
close OUTA; 

unless($TreeAOnly){

	open OUTB, ">$treeB.isotree"  or die $?;
	my $StringB = ExtractNewickSubtree($TreeBObject,$Broot,$branchesflag,$internalnodesflag);
	print OUTB $StringB."\n";
	close OUTB;
}

__END__
