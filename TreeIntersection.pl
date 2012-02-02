#! /usr/bin/env perl

=head1 NAME

TreeIntersection<.pl>

=head1 USAGE

 TreeIntersection.pl -t1 --tree1 <tree1.file> (-t2 --tree2 <tree2.file> | -l --list new_line_seperated_list_of_nodes.tx)

=head1 SYNOPSIS

A script to take two tree files and output two new treefiles, both of which contain the same nodes (i.e. they can only differ on topology and branch lengths, not on leaves). Another option
is to 
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
#use lib "$ENV{HOME}/bin/perl-libs-custom/";

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
#use Supfam::hgt;
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
           "internal|i:i" => \$internalnodesflag,
           "branches|b:i" => \$branchesflag,
           "verboseintersection|vi:i" => \$verboseintersection,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


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
