#! /usr/bin/perl -w

=head1 NAME

supfam2Hennig86<.pl>

=head1 USAGE

 supfam2Hennig86.pl [options -v,-d,-h] -t --tree <TreeFile in Newick> -s --traits <genome archs file, as per using RaxML>
 -o --output <outputfile name>

=head1 SYNOPSIS

A script to take the binary prescence/abscene trait vector file used in tree construction (under say RaxML) and convert it to the similar, but
frustratingly slightly different TNT format. The output will be by species, as specified in an inout newick tree.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

29-July-2011 Intiial Entry

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

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
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::SQLFunc;
use Bio::TreeIO;
use IO::String;
use Supfam::Utils;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $genome_archs_file;
my $outputfile;


# Sub definitions
#----------------------------------------------------------------------------------------------------------------

#Main Script
#----------------------------------------------------------------------------------------------------------------


#Set command line flags and parameters.n
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t:s" => \$TreeFile,
           "traits|s:s" => \$genome_archs_file, # -s is for RazML like compatibility
           "output|o:s" => \$outputfile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Get a collection of nodeIDs by which to root tree

my $input = new Bio::TreeIO(-file   => "$TreeFile",
                            -format => "newick") or die $!;
                            
my $tree = $input->next_tree;
my $root = $tree->get_root_node;
my @RootDescendents = $root->get_all_Descendents;

#Extract all the species in the treefile given
my @TreeTaxa = map{$_->id}grep{$_->is_Leaf}@RootDescendents;

#Create a hash of all the trait vectors per taxon

my $TraitHash = {};
#$TraitHash -> {taxon => binary traits}

open TRAITS, "<$genome_archs_file";

my $LineLength;

while (<TRAITS>){

my $line = $_;

	if($line =~ m/^\w{2,3}\s/i){
		
		chomp($line);
		my @SplitLine = split(/\W/,$line); # (Taxa, Trait list as string)
		$TraitHash->{$SplitLine[0]}=$SplitLine[1]; #$TraitHash -> {taxon => Trait list as string}
		$LineLength = length($SplitLine[1]); #Wastful way to store the length of a line I know, but meh
	}
}

close TRAITS;


#Wrtie only the records for species in the tree to file
	#Add some preamble to the output file

my $NoTaxa = scalar(keys(%$TraitHash));

open OUT, ">$outputfile";
	
print OUT "xread\n'This is a carachter trait file in Hennig86 format, using domain architectures (or supra domains, depending on the input file). This file was created using supfam2Hennig86'\n";
print OUT "$LineLength $NoTaxa\n";

foreach my $Taxon (keys(%$TraitHash)){
	
	print OUT $Taxon;
	print OUT "\t";
	print OUT $TraitHash->{$Taxon};
	print OUT "\n";
}

#print OUT "\n;";
close OUT;

__END__

