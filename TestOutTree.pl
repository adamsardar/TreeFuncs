#! /usr/bin/env perl

=head1 NAME

TestOutTree<.pl>

=head1 USAGE

 TestOutTree.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A simple script to aid in development of the Supfam::TreeFuncsNonBP module

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
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::TreeFuncsNonBP;
use Supfam::Utils;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $file;
my $LeafStatesFile;
my $NThreads = 10;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "file|f=s" => \$file,
           "help|h!" => \$help,
           "Nthreds|T:i" => \$NThreads, 
           "LeafTraits|lt=s" => \$LeafStatesFile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

open TREE, "<$file" or die;

my $TreeString = <TREE>;

my ($root,$TreeHash) = BuildTreeCacheHash($TreeString);

close TREE;

assignLeftRightIDs2TreeHash($TreeHash,$root);

EasyDump('Dump.oyt',$TreeHash);



#Convert to newick string and print to screen




__END__