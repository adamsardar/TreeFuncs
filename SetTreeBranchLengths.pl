#! /usr/bin/env perl

=head1 NAME

SetTreeBranchLengthsI<.pl>

=head1 USAGE

 SetTreeBranchLengthsI.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

Tiny script. Sets all branches in input tree to a value specified usign -bl flag. Prints tree to STDOUT

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2012 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom";


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
use Supfam::TreeFuncsNonBP;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $BranchLength;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$TreeFile,
           "branchlength|bl=f" => \$BranchLength,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

open TREE, "<$TreeFile" or die $!.$?;
my $TreeString = <TREE>;
close TREE;
	
my ($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

my @AllNodes = @{$TreeCacheHash->{$root}{'all_Descendents'}};
push(@AllNodes,$root);

map{$TreeCacheHash->{$_}{'branch_length'} = $BranchLength}@AllNodes;

my $NewickTree = ExtractNewickSubtree($TreeCacheHash, $root,1,0);
print STDOUT $NewickTree."\n";


__END__