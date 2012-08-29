#! /usr/bin/env perl

=head1 NAME

TidyUpTreeI<.pl>

=head1 USAGE

 TidyUpTree.pl [options -v,-d,-h] -t --tree TreeFile -b --branches -i --internal
 
=head1 SYNOPSIS

A simple clean up script for a newick tree. I wrote this because I was sick of RAxML pumping out errors for odly formatted trees.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use Modern::Perl;

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
use Carp;
use Supfam::Utils;
use Supfam::TreeFuncsNonBP;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $tree;
my $branchesflag = 0;
my $internalnodeflag = 0;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$tree,
           "internal|i!" => \$internalnodeflag,
           "branches|b!" => \$branchesflag,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

croak "You need to specify a treefile!\n" unless($tree);

open FH, "<$tree" or die $?;

my $NewickStringOfTree = <FH>;

close FH;

my ($root,$TreeCacheHash) = BuildTreeCacheHash($NewickStringOfTree);
print STDERR "Built TreeCacheHash\n";

sanitise_TreeHash($TreeCacheHash,$root);

my $Output = ExtractNewickSubtree($TreeCacheHash,$root,$branchesflag,$internalnodeflag);

print $Output."\n";


__END__