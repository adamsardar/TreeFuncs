#! /usr/bin/perl -w

=head1 NAME

DolloParsimonyNormaliseTreeI<.pl>

=head1 USAGE

 DolloParsimonyNormaliseTree.pl [options -v,-d,-h] (-rl root left | -rr root right)

=head1 SYNOPSIS

Extracts the SUEPRFAMILY tree beneath a certain root (defined by left_id or right_id) and then dumps out two trees - the first is simply the SUPERFAMILY tree and the second is that same tree with branch lengths equal to the number of domain architectures that have been lost.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
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

use Supfam::Utils;
use Supfam::TreeFuncsNonBP;
use Supfam::SQLFunc;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my $rootleft;  #Root node from the tree (left id in SQL table) to be used.
my $rootright; #Root node from the tree (right id in SQL table) to be used.

my $TreeByCreations;
my $TreeByDeletions;
my $TreeByDeletionsAndCreations;

my $ExcludeGenomeFile = 0;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
            "treebydeletions|tbc:i" => \$TreeByCreations,
           "treebycreations|tbd:i" => \$TreeByDeletions,
           "treebydelsandcreats|tbcd:i" => \$TreeByDeletionsAndCreations,
           "rootleft|rl:s" => \$rootleft,
            "rootright|rr:s" => \$rootright,
            "exludegenomes|eg:s" => \$ExcludeGenomeFile,   
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my $dbh = dbConnect();
my $sth;

if(defined($rootleft) && ! defined($rootright)){
	
	$sth = $dbh->prepare("SELECT tree.right_id FROM tree WHERE tree.left_id = ?;");
	$sth->execute($rootleft);
	$rootright = $sth->fetchrow_array();
	
}elsif(defined($rootright) && ! defined($rootleft)){
	
	$sth = $dbh->prepare("SELECT tree.left_id FROM tree WHERE tree.right_id = ?;");
	$sth->execute($rootright);
	$rootleft = $sth->fetchrow_array();
	
}else{
	
	die 'Need to provide either the left_id or right_id of root node, preferably both! ';	
}

$sth->finish;
dbDisconnect($dbh);


my $ExcludedGenomes = [];


if($ExcludeGenomeFile){
	
	die "Feature not supported at current!\n";
	
	open EXCGENOMES, "$ExcludeGenomeFile" or die $!;
	
	foreach my $line (<EXCGENOMES>){
		
		chomp($line);
		next if ($line =~ m/^#/); #Weed out comment lines
		
		if($line =~ m/\t/){push(@$ExcludedGenomes,split("\t",$line));}else{push(@$ExcludedGenomes,$line);} 
		#Input file can be any mixture of newlines and tab seperated entries, which this line catches and pushes onto @$SelectedGroupNodes
	}
	
	close EXCGENOMES;
}




my $SQLTreeCacheHash = supfamSQL2TreeHash($rootleft,$rootright);

my $root = $rootleft;

my $SupfamNewickTree = ExtractNewickSubtree($SQLTreeCacheHash,$root,1,0);

open RAWTREE, "> RawSupfamTree.nwk" or die $!;
print RAWTREE $SupfamNewickTree."\n";
close RAWTREE;


normailseSQLTreeHashBranchForDeletions($SQLTreeCacheHash);

my $DolloTree = ExtractNewickSubtree($SQLTreeCacheHash,$root,1,0);

open RAWTREE, "> DolloDeletionsNormalisedTree.nwk" or die $!;
print RAWTREE $DolloTree."\n";
close RAWTREE;



__END__