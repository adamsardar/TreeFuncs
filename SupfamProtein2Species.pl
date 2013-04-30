#! /usr/bin/env perl

=head1 NAME

I<SupfamProtein2Species.pl.pl>

=head1 USAGE

 SupfamCodes2FullSpecies.pl.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to take a tree, pull out all of it's leaves and then, if the leaf id is a protein id in superfamily, replace its name with the species or domain (of life) information

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2013 Gough Group, University of Bristol.

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
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
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use DBI;
use List::MoreUtils qw{uniq};

use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;

use Carp;
use Carp::Assert;
use Carp::Assert::More;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $treefile; #Traits file in phylip format
my $GenomeFlag = 0;
my $DomainFlag = 0;
my $NameFlag = 0;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$treefile,
           "genome|g!" => \$GenomeFlag,
           "domain|d!" => \$DomainFlag,
           "name|n!" => \$NameFlag,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------


open TREE, "<$treefile" or die;


my $dbh = dbConnect();
my $sth;

my $numberoftreeschanged = 0;

if($DomainFlag){

	$sth = $dbh->prepare("SELECT DISTINCT(genome.domain),protein.protein  FROM genome JOIN protein ON genome.genome = protein.genome WHERE protein.protein = ? AND (genome.include = 'y' OR genome.include = 's');");

}elsif($GenomeFlag){
	
	$sth = $dbh->prepare("SELECT DISTINCT(genome.genome),protein.protein  FROM genome JOIN protein ON genome.genome = protein.genome WHERE protein.protein = ? AND (genome.include = 'y' OR genome.include = 's');");

}elsif($NameFlag){
	
	$sth = $dbh->prepare("SELECT DISTINCT(genome.name),genome.domain  FROM genome JOIN protein ON genome.genome = protein.genome WHERE protein.protein = ? AND (genome.include = 'y' OR genome.include = 's');");
	
}else{
	
	die "You need to specify which mode you wish the script to run in";
}

my $dictionary = {};
#going to dump out a tab space mapping fiel, so that you can play with changing the input alignmenet using dictionary convert script


while (my $TreeString = <TREE>){

	my ($root,$TreeHash) = BuildTreeCacheHash($TreeString);
	
	my @TreeLeaveIDs = @{$TreeHash->{$root}{'Clade_Leaves'}};
	
	my @UniqTreeLeaveIDs = uniq(@TreeLeaveIDs);
	
	print STDERR "Number of  TreeLeaves =".scalar(@TreeLeaveIDs)."\n" if($verbose);
	print STDERR "Number of unique TreeLeaves =".scalar(@UniqTreeLeaveIDs)."\n" if($verbose);
	
	my $numberofnodeschanged = 0;
	
	foreach my $TreeLeafID (@TreeLeaveIDs){
	
			
		$sth->execute($TreeLeafID);
		my $nodename =	$TreeHash->{$TreeLeafID}{'node_id'};
			
		next unless($sth->rows > 0);
		#If SUPERFAMILY holds no record of the protein ID, simply move on to the next node id
		
		my @ProteinIDInfo;
		
		my ($ProteinDetail,$ProteinSecondaryDetails) = $sth->fetchrow_array();
			
		if($NameFlag){
			
			$ProteinDetail=~ m/<i[ ]*>(.*)<.i[ ]*>/i;
			$ProteinDetail = $1 if($1);
			$ProteinDetail = $ProteinDetail.'_'.$ProteinSecondaryDetails if($ProteinSecondaryDetails);
			$ProteinDetail = $ProteinDetail.'_'.$nodename;
			$ProteinDetail =~ s/ /_/g
		}
		
		
		
		push(@ProteinIDInfo,$ProteinDetail);
	
		my $String2ReplaceNodeName = join('-',@ProteinIDInfo);
	
		$TreeHash->{$TreeLeafID}{'node_id'}=$String2ReplaceNodeName;
		$dictionary->{$nodename} = $String2ReplaceNodeName;
	
		$numberofnodeschanged++;
		
	}
	
	$sth->finish;
	
	print STDERR "Number of nodes changed = ".$numberofnodeschanged."\n" if($verbose);
	
	my $BranchesFlag = 0;
	$BranchesFlag = 1 if($TreeHash->{$root}{'branch_length'});
	
	my $RenamedTree = ExtractNewickSubtree($TreeHash, $root,0,0);
	
	print STDOUT $RenamedTree."\n";	
	$numberoftreeschanged++;
}

print STDERR "Number of trees changed = ".$numberoftreeschanged."\n" if($verbose);
	
close TREE;

open FH,">ProteinID2ChangedValue.txt" or die $!."\t".$?;


while (my($key,$value)=each(%$dictionary)){
	
	print FH $key."\t".$value."\n";
}

close FH;


dbDisconnect($dbh) ; 

__END__

