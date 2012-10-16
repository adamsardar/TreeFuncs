#!/usr/bin/env perl

=head1 NAME

prepareTraitTreeOverlay<.pl>

=head1 USAGE

 prepareTraitTreeOverlay.pl (-si --supraid |-da --domarch) Input the domain architecture to study, with by DA string or as a comb_id [*-ss --supresssupra] supress supra domain colouring -t --tree tree_file.nwk

EXAMPLE: prepareTraitTreeOverlay.pl -t ./EukaryoteTreeDomaarc.tree -si 34 -ss
#This will find all genomes possesing supra id 34 and them prepare a css template for use with nw_utils.

You can then plot this as an svg tree using the command:

nw_display -sr -S -w 3000 -c treedisplayoptions.css ./tree_file.nwk > tree.svg

You may need to play around with -w to get a width that looks 'good'

Domain architectures are shown in red (with infered evoltionary history up to their most recent common ancestor). Supra domains are in blue (and, since they are not modular evolutionary units, there is no attempt at ancestral lineage being show).
Genomes that exists as an outgroup to those possesing the da are shown in grey. Genomes beneath the LCA but not in possesion of the DA or a supra domain are in black.

=head1 SYNOPSIS

A script for use with with newick utils. Given a domain architecture or a supra_id, this script will prepare a css
file for use with nw_display. This will be a tree with all clades where a dom arch/supra_id is present are colored red
and all that aren't black.

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
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::Utils;
use Supfam::SQLFunc;
use Supfam::TreeFuncsNonBP;
use Supfam::DolloParsmony;
use Carp;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my $TreeFile;
my $SupraID;
my $DomArch;
my $SupressSupra; #Flag for whter to show supra domain assignmenets or not ...

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$TreeFile,
           "domarch|da=s" => \$DomArch,
           "supraid|si=s" => \$SupraID,
           "supresssupra|ss!" => \$SupressSupra,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------

open TREE, "<$TreeFile" or die $!.$?;
my $TreeString = <TREE>;
close TREE;

my ($root,$TreeCacheHash) = BuildTreeCacheHash($TreeString);

my @TreeLeaves = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};

my $Leaf2NodeIDDictionary = {};

map{$Leaf2NodeIDDictionary->{$TreeCacheHash->{$_}{'node_id'}}=$_}@TreeLeaves;

open CSS, ">treedisplayoptions.css" or croak $!.$?;

my $FullTreeString = join(' ',map{$TreeCacheHash->{$_}{'node_id'}}@TreeLeaves);
print CSS '"stroke:grey" Clade '.$FullTreeString."\n";

my $dbh = dbConnect();


my $sth;


#Create a list of all genomes in which trait is present as a full on domain architecture

if($SupraID){

	$sth = $dbh->prepare("SELECT distinct(genome) FROM len_supra WHERE ascomb_prot_number > 0 AND supra_id = ?;");
	$sth->execute($SupraID);

}elsif($DomArch){
		
	$sth = $dbh->prepare("SELECT distinct(genome) FROM len_supra JOIN comb_index ON len_supra.supra_id = comb_index.id WHERE  len_supra.ascomb_prot_number > 0 AND comb_index.comb = ?;");
	$sth->execute($DomArch);
	
}else{
	
	die "You must choose to colour tree either by domain architecture or supra_id!\n";
}


my @TraitGenomes; #A list of genomes possesing the trait to colour
while(my $genome = $sth->fetchrow_array()){push(@TraitGenomes,$genome);}

my (undef,$Intersection,undef,undef) = IntUnDiff(\@TraitGenomes,[map{$TreeCacheHash->{$_}{'node_id'}}@TreeLeaves]);		
@TraitGenomes = @$Intersection;
die "No overlap between SUPERFAMILY and tree genomes \n" unless(scalar(@TraitGenomes));

#Create a list of all genomes in which trait is present as a supra domain but not as an architecture

if($SupraID){

	$sth = $dbh->prepare("SELECT distinct(genome) FROM len_supra WHERE ascomb_prot_number = 0 AND number > 0 AND supra_id = ?;");
	$sth->execute($SupraID);

}elsif($DomArch){
		
	$sth = $dbh->prepare("SELECT distinct(genome) FROM len_supra JOIN comb_index ON len_supra.supra_id = comb_index.id WHERE  len_supra.ascomb_prot_number = 0 AND number > 0 AND comb_index.comb = ?;");
	$sth->execute($DomArch);
	
}else{
	
	die "You must choose to colour tree either by domain architecture or supra_id!\n";
}


my @SupraGenomes; #A list of genomes possesing the domain architecture as just a supra domain
while(my $genome = $sth->fetchrow_array()){push(@SupraGenomes,$genome);}

(undef,$Intersection,undef,undef) = IntUnDiff(\@SupraGenomes,[map{$TreeCacheHash->{$_}{'node_id'}}@TreeLeaves]);
my (undef,undef,undef,$GenomesWithOnlySupra) = IntUnDiff(\@TraitGenomes,$Intersection); #Find the genomes unique to the SupraGenomes list		


my $TraitMRCA = FindMRCA($TreeCacheHash,$root,[@{$Leaf2NodeIDDictionary}{@TraitGenomes}]); #Find the dollo parsimony MRCA
		
my @TraitFullCladeLeaves = @{$TreeCacheHash->{$TraitMRCA}{'Clade_Leaves'}};

my $TraitCladeString = join(' ',map{$TreeCacheHash->{$_}{'node_id'}}@TraitFullCladeLeaves);
print CSS '"stroke-width:2; stroke:red" Clade '.$TraitCladeString."\n";

dolloTraitDecoration($TreeCacheHash,$TraitMRCA,\*CSS,[@{$Leaf2NodeIDDictionary}{@TraitGenomes}],'"stroke-width:2; stroke:black" Clade ');

my $SupraGenoemsString = join(' ',@$GenomesWithOnlySupra);

unless($SupressSupra){
	
	print CSS '"stroke-width:2; stroke:blue" Individual '.$SupraGenoemsString."\n" if(scalar(@$GenomesWithOnlySupra)); #only print out the format string for supra domains if there are any in the first place
}

close CSS;

dbDisconnect($dbh);


__END__