#! /usr/bin/env perl

=head1 NAME

supfam2Phylip<.pl>

=head1 USAGE

  supfam2TraitMatrix.pl [options -v,-d,-h] -t --tree <TreeFile in Newick> -o --output <outputfile name> -s --style <output style phylip|Hennig86|RAxML>

=head1 SYNOPSIS

A script to generate a carachter trait matrix file of domain architecture combinations copatible with Phylip|RAxML|Hennig86 format from a specified tree file 
(or list of trees with the same leaves - Phylip). This works purely with domaind combinations, not supradomains (yet!).

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

02-August-2011 Initial Entry

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "$ENV{HOME}/bin/perl-libs-custom/";

#CPAN Includes
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
my $OutputStyle = 'Phylip';
my $cachecombs = 0;

# Sub definitions
#----------------------------------------------------------------------------------------------------------------

sub Hennig86Output($$){
	
	my ($TraitHash,$outputfile) = @_;
	
	open OUT, ">$outputfile";
		
	my @TreeTaxa = keys(%$TraitHash);
	my $NoTreeTaxa = scalar(@TreeTaxa); 
	my $LineLength = length($TraitHash->{$TreeTaxa[0]});
		
	print OUT "xread\n'This is a carachter trait file in Hennig86 format, using domain architectures (or supra domains, depending on the input file). This file was created using supfam2Phyliptraits'\n";
	print OUT "$LineLength $NoTreeTaxa\n";
	
	foreach my $Taxon (keys(%$TraitHash)){
		
		print OUT $Taxon;
		print OUT "\t";
		print OUT $TraitHash->{$Taxon};
		print OUT "\n";
	}
	close OUT;
}

sub PhylipOutput($$$){

	my ($TraitHash,$outputfile,$NewickTrees) = @_;
	
	my @TreeTaxa = keys(%$TraitHash);
	my $NoTreeTaxa = scalar(@TreeTaxa); 
	my $LineLength = length($TraitHash->{$TreeTaxa[0]});
	
	open OUT, ">$outputfile";
		
	print OUT "$NoTreeTaxa\t$LineLength\n";
	
	foreach my $Taxon (@TreeTaxa){
	
		print OUT "Taxon".$Taxon."      "; #Phylip sucks SO much that the minimu taxon name length is five carachters. URGH!
		my $TraitString = $TraitHash->{$Taxon};
		print OUT $TraitString;
		print OUT "\n";
	}
	
	print OUT scalar(@$NewickTrees)."\n";
	
	open OUTTREE, ">intree"; #This is purely for PHYLIP. I feel dirty writing this stupid step.
	
	foreach my $Tree (@$NewickTrees){
			
		$Tree =~ s/:[\.\d\-\+E]*//g;
		$Tree =~ s/(\))([\.\d\-\+E]+)(;$)/$1$3/;
		$Tree =~ s/\((Taxon\w{2,3})\)/$1/g; #Frustratingly, bioperl will do things liek output a clade like this (A,B,(C)) rather than (A,B,C)
		
		#$Tree =~ s/^\(//;
		#$Tree =~ s/\);$/;/; #Phylip complains that there is one too many parenthases when bioperl outputs the tree. Go figure
		
		print OUT $Tree."\n";
		print OUTTREE $Tree."\n";
	}
	
	close OUTTREE;
	close OUT;
}

sub RAxMLOutput($$){

	my ($TraitHash,$outputfile) = @_;
	
	my @TreeTaxa = keys(%$TraitHash);
	my $NoTreeTaxa = scalar(@TreeTaxa); 
	my $LineLength = length($TraitHash->{$TreeTaxa[0]});
	
	open OUT, ">$outputfile";
		
	print OUT "$NoTreeTaxa\t$LineLength\n";
	
	foreach my $Taxon (@TreeTaxa){
	
		print OUT $Taxon."      "; #Phylip sucks SO much that the minimu taxon name length is five carachters. URGH!
		my $TraitString = $TraitHash->{$Taxon};
		print OUT $TraitString;
		print OUT "\n";
	}
	
	close OUT;
}


#Main Script
#----------------------------------------------------------------------------------------------------------------


#Set command line flags and parameters.n
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$TreeFile,
           "output|o=s" => \$outputfile,
           "style|s:s" => \$OutputStyle,
           "cachecombs|c:i" => \$cachecombs,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Get a collection of nodeIDs by which to root tree

my $input = new Bio::TreeIO(-file   => "$TreeFile",
                            -format => "newick") or die $!;
                            
my $tree = $input->next_tree;
my $root = $tree->get_root_node;
my @RootDescendents = $root->get_all_Descendents;

#Extract all the species in the treefile given
my @TreeTaxa = map{$_->id}grep{$_->is_Leaf}@RootDescendents;
my $NoTreeTaxa = scalar(@TreeTaxa);

#Check if there are any other trees in the input file and add them to an array of strings. Also, ensure that they have the same taxa

my $NewickTrees = []; #A list of newick trees from input file
map{my $renamenode = $tree->find_node(-id => $_); $renamenode->id ("Taxon".$_);}@TreeTaxa if($OutputStyle =~ m/Phylip/i); # Phylip is so shit that it has a minimum taxon label of 5 carachters. Urgh. I muck around with it here to make the output tree fit
my $NewickTree = $tree->as_text('newick');

push (@$NewickTrees,$NewickTree);

while (my $OtherTree = $input->next_tree){ #Alt is a prefix to show that it refers to another tree from the same input file
	
	unless (scalar($OtherTree->get_leaf_nodes)){
		
		my $Altroot = $tree->get_root_node;
		my @AltTreeTaxa = map{$_->id}grep{$_->is_Leaf}($Altroot->get_all_Descendents);
		
		map{my $renamenode = $OtherTree->find_node(-id => $_); $renamenode->id("Taxon".$_);}@AltTreeTaxa if($OutputStyle =~ m/Phylip/i);
		
		my $AltNewickTree = $OtherTree->as_text('newick');
		push (@$NewickTrees,$AltNewickTree) unless ($AltNewickTree =~ m/^;/);
		
		my ($Union,$Intersection,$ListAExclusive,$ListBExclusive) = IntUnDiff(\@AltTreeTaxa,\@TreeTaxa);	#Test to see if there are any difference between the two lists of taxa
		
		if (scalar(@$ListAExclusive)){
			
			print STDERR 'Master Tree -> '.$NewickTree."\n";
			print STDERR 'Current Tree -> '.$AltNewickTree."\n";
			die "Different genomes in two trees of input file (see above) \n";		
		}
	}
}# Quick while loop to gather all the other trees in the input file and check that they posses the same genomes. These shall be written to file in the Phylip (default) use case.



#Create a hash of all the trait vectors per taxon
my $TraitHash = {};
#$TraitHash -> {taxon => binary traits}
my $FullSpeciesTraitsHash = {};
#$FullSpeciesTraitsHash -> {taxon => [binary prescnce/absecene data]}

my $lensupraquery = join ("' or len_supra.genome='", @TreeTaxa); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run - as there is no way to input a list of items explicitly into SQL, I'm just concatenating a string of truth statements

my $dbh = dbConnect();
my $sth = $dbh->prepare("SELECT DISTINCT(supra_id) FROM len_supra WHERE len_supra.ascomb_prot_number > 0 AND $lensupraquery;");
$sth->execute();

my @comb_ids;


#unless(-e './.Combs' && $cachecombs){ #Script allows for caching the comb list - or indeed feeding in a strangely tailored combs list
	while (my $CombID = $sth->fetchrow_array() ) {
	
		push(@comb_ids,$CombID);
	}
	#EasyDump('.Combs',\@comb_ids);
#}else{
	
#	my $combspoint = EasyUnDump('.Combs');
	#@comb_ids = @{$combspoint};
#}



my %CombHash;

@{\%CombHash}{@comb_ids}=((0)x scalar(@comb_ids));

$sth = $dbh->prepare("SELECT supra_id FROM len_supra WHERE ascomb_prot_number > 0 AND genome = ?;");

my %ModelCombHash = %CombHash;



foreach my $taxa (@TreeTaxa){
	
	my %SpeciesCombsHash = %ModelCombHash; #Create a duplicate of %CombHash
	
	$sth->execute($taxa);
	
	while (my $SpeciesCombID = $sth->fetchrow_array() ) {
	
		$SpeciesCombsHash{$SpeciesCombID}=1; #Per species presence/abscence
		$CombHash{$SpeciesCombID}++; #Global total sightings
	}
		
	my @SpeciesCombs = @SpeciesCombsHash{sort(@comb_ids)}; #Sorted by comb_id -> presences absece matrix 000101 etc
	
	$FullSpeciesTraitsHash->{$taxa}=join(',',@SpeciesCombs);
}

dbDisconnect($dbh) ; 

#Calculate the informative sites and exclude the others
my $index=0;
my @InformativeSites;

foreach my $comb_id (sort(@comb_ids)){
	
	push (@InformativeSites,$index) if($CombHash{$comb_id} != $NoTreeTaxa && $CombHash{$comb_id} != 0);
	$index++;
}

#Selecting only the informative sites, create the trait strings which shall be outputted to file
foreach my $taxa (@TreeTaxa){
	
	my @Traits = split(',',$FullSpeciesTraitsHash->{$taxa}); #Full combs
	my $TraitString = join('',@Traits[@InformativeSites]);
	$TraitHash->{$taxa}=$TraitString;
}

#Wrtie only the records for species in the tree to file

if($OutputStyle =~ m/Hennig86/i){
	
	Hennig86Output($TraitHash,$outputfile);
	
}elsif($OutputStyle =~ m/Phylip/i){
	
	PhylipOutput($TraitHash,$outputfile,$NewickTrees);
	
}elsif($OutputStyle =~ m/RAxML/i){
	
	RAxMLOutput($TraitHash,$outputfile);
	
}else{
	
	PhylipOutput($TraitHash,$outputfile,$NewickTrees);
	print STDERR 'No Appropriate Output chosen, outputted Phylip format instead';
}

__END__

