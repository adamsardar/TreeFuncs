#! /usr/bin/env perl

=head1 NAME

supfam2TraitMatrix<.pl>

=head1 USAGE

  supfam2TraitMatrix.pl [options -v,-d,-h] -t --tree <TreeFile in Newick> -o --output <outputfile name> -s --style <output style phylip|Hennig86|RAxML> (-T -- traitstyle supra|comb|multi  | -s --genomelist new line seperated lsit of genomes to produce traits for )

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
use Supfam::Utils;
use Supfam::TreeFuncsNonBP;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TreeFile;
my $GenomeListFile;
my $genome_archs_file;
my $outputfile = 'output';
my $OutputStyle = 'RAxML';
my $TraitStyle = 'comb';

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
	
	open OUT, ">$outputfile" or die $?;
		
	print OUT "$NoTreeTaxa\t$LineLength\n";
	
	foreach my $Taxon (@TreeTaxa){
	
		print OUT $Taxon."      "; #Phylip format is quite particular in the way it seperate tax names and state details.
		my $TraitString = $TraitHash->{$Taxon};
		print OUT $TraitString;
		print OUT "\n";
	}
	
	close OUT;
}

sub generateDomArchTraits($){
	
	my (@TreeTaxa);
	@TreeTaxa = @{$_[0]};
	
	my $NoTreeTaxa = scalar(@TreeTaxa);
	
	#Create a hash of all the trait vectors per taxon
	my $TraitHash = {};
	#$TraitHash -> {taxon => binary traits}
	my $FullSpeciesTraitsHash = {};
	#$TraitHash -> {taxon => binary traits}, but crucially, this will still contain sites which are identical throughout the whole sample of taxa
			
	my $lensupraquery = join ("' or len_supra.genome='", @TreeTaxa); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run - as there is no way to input a list of items explicitly into SQL, I'm just concatenating a string of truth statements
	
	my $dbh = dbConnect();
	my $sth = $dbh->prepare("SELECT DISTINCT(supra_id) FROM len_supra WHERE len_supra.ascomb_prot_number > 0 AND $lensupraquery;");
	$sth->execute();
	
	my @comb_ids;
	
	while (my $CombID = $sth->fetchrow_array() ) {
		
		push(@comb_ids,$CombID);
	}
	
	my %CombHash;
	
	@CombHash{@comb_ids}=((0)x scalar(@comb_ids));#Preallocate
	
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
		
	return($TraitHash);
		
}

sub generateSupraTraits($){
	
	my (@TreeTaxa);
	@TreeTaxa = @{$_[0]};
	
	my $NoTreeTaxa = scalar(@TreeTaxa);
	
	#Create a hash of all the trait vectors per taxon
	my $TraitHash = {};
	#$TraitHash -> {taxon => binary traits}
	my $FullSpeciesTraitsHash = {};
	#$TraitHash -> {taxon => binary traits}, but crucially, this will still contain sites which are identical throughout the whole sample of taxa
			
	my $lensupraquery = join ("' or len_supra.genome='", @TreeTaxa); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run - as there is no way to input a list of items explicitly into SQL, I'm just concatenating a string of truth statements
	
	my $dbh = dbConnect();
	my $sth = $dbh->prepare("SELECT DISTINCT(supra_id) FROM len_supra WHERE $lensupraquery;");
	$sth->execute();
	
	my @comb_ids;
	
	while (my $CombID = $sth->fetchrow_array() ) {
		
		push(@comb_ids,$CombID);
	}
	
	my %CombHash;
	
	@CombHash{@comb_ids}=((0)x scalar(@comb_ids));#Preallocate
	
	$sth = $dbh->prepare("SELECT supra_id FROM len_supra WHERE genome = ?;");
	
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
		
	return($TraitHash);
		
}

sub generateMultistateTraits($){
	
	#Multistate refers to a system which, rather than simply treating traits as present or ansent, gives several levels. 0 = absent, 1 = present as a supradomain and 2 = present as a domain architecture
	
	my (@TreeTaxa);
	@TreeTaxa = @{$_[0]};
	
	my $NoTreeTaxa = scalar(@TreeTaxa);
	
	#Create a hash of all the trait vectors per taxon
	my $TraitHash = {};
	#$TraitHash -> {taxon => binary traits}
	my $FullSpeciesTraitsHash = {};
	#$TraitHash -> {taxon => binary traits}, but crucially, this will still contain sites which are identical throughout the whole sample of taxa
			
	my $lensupraquery = join ("' or len_supra.genome='", @TreeTaxa); $lensupraquery = "(len_supra.genome='$lensupraquery')";# An ugly way to make the query run - as there is no way to input a list of items explicitly into SQL, I'm just concatenating a string of truth statements
	
	my $dbh = dbConnect();
	
	my $sth = $dbh->prepare("SELECT DISTINCT(supra_id) FROM len_supra WHERE $lensupraquery;");
	$sth->execute();
	
	my @supra_ids;
	
	while (my $supraID = $sth->fetchrow_array() ) {
		
		push(@supra_ids,$supraID);
	}
	#@supra_ids is a list of all supra ids in all taxa given
	
	my %SupraHash;
	@SupraHash{@supra_ids}=((0)x scalar(@supra_ids));#Preallocate
	
	my %ModelSupraHash = %SupraHash;
	my %CombHash = %SupraHash;
	
	foreach my $taxa (@TreeTaxa){
		
		$sth = $dbh->prepare_cached("SELECT supra_id FROM len_supra WHERE genome = ?;");
		$sth->execute($taxa);
		
		my %SpeciesSupraHash = %ModelSupraHash; #Create a duplicate of %SupraHash
				
		while (my $SpeciesSupraID = $sth->fetchrow_array() ) {
		
			$SpeciesSupraHash{$SpeciesSupraID}=1; #Per species supra domain pres/abs
			$SupraHash{$SpeciesSupraID}++; #Global total sightings
		}
		
		$sth = $dbh->prepare_cached("SELECT supra_id FROM len_supra WHERE ascomb_prot_number > 0 AND genome = ?;");
		$sth->execute($taxa);
		
		while (my $SpeciesCombID = $sth->fetchrow_array() ) {
		
			$SpeciesSupraHash{$SpeciesCombID}=2; #Per species supra domain pres/abs
			$CombHash{$SpeciesCombID}++; #Global total sightings of comb
			$SupraHash{$SpeciesCombID}--; #Correction to the supra hash (else it will be overcounted)
		}
			
		my @SpeciesFullSupras = @SpeciesSupraHash{sort(@supra_ids)}; #@SpeciesCombs is a hash slice , sorted by supra_id -> multstate  matrix string 0202210101 etc.
		$FullSpeciesTraitsHash->{$taxa}=join(',',@SpeciesFullSupras); #Keep as a string to remain effcient with memory
	}
	
	dbDisconnect($dbh) ; 
	
	#Calculate the informative sites and exclude the others
	my $index=0;
	my @InformativeSites;
	
	foreach my $supra_id (sort(@supra_ids)){
		{
	#		no warnings 'uninitialized'; #Stop perl complaining about unitialized hash entries
			push (@InformativeSites,$index) if(($CombHash{$supra_id} != $NoTreeTaxa) && ($SupraHash{$supra_id} != $NoTreeTaxa)  && ($CombHash{$supra_id} != 0));
			#i.e all or none of the taxa possess a comb or a supradomain, then exclude it from the trait hash
		}
		$index++;
	}
	
	#Selecting only the informative sites, create the trait strings which shall be outputted to file
	foreach my $taxa (@TreeTaxa){
		
		my @Traits = split(',',$FullSpeciesTraitsHash->{$taxa}); #All supras
		my $TraitString = join('',@Traits[@InformativeSites]); #Only the infromative supras
		$TraitHash->{$taxa}=$TraitString; #Join into string, so as to save on memory overhead
	}
		
	return($TraitHash);
		
}

#Main Script
#----------------------------------------------------------------------------------------------------------------


#Set command line flags and parameters

GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "tree|t=s" => \$TreeFile,
           "genomelist|g=s" => \$GenomeListFile,
           "output|o=s" => \$outputfile,
           "style|s:s" => \$OutputStyle,
           "traitstyle|T=s" => \$TraitStyle,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

my $NewickTrees = [];
my @TreeTaxa;

if($TreeFile){
	
	#Take as input a tree file. Extract the leaf names and push the newick strings onto $NewickTrees
	
	open FH, "<$TreeFile" or die $?;
	
	my $NewickStringOfTree = <FH>;
	
	my ($root,$TreeCacheHash) = BuildTreeCacheHash($NewickStringOfTree);
	print "Built TreeCacheHash\n";
	
	my @RootDescendents = @{$TreeCacheHash->{$root}{'all_Descendents'}};
	@TreeTaxa = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$root}{'Clade_Leaves'}};
		
	push(@$NewickTrees,$NewickStringOfTree);
	
	
	while(my $OtherTree = <FH>){
		
		next if($OtherTree =~ m/^\s*$/);
		
		chomp($OtherTree);
		push(@$NewickTrees,$OtherTree);
		
		my ($Altroot,$AltTreeCacheHash) = BuildTreeCacheHash($NewickStringOfTree);
		my @AltTreeTaxa = map{$AltTreeCacheHash->{$_}{'node_id'}}@{$AltTreeCacheHash->{$Altroot}{'Clade_Leaves'}};
		
		my (undef,undef,$ListAExclusive,$ListBExclusive) = IntUnDiff(\@TreeTaxa,\@AltTreeTaxa);
		die "Input tree taxa don't agree with one another\n" if(scalar(@$ListAExclusive) || scalar(@$ListBExclusive));
		#Error check the input trees
		
		undef($Altroot);
		undef($AltTreeCacheHash); #Hopeful lower the memory footprint
	}
	
	close FH;
	
	#TODO Rename leaf names in taxa for Phylip and push onto an array for output in the phylip file. currently the script just warns
	print STDERR "Phylip has a minimum taxon size\n" if(grep{length($_)< 3}@TreeTaxa  && $TraitStyle =~ m/phylip/i);	
	
}elsif($GenomeListFile){
	
	open FH, "<$GenomeListFile" or die $?;
	
	while (my $Genome = <FH>){
		
		chomp($Genome);
		push(@TreeTaxa,$Genome)
	}
	
	close FH;
	
}else{
	
	die "You need to specify a tree or list of genomes to study\n";
}


#error check for supfam codes

my $dbh = dbConnect();

my $sth = $dbh->prepare_cached("SELECT genome FROM genome WHERE genome = ?;");

foreach my $Taxa (@TreeTaxa){
	
	$sth->execute($Taxa);
	die "No entry in genome table of SUPEERFAMILY for $Taxa\n" unless($sth->rows());
	$sth->finish;
}

dbDisconnect($dbh) ; 


#Generate the appropriate set of traits

my $TraitHash;

if($TraitStyle =~ m/comb/i){
	
	$TraitHash = generateDomArchTraits(\@TreeTaxa);
	
}elsif($TraitStyle =~ m/supra/i){
	
	$TraitHash = generateSupraTraits(\@TreeTaxa);
	
}elsif($TraitStyle =~ m/multi/i){
	
	$TraitHash = generateMultistateTraits(\@TreeTaxa);
	
}else{
	
	$TraitHash = generateDomArchTraits(\@TreeTaxa);
	print STDERR "No Appropriate Output chosen, generating domain architecture traits instead \n";
}


#Wrtie only the records for species in the tree to file

if($OutputStyle =~ m/Hennig86/i){
	
	Hennig86Output($TraitHash,$outputfile);
	
}elsif($OutputStyle =~ m/Phylip/i){
	
	die "You must provide a tree input file if you want a Phylip output" unless(scalar(@$NewickTrees));
	PhylipOutput($TraitHash,$outputfile,$NewickTrees);
	
}elsif($OutputStyle =~ m/RAxML/i){
	
	RAxMLOutput($TraitHash,$outputfile);
	
}else{
	
	RAxMLOutput($TraitHash,$outputfile);
	print STDERR 'No Appropriate Output chosen, outputted RAxML format instead';
}

__END__

