#! /usr/bin/env perl

=head1 NAME

checkAncestralSQLTable<.pl>

=head1 USAGE

checkAncestralSQLTable.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A simple script to check from errors in assignment of dollo parsimony ancestral states in SQL ancestral_len_supra

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
use DBI;

use Supfam::TreeFuncsNonBP;
use Supfam::Utils;
use Supfam::SQLFunc;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $rootleft;  #Root node from the tree (left id in SQL table) to be used.
my $rootright; #Root node from the tree (right id in SQL table) to be used.

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "rootleft|rl:i" => \$rootleft,
           "rootright|rr:i" => \$rootright,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------

my $dbh = dbConnect();
my $sth;

#Check that root is defined and extract left/right id

my $cladename;

if(defined($rootleft) && ! defined($rootright)){

	$sth = $dbh->prepare("SELECT tree.right_id,ncbi_taxonomy.name FROM tree JOIN ncbi_taxonomy ON ncbi_taxonomy.taxon_id = tree.taxon_id WHERE tree.left_id = ?;");
	$sth->execute($rootleft);
	($rootright,$cladename) = $sth->fetchrow_array();
	
}elsif(defined($rootright) && ! defined($rootleft)){
	
	$sth = $dbh->prepare("SELECT tree.left_id,ncbi_taxonomy.name FROM tree JOIN ncbi_taxonomy ON ncbi_taxonomy.taxon_id = tree.taxon_id WHERE tree.right_id = ?;");
	$sth->execute($rootright);
	($rootleft,$cladename) = $sth->fetchrow_array();
	
}else{
	
	die 'Need to provide either the left_id or right_id of root node, preferably both! ';	
}

print STDERR $rootleft." Root Left ".$rootright." Root Right    - Clade Name: ".$cladename."\n";

#Get all full domain architectures 

my $DomArch2leftRightIDs = {};
#This will be a hash of $hash->{DomArch}=[left_id,right_id]

$sth = $dbh->prepare("SELECT len_supra.supra_id, MIN(tree.left_id)
					FROM len_supra 
					JOIN genome 
					ON genome.genome = len_supra.genome
					JOIN tree
					ON genome.taxon_id = tree.taxon_id
					WHERE len_supra.ascomb_prot_number > 0
					AND (genome.include = 'y' OR genome.include = 's') 
					GROUP BY len_supra.supra_id;");
					
$sth->execute();

while (my ($DomArch,$Min_Left_id) = $sth->fetchrow_array()){
	
	$DomArch2leftRightIDs->{$DomArch}=[$Min_Left_id,undef];
}

$sth->finish;

$sth = $dbh->prepare("SELECT len_supra.supra_id, MAX(tree.right_id)
					FROM len_supra 
					JOIN genome 
					ON genome.genome = len_supra.genome
					JOIN tree
					ON genome.taxon_id = tree.taxon_id
					WHERE len_supra.ascomb_prot_number > 0
					AND (genome.include = 'y' OR genome.include = 's') 
					GROUP BY len_supra.supra_id;");

$sth->execute();

while (my ($DomArch,$Max_right_id) = $sth->fetchrow_array()){
	
	$DomArch2leftRightIDs->{$DomArch}[1]=$Max_right_id;
}

$sth->finish;

#Test that they were created at their MRCA

#Prepare SQL statement
$sth = $dbh->prepare_cached("SELECT ancestral_len_supra.ascomb_deleted,ancestral_len_supra.left_id,tree.right_id
					FROM ancestral_len_supra
					JOIN tree
					ON tree.left_id = ancestral_len_supra.left_id 
					WHERE ancestral_len_supra.left_id IN (SELECT MAX(left_id) FROM tree WHERE left_id < ? AND right_id > ?)
					AND ancestral_len_supra.supra_id = ?;");


open CREATED, ">./MissassignedCreations.txt" or die $!.$?;

my $DomArchMRCA={};

my $count = 0;

foreach my $DomArch (keys(%$DomArch2leftRightIDs)){
	
	my ($minleft,$maxright) = @{$DomArch2leftRightIDs->{$DomArch}};
	
	next if ($minleft == $maxright-1); #i.e. the architecture is present in just one genome, which is it's own MRCA
	

	$sth->execute($minleft,$maxright,$DomArch);

	my ($SupraStatus,$MRCAleft_id,$MRCAright_id) = $sth->fetchrow_array();
	
	if ($minleft <= $rootleft || $maxright >= $rootright ){
		delete($DomArchMRCA->{$DomArch});
		next ;
	}	#If the dom arch has a MRCA earlier than the desired root of study
	
	$DomArchMRCA->{$DomArch}=[$MRCAleft_id,$MRCAright_id];
	print STDERR "$DomArch not assigned correctly - min left_id: $minleft and max right: $maxright MRCA left_id: $MRCAleft_id MRCA status\n" if($SupraStatus ~~ undef);
	
	print CREATED "Domain Architecture with supra_id: $DomArch not assigned correctly - min left_id: $minleft and max right: $maxright MRCA left_id: $MRCAleft_id MRCA status: $SupraStatus (should be created!) \n" unless ($SupraStatus =~ m/created/i);	

	#last if ($count++ > 100);
}

close CREATED;

print STDERR "Checked all creation events\n";

#Test that they were deleted at the appropriate point


open DELETED, ">./MissassignedDeletions.txt" or die $!.$?;

#For each architecture
foreach my $DomArch (keys(%$DomArchMRCA)){
	
	my ($MRCAleft_id,$MRCAright_id) = @{$DomArchMRCA->{$DomArch}};
	#Find all the leaf nodes beneath the current MRCA that DONT have a dom arch assignment of the current supra id

	$sth = $dbh->prepare_cached("SELECT nodename,left_id,right_id 
								FROM tree 
								WHERE left_id >= ? 
								AND right_id <= ?
								AND left_id = right_id-1 
								AND nodename NOT IN (SELECT tree.nodename FROM tree JOIN genome ON tree.taxon_id = genome.taxon_id JOIN len_supra ON len_supra.genome=genome.genome WHERE tree.left_id >= ? AND tree.right_id <= ? AND len_supra.supra_id = ? AND genome.include = 'y' AND len_supra.ascomb_prot_number > 0);");	
	
	$sth->execute($MRCAleft_id,$MRCAright_id,$MRCAleft_id,$MRCAright_id,$DomArch);
	
	my $DomArchNonAssignments = {};
	
	while (my ($nodename,$lid,$rid) = $sth->fetchrow_array()){
		
		$DomArchNonAssignments->{$nodename}=[$lid,$rid];
	}
	
	$sth->finish();
		
	#Next if the dom arch is possesed by everything
	next if (scalar(keys(%$DomArchNonAssignments)) == 0);
	
	#Discard all the leaf nodes with dom arch assignments

	#For each unassigned leaf
	my @UnassignedGenomes = keys(%$DomArchNonAssignments);
				
		my $ProcessedGenomesHash = {};
				
		while (my $UnassignedGenome = pop(@UnassignedGenomes)){
			
			next if(exists($ProcessedGenomesHash->{$UnassignedGenome}));
			
			my ($genome_lid,$genome_rid) = @{$DomArchNonAssignments->{$UnassignedGenome}};
			
			my $indicator = 1;
			
				while($indicator){
				
				$sth = $dbh->prepare_cached("SELECT left_id,right_id FROM tree WHERE left_id+1 = ? OR right_id-1 = ?;");	
				#Find ascestor	
				$sth->execute($genome_lid,$genome_rid);
				
				die "Error here - SQL implies that there is more than one ancestor for left_id: $genome_lid right_id: $genome_rid\n" if($sth->rows != 1);
				
				my ($Ancestor_lid,$Ancestor_rid) = $sth->fetchrow_array();
				
				$sth->finish();
			
				#Test if any of the ancestors descendents posses the trait
				$sth = $dbh->prepare_cached("SELECT tree.nodename 
											FROM tree 
											JOIN genome ON tree.taxon_id = genome.taxon_id 
											JOIN len_supra ON len_supra.genome=genome.genome 
											WHERE tree.left_id >= ? 
											AND tree.right_id <= ? 
											AND len_supra.supra_id = ? 
											AND len_supra.ascomb_prot_number > 0;");	
			
				$sth->execute($Ancestor_lid,$Ancestor_rid,$DomArch);
				
				#If any descendents posses the trait
				if($sth->rows >= 1){
					
					my $ExampleOutGroup = $sth->fetchrow_array();
									
					$indicator = 0;
					$sth->finish();
					
					$sth = $dbh->prepare_cached("SELECT nodename,left_id,right_id 
								FROM tree 
								WHERE left_id <= ? 
								AND right_id >= ?
								AND left_id = right_id-1 
								AND nodename NOT IN (SELECT tree.nodename FROM tree JOIN genome ON tree.taxon_id = genome.taxon_id JOIN len_supra ON len_supra.genome=genome.genome WHERE tree.left_id >= ? AND tree.right_id <= ? AND len_supra.supra_id = ? AND genome.include = 'y' AND len_supra.ascomb_prot_number > 0);");					
					#Same query as earlier, but this time we're just going to get all of the genomes that are beneath the current node (so that we can remove them from @UnassignedGenomes)
					
					$sth->execute($genome_lid,$genome_rid,$genome_lid,$genome_rid,$DomArch);
					
					my $TempArrayOfDeletedGenome = [];
					
					while(my $UnassGenome = $sth->fetchrow_array()){		
						
						$ProcessedGenomesHash->{$UnassGenome}=undef;
					}
					
					$sth->finish();
					
					#Check for deletion at point in the tree:
					
					$sth = $dbh->prepare("SELECT ancestral_len_supra.ascomb_deleted,ancestral_len_supra.left_id,tree.right_id
							FROM ancestral_len_supra
							JOIN tree
							ON tree.left_id = ancestral_len_supra.left_id
							WHERE ancestral_len_supra.left_id = ?
							AND ancestral_len_supra.supra_id = ?;");
					
					$sth->execute($genome_lid,$DomArch);
							
					my ($SupraStatus,$del_lid,$del_rid) = $sth->fetchrow_array();			
					
					print STDERR "Domain Architecture with supra_id: $DomArch not assigned correctly - left_id: $del_lid and right_id: $del_rid" if($SupraStatus ~~ undef);
	
					
					print DELETED "Domain Architecture with supra_id: $DomArch not assigned correctly - left_id: $del_lid and right_id: $del_rid: MRCA status: $SupraStatus (should be deleted!) - Example OutGroup: $ExampleOutGroup \n" unless ($SupraStatus =~ m/deleted/i);	
										
					
				}else{
					
					($genome_lid,$genome_rid)=($Ancestor_lid,$Ancestor_rid);#Move up the tree
				}
				
				
			#If yes, then a deletion has occured along that branch
			#If no, move up to the next ancestor, repeat
			#Die with an error if you hit the MRCA 
			
			#Once you find a deletion, see whehter you agree with the ance SQL assignment
			# Dump to a file if you don't
			
				$sth->finish();
				
			}
		}
}

close DELETED;

print STDERR "Checked all deletion events\n";
__END__