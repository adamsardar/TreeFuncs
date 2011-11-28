#! /usr/bin/perl -w

=head1 NAME

AncestralSQL2Newick<.pl>

=head1 USAGE

 AncestralSQL2Newick.pl [options -v,-d,-h] <ARGS>
 
 Example usage:
 
 AnceSQL2Newick.pl -rl 3 -rr 1364 -o Ancestral
 
 To extract genomes from the root of all eukaryotes.
 Don't forget! Bacterial dollo parsmiony results are not present in SUPERFAMILY, as dollo parsimony is a very poor model for bacteria.

=head1 SYNOPSIS

A script to parse an the ancestral_info table from SUPERFAMILY and output a tree with bracnh lengths proportional to deletions in newick format. Simply specifiy the root of the tree and this script will
sort out the rest.

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "/home/sardar/bin/perl-libs-custom/";


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
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use IO::String;
use Supfam::SQLFunc;
# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $rootleft;  #Root node from the tree (left id in SQL table) to be used.
my $rootright; #Root node from the tree (right id in SQL table) to be used.
my $deletiontreeflag; #Produce a tree with branch lengths equal to the number of domain architecture deltions along that edge
my @excluded; # A list of genomes not wanted in this tree
my $fileout = 'fileout'; #Output file
my $include = 1; #use the 'include=y' option from the genome table?

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "rootleft|rl=i" => \$rootleft,
           "rootright|rr=i" => \$rootright,
           "deletions|del:i" => \$deletiontreeflag,
           "output|o=s" => \$fileout,
           "exclude|ex:s" => \@excluded,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

# Main Script Content
#----------------------------------------------------------------------------------------------------------------


my $RightTreeHash = {}; #A hash of the tree ordered by right id
my $LeftTreeHash = {}; #A hash of the tree ordered by left id

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

$sth = $dbh->prepare("SELECT tree.left_id,tree.right_id,ancestral_info.comb_deleted,tree.nodename FROM tree JOIN ancestral_info ON tree.left_id=ancestral_info.left_id WHERE tree.left_id >= ? AND right_id <= ?;");
$sth->execute($rootleft,$rootright);

#Initialise two hashes, keyed respectively by right and left ids, containing node data.
while (my ($leftid,$rightid,$edge_length,$nodename) = $sth->fetchrow_array()){
	
	$edge_length=($edge_length+0.1)/33000; #33000 is a normalising factor
	
	my $NewickClade = 'NULL';#The calde below this node in newick format - initialise with a default value
	
	my $NodeData = [$nodename,$edge_length,$NewickClade];
	
	$RightTreeHash->{$rightid} = [$leftid,$NodeData];
	$LeftTreeHash->{$leftid} = [$rightid,$NodeData];
}

print join('-',sort(keys(%$LeftTreeHash)));
print "\n";

my $NextGenLeftIDs = {};
my $CurrentGenLeftIDs = [];

foreach my $leftid (keys(%$LeftTreeHash)){
	
	my $rightid = $LeftTreeHash->{$leftid}[0];
	
	if ($rightid == ($leftid+1)){ #if node is a leaf
			
		no warnings 'uninitialized';	
			
		my ($nodename,$edge_length,$NewickClade) = @{$LeftTreeHash->{$leftid}[1]};

		my $NewickString = $nodename.":".$edge_length;
		$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string from 'NULL' to correct value
		
		my $ParentLeftID;
		$ParentLeftID = ($leftid-1) if(exists($LeftTreeHash->{($leftid-1)}));
		$ParentLeftID = $RightTreeHash->{$rightid+1}[0] if(exists($RightTreeHash->{$RightTreeHash+1}));
		$NextGenLeftIDs->{$ParentLeftID}++;
	}	
}
#Collect a list of all leaf nodes (possibility of adding  additional functionality here - exclude a list of genomes)

print $rootleft."  <- root left  ".$rootright."  <- root right\n";

my $count =0;

while ($LeftTreeHash->{$rootleft}[1][2] eq 'NULL'){
	
	#print "Iteration: ".$count++;
	
	no warnings;
	
	@$CurrentGenLeftIDs = keys(%$NextGenLeftIDs);
	$NextGenLeftIDs = {};
	#print "     Nodes in curent gen: ".scalar(@$CurrentGenLeftIDs)."\n";
	
	foreach my $leftid (@$CurrentGenLeftIDs) {
		
		my $rightid = $LeftTreeHash->{$leftid}[0];
		
		my ($LeftNewick,$RightNewick) = ('NULL','NULL');
		$LeftNewick = $LeftTreeHash->{$leftid+1}[1][2] if (exists($LeftTreeHash->{$leftid+1}));
		$RightNewick = $RightTreeHash->{$rightid-1}[1][2] if (exists($RightTreeHash->{$rightid-1}));
								
		if($LeftNewick ne 'NULL' && $RightNewick ne 'NULL'){
			
			#Aggregate two nodes below this one into a newick format
			my $edge_length = $LeftTreeHash->{$leftid}[1][1];
			my $NewickString = "(".$LeftNewick.",".$RightNewick."):".$edge_length;
			#print $NewickString."\n";
			$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string
											
		}else{

			$NextGenLeftIDs->{$leftid}++;
		}
				
		if($LeftTreeHash->{$leftid}[1][2] ne 'NULL'){
		
		$LeftTreeHash->{($leftid+1)}[1][2] = '' if(exists($LeftTreeHash->{($leftid+1)}));
		$RightTreeHash->{$rightid-1}[1][2] = '' if(exists($RightTreeHash->{$rightid-1}));
		
		#Boring way to find the parent of the node under study if we have already found the newick format representation of the clade below (parent is ether left_id-1 or right_id+1)			
			my $ParentLeftID;
			$ParentLeftID = ($leftid-1) if(exists($LeftTreeHash->{($leftid-1)}));
			$ParentLeftID = $RightTreeHash->{$rightid+1}[0] if(exists($RightTreeHash->{$RightTreeHash+1}));
			$NextGenLeftIDs->{$ParentLeftID}++;
		}	
	}	
}


#Strange way to output the newick format tree - this is simply so as to be sure that the outputted tree is bioperl compatible
my $FullTree = $LeftTreeHash->{$rootleft}[1][2];
my $TreeString = "(".$FullTree.");\n";

my $io = IO::String->new($TreeString);

my $input = Bio::TreeIO->new(-fh => $io,
                              -format => 'newick') or die $!;

my $tree = $input->next_tree;

my $NewickTree = $tree->as_text('newick');

my $treeout = new Bio::TreeIO(-format => 'newick', -file   => ">$fileout");
$treeout->write_tree($tree);

__END__

