#! /usr/bin/env perl

=head1 NAME

TreeSQL2Newick<.pl>

=head1 USAGE

TreeSQL2Newick.pl [options -v,-d,-h] <ARGS> -rl|--rootleft node_left_id -rr|--rootright node_right_id 

=head1 SYNOPSIS

A script to parse an sql table from SUPERFAMILY and output a tree in newick format. Simply specifiy the root of the tree and this script will
sort out the rest. Specify the root by passing in eith the left or the right id of the root node, or better still both.

Outputs tree as a newick file.

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
use Supfam::TreeFuncsNonBP;
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
           "rootleft|rl:s" => \$rootleft,
           "rootright|rr:s" => \$rootright,
        ) or die "Fatal Error: Problem parsing command-line ".$!;
        
# Main Script Content
#----------------------------------------------------------------------------------------------------------------

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

print STDERR $rootleft." Root Left ".$rootright." Root Right\n";

my $SQLTreeCacheHash = supfamSQL2TreeHash($rootleft,$rootright);

my $FullDeletionsNormailsedTree = ExtractNewickSubtree($SQLTreeCacheHash,$rootleft,1,0);

print $FullDeletionsNormailsedTree."\n";

__END__

