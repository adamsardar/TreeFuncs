#! /usr/bin/perl -w

=head1 NAME

I<SupfamCodes2FullSpecies.pl.pl>

=head1 USAGE

 SupfamCodes2FullSpecies.pl.pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to take a supfam trait matrix and replace the superfamily ID codes with full taxon names using NCBI taxonomy (or leave them the same if there is no mapping)

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
use DBI;
use Supfam::SQLFunc;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $TraitsFiles; #Traits file in phylip format
my $output = 'ModifiedTraits.phy';

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "traitfile|tf=s" => \$TraitsFiles,
           "output|o:s" => \$output,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------


my $dbh = dbConnect();
my $sth = $dbh->prepare("SELECT ncbi_taxonomy.name FROM tree JOIN ncbi_taxonomy ON tree.taxon_id = ncbi_taxonomy.taxon_id WHERE tree.nodename=?");


open TRAITFILE, "<$TraitsFiles";
open MODIFIEDTRAITS, ">$output";

my $temp = <TRAITFILE>; # Skip past the firs (useless) line
print MODIFIEDTRAITS $temp;


while (my $line = <TRAITFILE>){
	
	
	my ($SFspeciesCode,$traits) = split(/\s+/,$line);

	$sth->execute($SFspeciesCode);	
	my $NCBISpecies = $sth->fetchrow_array();
	
	$NCBISpecies =~ s/\s+/\_/g;

	print MODIFIEDTRAITS $NCBISpecies."\t".$traits."\n";	
}

$sth->finish;

dbDisconnect($dbh) ; 

close TRAITFILE;
close MODIFIEDTRAITS;

__END__

