#! /usr/bin/env perl

=head1 NAME

ExtractHomologousSequencesFromSuperfamilyI<.pl>

=head1 USAGE

ExtractHomologousSequencesFromSuperfamilyI.pl [options -v,-d,-h] <ARGS> (-f --family family_id|-sf --superfamily superfamily_id) [-c --cut flag_to_extract_assigned_sequences]

=head1 SYNOPSIS

A script to take as input a family or superfamily id and then return all sequences that have a model hitting to that classification within superfamily.
An optional flag allows you to cut only the sequences that have model assignments (so a protein of length 500 with hits at 34-212 and 367-457 will have a sequence of 34-457 from the original sequence)

This is a SLOW process (ca. 20min for a family)

=head1 AUTHOR

B<Adam Sardar> - I<adam.sardar@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2012 Gough Group, University of Bristol.

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
use Supfam::SQLFunc;
use Supfam::Utils;
use List::Util qw(max min);
use Time::HiRes;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__

my $output = 'HomologousSequences.fasta';
my $superfamily_id;
my $family_id;
my $cut_sequences = 0;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "help|h!" => \$help,
           "output|o:s" => \$output,
           "superfamily|sf=i" => \$superfamily_id,
           "family|f=i" => \$family_id,
           "cut|c!" => \$cut_sequences,
        ) or die "Fatal Error: Problem parsing command-line ".$!;


# Main Script Content
#----------------------------------------------------------------------------------------------------------------


my $dbh = dbConnect();


my $sth; 

my $tic = Time::HiRes::time;
print STDERR "Fetching sequences from SUPERFAMILY ...... ";

if($superfamily_id){

		$sth=$dbh->prepare("SELECT protein.seqid, genome_sequence.sequence, ass.region 
		FROM ass JOIN genome_sequence ON ass.protein = genome_sequence.protein
		JOIN protein ON ass.protein = protein.protein
		WHERE ass.sf = ?
		AND ass.evalue <= 0.0001
		;");
	$sth->execute($superfamily_id);
	
}elsif($family_id){
	
	$sth=$dbh->prepare("SELECT protein.seqid, genome_sequence.sequence, ass.region 
		FROM ass JOIN genome_sequence ON ass.protein = genome_sequence.protein 
		JOIN family ON family.auto = ass.auto
		JOIN protein ON ass.protein = protein.protein
		WHERE family.fa = ?
		AND family.evalue <= 0.0001
		;");

	$sth->execute($family_id);
	
}else{
	
	die "You must input a family or superfamily level id!\n";
}


my $toc = Time::HiRes::time;
my $Time = ($toc-$tic);

print STDERR "Done! ($Time seconds)";


my $SequencesHash = {};

while(my ($seqid,$sequence,$model_regions) = $sth->fetchrow_array()){
	
	$SequencesHash->{$seqid}=[$sequence,$model_regions];
}



if($cut_sequences){
	
	#If you only want regions of sequence that have hits to an HMM model.
	
	foreach my $seqid (keys(%$SequencesHash)){
		
		my $sequence = $SequencesHash->{$seqid}[0];
		
		my @Hit_indicies = map{split('-',$_)}split(',',$SequencesHash->{$seqid}[1]);
		my $MaxHit = max(@Hit_indicies);
		my $MinHit = min(@Hit_indicies);
		my $ExtractLength = ($MaxHit-$MinHit)+1;
		my $SeqStart = $MinHit-1; #Due to perl starting at 0
		
		
		if(length($sequence) < (length($ExtractLength) + length($SeqStart)) ){
			
			die "Something wrong with protein id: $seqid. Extract length: $ExtractLength with offset: $SeqStart. Length of sequence is:".length($sequence)." \n  "
			
		}
		
		my $ExtractedSeq = substr $sequence, $SeqStart, $ExtractLength;
		
		$SequencesHash->{$seqid}[0] = $ExtractedSeq;
	}
}

#EasyDump('ProtID2Sequence.dump',$SequencesHash);

open SEQOUT, ">$output" or die "Unable to open outfile $output".$!;

foreach my $seqid (keys(%$SequencesHash)){
	
	die "$seqid not found!\n" if($SequencesHash->{$seqid}[0] ~~ undef);
	
	print SEQOUT ">".$seqid."\n";
	print SEQOUT $SequencesHash->{$seqid}[0]."\n";
}

close SEQOUT;




__END__

