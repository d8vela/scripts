#!/usr/bin/perl

# Written by David La
# Updated: Thu Dec  4 17:17:51 PST 2014

# Description:
# Reads in FASTA and outputs sequences of the desired length range.

my $usage = "Usage: lengthFASTAFilter.pl <file.fasta> <min_length> <max_length>\n";

use strict;

my $file = $ARGV[0] or die $usage;
my $min_len = $ARGV[1] or die $usage;
my $max_len = $ARGV[2] or die $usage;

my ($count,$len,$fasta_ref) = fasta2hash($file);

my ($seq_id,$seq,$seq_len);
for $seq_id (sort {$a <=> $b} keys %{$fasta_ref}) {
	$seq = ${$fasta_ref}{$seq_id};
	$seq_len = length($seq);
	if ($seq_len >= $min_len and $seq_len <= $max_len) {
		print ">$seq_id | length=$seq_len\n$seq\n";
	}
}

# ---- Subs ----

sub fasta2seq {
	
	my $fasta_ref = shift;
	my $seqs;

	my $id;
	for (sort {$a <=> $b} keys %{$fasta_ref}) {
		$id = $_;
		$id =~ s/^\d+_(.*)/$1/g;
		$seqs .= "$id\t";
		$seqs .= "${$fasta_ref}{$_}\n";
	}

	return $seqs;

}

sub fasta2hash {
	
	my $msa_file = shift;
	
	my $count;
	my $id;
	my %fasta; 
	my $len;
	#my $buffer_len;
	
	open(FILE,$msa_file) or die "$msa_file: $!";
	while (<FILE>) {
		tr/\r\n//d;
		if (/^>(\S+)/) { 
			$count++;
			$id = $1;
			#die "Not all MSA lengths are equal!\n" if $buffer_len && $len != $buffer_len;
			#$buffer_len = $len;
			$len = 0;
		}
		else { 
			$fasta{"$count\_$id"} .= $_;
			$len += length($_);
		}
	}
	close(FILE);
	
	return ($count,$len,\%fasta);

}


