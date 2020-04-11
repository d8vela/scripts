#!/usr/bin/perl

# Reads in FASTA of protein sequences and try to add GS-linker until desired length is reached

my $usage = "gslink_augment.pl <fasta_file.fasta> <desired_length>\n";

use strict;

my $file = $ARGV[0] or die $usage;
my $desired_len = $ARGV[1] or die $usage;

my $gs_link = 'GGGGS';

my ($count,$len,$fasta_ref) = fasta2hash($file);

my ($seq_id,$seq,$seq_len);
for $seq_id (sort {$a <=> $b} keys %{$fasta_ref}) {
	$seq = ${$fasta_ref}{$seq_id};
	$seq_len = length($seq);

	while ($seq_len < $desired_len) {
		$seq .= $gs_link;
		$seq_len = length($seq);
	}

	$seq_id =~ s/^\d+\_//;

	print ">$seq_id\n";
	print "$seq\n";
}

# ---- Subs ----

sub fasta2hash {
	
	my $msa_file = shift;

	my $count;
	my $id;
	my %fasta; 
	my $len;
	
	open(FILE,$msa_file) or die "$msa_file: $!";
	while (<FILE>) {
		tr/\r\n//d;
		if (/^>(\S+)/) { 
			$count++;
			$id = $1;
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


