#!/usr/bin/perl

# Written by David La
# Updated March 18, 2013

my $usage = "what_mutation.pl <wildtype_seq> <mutant_seq>\n";

my $seq1 = $ARGV[0] or die $usage;
my $seq2 = $ARGV[1] or die $usage;

# Get Pairwise Alignment
my ($seq1_aln,$seq2_aln) = align_pair($seq1,$seq2);

# Get Mutation Notations
my @muta = what_mutation($seq1_aln,$seq2_aln);

# Output Mutation Notations
my $notation;
for $notation (@muta) {
	print "$notation\n";
}

# ---- Subs ----

sub what_mutation {
	# Generate the X2Y mutation notation, given the original protein sequence and the mutated protein sequence.
	# Assumes that both protein sequences are the same length (no indels).
	my $protein_seq = shift;
	my $mutated_protein_seq = shift;

	# Check if original protein sequence is the same length as the mutated one
	my $protein_seq_length = length($protein_seq);
	my $mutated_protein_seq_length = length($mutated_protein_seq);
	die "ERROR: Protein Sequence is not the same length as the mutated protein sequence!\n" if $protein_seq_length != $mutated_protein_seq_length;
	
	# Find the protein sequence mutation
	my $aa;
	my $mutated_aa;
	my $index = 0;
	my $aa_pos;
	my @mut_notation;
	for $aa (split //, $protein_seq) {
		$mutated_aa = substr($mutated_protein_seq,$index,1);
		
		if ($aa ne $mutated_aa) {
			$aa_pos = $index + 1;
			push @mut_notation, "$aa$aa_pos$mutated_aa";
		}

		$index++;

	}

	return @mut_notation;
	
}

sub align_pair {
	# Use MUSCLE for Pairwise Sequence Alignment
	
	my $seq1 = shift;
	my $seq2 = shift;
	my $matrix = shift;
		
	my $muscle_exe = "/work/davela/bin/muscle";
	my $temp_dir = "/tmp/muscle_$$";

	mkdir($temp_dir);
	
	my $in_file = "$temp_dir/in.txt";
	my $out_file = "$temp_dir/out.txt";
	
	open(FILE,">$in_file") or die "Cannot open file $in_file\n";
	# Sequence 1
	print FILE ">seq1\n$seq1\n";
	# Sequence 2
	print FILE ">seq2\n$seq2\n";
	close(FILE);
	
	# Run MUSCLE
	if ($matrix) {
		# Use specified matrix
		`$muscle_exe -in $in_file -out $out_file -matrix $matrix 2> /dev/null`;
	}
	else {
		# Assume default matrix
		`$muscle_exe -in $in_file -out $out_file 2> /dev/null`;
	}
	
	# Read MUSCLE Output
	open(FILE,"$out_file") or die "Cannot open file $out_file\n";
	
	my $out;
	my $read;
	my ($first,$second);
	my ($aln_seq1,$aln_seq2);
	
	while (<FILE>) {
		chomp;
		if (/^>seq1/) {
			$first++;
		}
		elsif (/^>seq2/) {
			$first = 0;
			$second++;
		}
		elsif ($first) {
			$aln_seq1 .= $_;
		}
		elsif ($second) {
			$aln_seq2 .= $_;
		}
	}
	close(FILE);
	
	print `rm -r $temp_dir`;
	
	return ($aln_seq1,$aln_seq2);
	
}

