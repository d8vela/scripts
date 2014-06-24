#!/usr/bin/perl

# Written by David La
# Updated Wed May 21 01:49:29 PDT 2014

# Description:
# This script will take an input SSM data (*.tab) and calculate
# the sequence conservation and map them on to the B-factors of
# a PDB file of your design

my $usage = "ssm2pdb.pl <ssm_file.tab> <pdb_file.pdb> [pdb_chain]\n";

use strict;

my $ssm_file = $ARGV[0] or die $usage;
my $pdb_file = $ARGV[1] or die $usage;
my $pdb_chain = $ARGV[2] || "B"; # Default for protein interface design

# ---- Read PDB File ----
my %pdb_atom = readPDBAtom($pdb_file,$pdb_chain);

# ---- Get the Design PDB Sequence ----
my $design_pdb_seq = atom2seq(\%pdb_atom);

# ---- Get the Design SSM Sequence (possibly with residues of missing densities) ----
my $design_ssm_seq = ssm2seq($ssm_file);


# ---- Align PDB to SSM WT Sequence ----
my ($design_pdb_seq_aln,$design_ssm_seq_aln) = align_pair($design_pdb_seq,$design_ssm_seq);

my %design_aln = ('PDB' => $design_pdb_seq_aln,
				  'SSM' => $design_ssm_seq_aln);

# ---- Reference PDB and Alignment Positions to Real Residue Positions ----

# Hash for Referencing PDB Residue Number to Real Residue Number
my %aa_num = pdbSeqNum(\%pdb_atom);

# Hash for Referencing Real Residue Number to Alignment Position Number
my %aln_num = alnSeqNum(\%design_aln,$design_pdb_seq);

# ---- Calculate Sequence Conservation Score ----

#my %score_ssm = ssm_entropy($ssm_file);
my %score_ssm = ssm_ratio($ssm_file);

# Convert Alignment Positions to PDB Positions
my %aln2aa = reverse %aln_num;
my %aa2pdb = reverse %aa_num;

my %score_pdb;
my $aa_pos;
my $pdb_pos;
for my $pos (keys %score_ssm) {
	$aa_pos = $aln2aa{$pos}; # Get Real AA position
	$pdb_pos = $aa2pdb{$aa_pos}; # Get PDB AA position
	$score_pdb{$pdb_pos} = $score_ssm{$pos};
}

# ---- Output PDB with Scores Mapped to B-Factors ----

#my $map_out = map2bfactor($pdb_file,$pdb_chain,\%score_pdb,1,365,420);
my $map_out = map2bfactor($pdb_file,$pdb_chain,\%score_pdb,10);
print $map_out;


# ---- Subs ----

sub readPDBAtom {
	
	my $pdb_file = shift;
	my $pdb_chain = shift;

	# -- Read PDB File

	my ($index,$atom,$aa,$mer,$res_num);
	my ($x,$y,$z);
	my %pdb_atom;

	open(PDB,$pdb_file) or die "Cannot open PDB file: $pdb_file\n";
	while (<PDB>) {
		# Skip non-atom lines
		next unless /^ATOM/;

		$index = substr($_,6,5);
		$index =~ s/ //g;

		# Parse necessary PDB data
		$atom = substr($_,13,3);
		$atom =~ s/ //g;
		
		$aa = substr($_,17,3);
		$mer = substr($_,21,1);

		# Skip chain if not needed unless it's not available (NA)
		unless ($pdb_chain eq 'NA') {
			next if $mer ne $pdb_chain;
		}
		
		$res_num = substr($_,22,8);
		$res_num =~ s/ //g;
		
		$x = substr($_,31,7);
		$x =~ s/ //g;
		
		$y = substr($_,39,7);
		$y =~ s/ //g;
		
		$z = substr($_,47,7);
		$z =~ s/ //g;
		
		# Store PDB Coordinates
		$pdb_atom{$index}{'atom'} = $atom;
		$pdb_atom{$index}{'aa'} = $aa;
		$pdb_atom{$index}{'res_num'} = $res_num;
		$pdb_atom{$index}{'x'} = $x;
		$pdb_atom{$index}{'y'} = $y;
		$pdb_atom{$index}{'z'} = $z;
		
	}
	close(PDB);

	return %pdb_atom;
	
}

sub atom2seq {
	
	my $pdb_atom_ref = shift;

	my $index;
	my $seq;
	my $res;
	my $aa;
	my $res_num;
	my $res_num_rem;

	my %aa_code = (
		'ALA' => 'A',
		'CYS' => 'C',
		'ASP' => 'D',
		'GLU' => 'E',
		'PHE' => 'F',
		'GLY' => 'G',
		'HIS' => 'H',
		'ILE' => 'I',
		'LYS' => 'K',
		'LEU' => 'L',
		'MET' => 'M',
		'ASN' => 'N',
		'PRO' => 'P',
		'GLN' => 'Q',
		'ARG' => 'R',
		'SER' => 'S',
		'THR' => 'T',
		'VAL' => 'V',
		'TRP' => 'W',
		'TYR' => 'Y'
	);
	
	for $index (sort {$a <=> $b} keys %{$pdb_atom_ref}) {
		$res_num = ${$pdb_atom_ref}{$index}{'res_num'};
		if ($res_num_rem ne $res_num) {
			$res = ${$pdb_atom_ref}{$index}{'aa'};
			$aa = $aa_code{$res};
			$seq .= $aa;
		}
		$res_num_rem = $res_num;
	}

	return $seq;
	
}

sub pdbSeqNum {
	
	my $pdb_atom_ref = shift;

	my $index;
	my $res_num;
	my $res_num_rem;
	my $count;
	my %aa_num;
	
	for $index (sort {$a <=> $b} keys %{$pdb_atom_ref}) {
		$res_num = ${$pdb_atom_ref}{$index}{'res_num'};
		if ($res_num_rem != $res_num) {
			$count++;
			$aa_num{$res_num} = $count;
		}
		$res_num_rem = $res_num;
	}

	return %aa_num;
	
}

sub ssm2seq {
	
	my $ssm_file = shift;
	
	open(FILE,"$ssm_file") or die "Cannot open file $ssm_file\n";

	my @line;
	my $aa;
	my $wt;
	my $seq;
	while (<FILE>) {
		chomp;
		
		# print;
		@line = split /\t/, $_;
		$wt = $line[8];
		if ($wt eq 'WT') {
			$aa = $line[2];
			$seq .= $aa;
		}
	}
	close(FILE);
	
	return $seq;
	
}

sub alnSeqNum {
	
	my $aln_ref = shift;
	my $seq = shift;

	my $seq_id;
	my $found_id;
	my $aln_seq;
	my $aln_seq_degap;
	my $aln_aa;
	my $aln_count;
	my $seq_count;
	my %aln_num;
	
	for $seq_id (sort {$a <=> $b} keys %{$aln_ref}) {
		$aln_seq = ${$aln_ref}{$seq_id};
		$aln_seq_degap = $aln_seq;
		$aln_seq_degap =~ s/-//g;
		if ($seq eq $aln_seq_degap) {
			for $aln_aa (split //, $aln_seq) {
				$aln_count++;
				unless ($aln_aa eq '-') {
					$seq_count++;
					$aln_num{$seq_count} = $aln_count;
				}
			}
			last;
		}
	}
	
	die "Cannot find PDB sequence ($seq) in aln\n" if !$aln_count;

	return %aln_num;
	
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

sub ssm_entropy {
	# Calculate the sequence entropy of every residue position in the SSM "tab" input file
	
	my $ssm_file = shift;
	my $mutation_th = shift || 1;
	
	open(FILE,"$ssm_file") or die "Cannot open file $ssm_file\n";

	my %a_freq; # Enrichment (Sorts)
	my %b_freq; # Background (Naive Library)
	
	my $pos;
	my $aa;
	my @line;
	
	while (<FILE>) {
		
		next if $. == 1;
		chomp;
		
		@line = split "\t", $_;
		
		# Skip WT amino acids
		#next if $line[8] eq 'WT';
		
		# Amino acid sequence position
		$pos = $line[1];
		# Start with 1
		$pos++;
		
		$aa = $line[2];
		$a_freq{$pos}{$aa} = $line[3]; # Total Mutation Frequency for Enrichment Sorts
		$b_freq{$pos}{$aa} = $line[4]; # Total Mutation Frequency for Naive Library
		
		#print "$pos $aa $a_freq{$pos}{$aa} $b_freq{$pos}{$aa}\n";
		
	}
	close(FILE);
	
	# All 20 standard amino acids
	my @aa_all = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
	
	# p := probability
	# f := frequency of amino acid
	# a := Amino acid type
	# N := total frequency of all amino acids
	# K := Number of amino acid types (the standard 20)
	
	# Given p_a = f_a / N
	# Shannon Entropy = -1 * K_SUM_a (p_a log_2 p_a)
	
	my $freq_type;
	my $N;
	my $pa;
	my $sum_log_pa;
	
	my $entropy_score;
	my %entropy;
	
	my $missing_data;
	
	for $pos (sort {$a <=> $b} keys %a_freq) {
		
		for $aa (sort {$a cmp $b} @aa_all) {
			$N += $a_freq{$pos}{$aa};
			# Keep track of missing mutations
			$missing_data++ if $a_freq{$pos}{$aa} < $mutation_th;
		}
		
		# Skip position if there is missing data (require all 20 minus WT mutations to be present)
		if ($missing_data > 1) {
			$missing_data = 0;
			next;
		}
		
		for $aa (sort {$a cmp $b} @aa_all) {
			$freq_type = $a_freq{$pos}{$aa};
			$pa = $freq_type / $N;
			$sum_log_pa += $pa * ( log($pa) / log(2) );
		}
		
		$entropy_score = -1 * $sum_log_pa;
		#print "$pos $entropy_score\n";
		
		$entropy{$pos} = $entropy_score;
		# Normalized
		#$entropy{$pos} = $entropy_score / ( log(20) / log(2) );
		
		$N = 0;
		$sum_log_pa = 0;
	}
	
	return %entropy;
}

sub ssm_ratio {
	# Calculate the sequence entropy of every residue position in the SSM "tab" input file
	
	my $ssm_file = shift;
	my $ratio_th = shift || 0;
	
	open(FILE,"$ssm_file") or die "Cannot open file $ssm_file\n";

	my %ratio_aa;
	
	my $pos;
	my $aa;
	my @line;
	
	while (<FILE>) {
		
		next if $. == 1;
		chomp;
		
		@line = split "\t", $_;
		
		# Skip WT amino acids
		#next if $line[8] eq 'WT';
		
		# Amino acid sequence position
		$pos = $line[1];
		# Start with 1
		$pos++;
		
		$aa = $line[2];
		$ratio_aa{$pos}{$aa} = $line[11];
		
	}
	close(FILE);
	
	# All 20 standard amino acids
	my @aa_all = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
	
	my %ratio;
	my $total;
	my $aa_num;
	
	for $pos (sort {$a <=> $b} keys %ratio_aa) {
		
		$total = 0;
		#$aa_num = 0;
		
		for $aa (sort {$a cmp $b} @aa_all) {
			$total += $ratio_aa{$pos}{$aa};
			#$aa_num++ if $ratio_aa{$pos}{$aa} > $ratio_th;
		}
		
		#print STDERR "$pos\t$total\n";
		
		$ratio{$pos} = $total;
		
	}
	
	return %ratio;
}

sub map2bfactor {
	
	# Map all scores to a pdb structure
	
	my $pdb_file = shift;
	my $chain = shift;
	my $pdist_ref = shift;
	my $mult = shift || 1;
	my $min_res_num = shift;
	my $max_res_num = shift;
	
	my %chain_id;
	if ($chain eq 'NA') {
		$chain_id{'NA'}++;
	}
	else {
		$chain_id{$_}++ for split //, $chain;
	}
	
	my $mer_id;
	my $res_num;
	my $score;
	my $default_score = "      "; # Formatted Blank
	my $map_score;
	my ($min_score, $max_score);
	
	open(PDB,"$pdb_file");
	my $map_out;
	
	while (<PDB>) {
		# If ATOM
		if (/^ATOM/) {
			if ($chain eq 'NA') {
				$mer_id = 'NA';
			}
			else {
				$mer_id = substr($_,21,1);
				$mer_id =~ s/ //g;
			}
			
			# If specificed chain
			if ($chain_id{$mer_id}) {
			
				$res_num = substr($_,22,8);
				$res_num =~ s/ //g;
				
				# If not within residue range
				if ($res_num < $min_res_num or $res_num >= $max_res_num && $min_res_num && $max_res_num) {
					$score = $default_score;
					substr($_,60,6,$score);
					$map_out .= $_;
				}
				# If scoring residue
				elsif (defined ${$pdist_ref}{$res_num}) {
					$map_score = ${$pdist_ref}{$res_num};
					$map_score *= $mult;
					#$map_score = 2**$map_score; # Make linear (from log_2)
					$score = sprintf("%6.2f",$map_score);
					substr($_,60,6,$score);
					$map_out .= $_;
					$min_score = $score if $min_score > $score;
					$max_score = $score if $max_score < $score;
				}
				# If not scoring residue
				else {
					$score = $default_score;
					substr($_,60,6,$score);
					$map_out .= $_;
				}
			}
			# If not specified chain
			else {
				$score = $default_score;
				substr($_,60,6,$score);
				$map_out .= $_;
			}
		}
		# If not ATOM
		else {
			$map_out .= $_;
		}
	}
	
	close(PDB);
	
	return $map_out;
}
