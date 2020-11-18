#!/usr/bin/perl

# Written by David La
# Updated: Wed Jul 29 15:43:07 PDT 2015

# Description:
# This script renumbers the PDB structure by its original 
# full length sequence residue numbering order!

my $usage = "renumPDBbySeq.pl <pdb_file> <pdb_chain> <original_sequence>\n";

use strict;

my $pdb_file = $ARGV[0] or die $usage;
my $pdb_chain = $ARGV[1] or die $usage;
my $full_seq = $ARGV[2] or die $usage;

# ---- Read PDB File ----
my %pdb_atom = readPDBAtom($pdb_file,$pdb_chain);

# ---- Find the Target PDB Sequence ----
my $pdb_seq = atom2seq(\%pdb_atom);

# Hash for Referencing PDB Residue Number to
# Chronologically Renumbered PDB Residue Number
my %pdb_num_2_pdb_renum = pdbSeqNum(\%pdb_atom);

# ---- Generate Sequence Alignment ----
my ($seq_aln1,$seq_aln2) = align_pair($full_seq,$pdb_seq);
#print "\n\n";
#print "$seq_aln1\n";
#print "$seq_aln2\n";

my %msa;
$msa{'full_seq_aln'} = $seq_aln1;
$msa{'pdb_seq_aln'} = $seq_aln2;

# Hash for Referencing Renumbered PDB Residue Number to MSA Position Number
my %pdb_renum_2_msa_num = msaSeqNum(\%msa,$pdb_seq);

# Hash for Referencing Full Length Residue Number to MSA Position Number
my %full_num_2_msa_num = msaSeqNum(\%msa,$full_seq);

# Hash for Referencing MSA position to Full Length Sequence Number
my %msa_num_2_full_num = reverse %full_num_2_msa_num;

# Convert PDB Residue Number to Full Length Residue Number
my ($pdb_num,$pdb_renum);
my ($msa_num,$full_num);
my %pdb_num_2_full_num;

for $pdb_num (sort {$a <=> $b} keys %pdb_num_2_pdb_renum) {
	# PDB Num => PDB Renum
	$pdb_renum = $pdb_num_2_pdb_renum{$pdb_num};
	# PDB Renum => MSA Num
	$msa_num = $pdb_renum_2_msa_num{$pdb_renum};
	# MSA Num => Full Num
	$full_num = $msa_num_2_full_num{$msa_num};

	# Hash Referencing PDB Number to Full Length Sequence Number
	$pdb_num_2_full_num{$pdb_num} = $full_num;
}

# ---- Output PDB with Scores Mapped to B-Factors ----
my $map_out = map2resnum_all_out($pdb_file,$pdb_chain,\%pdb_num_2_full_num);
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

sub fasta2hash {
	
	my $msa_file = shift;
	
	my $count;
	my $id;
	my %fsa; 
	my $len;
	#my $buffer_len;
	
	open(FSA,$msa_file) or die "$msa_file: $!";
	while (<FSA>) {
		tr/\r\n//d;
		if (/^>(.*)/) { 
			$count++;
			$id = $1;
			$id =~ s/ /\_/;
			#die "Not all MSA lengths are equal!\n" if $buffer_len && $len != $buffer_len;
			#$buffer_len = $len;
			$len = 0;
		}
		else { 
			$fsa{"$count\_$id"} .= $_;
			$len += length($_);
		}
	}
	close(FSA);
	
	return ($count,$len,\%fsa);

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

	rmdir($temp_dir);

	return ($aln_seq1,$aln_seq2);                                                               

}

sub msaSeqNum {
	
	my $msa_ref = shift;
	my $seq = shift;

	my $seq_id;
	my $found_id;
	my $msa_seq;
	my $msa_seq_degap;
	my $msa_aa;
	my $msa_count;
	my $seq_count;
	my %msa_num;
	
	for $seq_id (sort {$a <=> $b} keys %{$msa_ref}) {
		$msa_seq = ${$msa_ref}{$seq_id};
		$msa_seq_degap = $msa_seq;
		$msa_seq_degap =~ s/-//g;
		if ($seq eq $msa_seq_degap) {
			for $msa_aa (split //, $msa_seq) {
				$msa_count++;
				unless ($msa_aa eq '-') {
					$seq_count++;
					$msa_num{$seq_count} = $msa_count;
				}
			}
			last;
		}
	}
	
	die "Cannot find PDB sequence ($seq) in MSA\n" if !$msa_count;

	return %msa_num;
	
}

sub map2resnum_all_out {
	
	# Map all res_num_new to pdb structure
	
	my $pdb_file = shift;
	my $chain = shift;
	my $pdist_ref = shift;
	my $adder = shift || 0;
	
	my %chain_id;
	if ($chain eq 'NA') {
		$chain_id{'NA'}++;
	}
	else {
		$chain_id{$_}++ for split //, $chain;	
	}
	
	my $mer_id;
	my $res_num;
	my $res_num_new;
	my $map_res_num_new;
	
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
				
				# If scoring patch
				if (defined ${$pdist_ref}{$res_num}) {
					$map_res_num_new = ${$pdist_ref}{$res_num};
					$map_res_num_new += $adder;
					$res_num_new = sprintf("%4i",$map_res_num_new);
					substr($_,22,4,$res_num_new);
					$map_out .= $_;
				}
				# If not high scoring patch
				else {
					# Default Output
					$map_out .= $_;
				}
			}
			# If not specified chain
			else {
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
