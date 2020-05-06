#!/usr/bin/perl

# Written by David La
# Updated: Fri Jul 13 13:03:32 PDT 2018

# Description:
# This script will transfer residue numbering from one pdb (pdb_file_resnum_want) to another (pdb_file_resnum_get).
# Corresponding residues will be determined by pairwise sequence alignment.
# This is useful in case you want to renumber residues in your models in the same way as the xtal structure!

# NOTE: Sequence Trunction: Do not output if residue does not exist in the reference structure
# We expect these possible trunctions to be located in the N' or C' terminal ends and/or other missing densities
# in the xtal structure.  See the subroutine 'outputPDBAtom_resnum_conv()' for more info.

my $usage = "Usage: transfer_resnum.pl <pdb_file_resnum_want> <pdb_file_resnum_get>\n";

use strict;

my $pdb_num_want = $ARGV[0] or die $usage; # pdb_file_num_want (e.g. Xtal)
my $pdb_num_get = $ARGV[1] or die $usage; # pdb_file_num_get (e.g. Model)

#%pdb_atom = readPDBAtom($native_pdb,$nat_chain);
#$native_seq = atom2seq(\%pdb_atom);

#my ($native_aln2,$native_full_aln) = align_pair($native_seq,$native_full_seq);

my %chains_num_want = getChains($pdb_num_want);
my %chains_num_get = getChains($pdb_num_get);

my %pdb_atom_want;
my $pdb_seq_want;

my %pdb_atom_get;
my $pdb_seq_get;

my $out;
my $renum;

my ($pdb_aln_want, $pdb_aln_get);

my %msa;

my %seq2msa_want;
my %seq2msa_get;

my %msa2seq_want;
my %msa2seq_get;

my %pdb2seq_want;
my %pdb2seq_get;

my %seq2pdb_want;
my %seq2pdb_get;

my %pdb2msa_want;
my %pdb2msa_get;

my %msa2pdb_want;
my %msa2pdb_get;

my %get2want;
my %pdb_atom_get2want;

for my $mer (sort {$a cmp $b} keys %chains_num_get) {

	# Read PDB File (this also renumbers)
	%pdb_atom_want = readPDBAtom($pdb_num_want,$mer);
	$pdb_seq_want = atom2seq(\%pdb_atom_want);

	%pdb_atom_get = readPDBAtom($pdb_num_get,$mer);
	$pdb_seq_get = atom2seq(\%pdb_atom_get);

	# Pairwise Sequence Alignment
	#print "WANT SEQ: $pdb_seq_want\n";
	#print "GET SEQ: $pdb_seq_get\n";
	($pdb_aln_want, $pdb_aln_get) = align_pair($pdb_seq_want, $pdb_seq_get);
	#print "WANT ALN: $pdb_aln_want\n";
	#print "GET ALN: $pdb_aln_get\n";

	# MSA position number to real residue number
	%msa = ();
	$msa{'seq1'} = $pdb_aln_want; 
	$msa{'seq2'} = $pdb_aln_get; 

	%seq2msa_want = msaSeqNum(\%msa, $pdb_seq_want);
	%msa2seq_want = reverse %seq2msa_want;

	%seq2msa_get = msaSeqNum(\%msa, $pdb_seq_get);
	%msa2seq_get = reverse %seq2msa_get;

	# PDB residue number to real residue number
	%pdb2seq_want = pdbSeqNum(\%pdb_atom_want); # Xtal
	%seq2pdb_want = reverse %pdb2seq_want;

	%pdb2seq_get = pdbSeqNum(\%pdb_atom_get); # Model
	%seq2pdb_get = reverse %pdb2seq_get;

	# Get from residue num of PDB_want to residue_num of PDB_get
	%pdb2msa_want = transfer_pos(\%pdb2seq_want,\%seq2msa_want);
	#print "$mer: PDB2MSA WANT ------------------------------------------------------\n";
	%pdb2msa_get = transfer_pos(\%pdb2seq_get,\%seq2msa_get);
	#print "$mer: PDB2MSA GET ------------------------------------------------------\n";

	%msa2pdb_want = reverse %pdb2msa_want;
	%msa2pdb_get = reverse %pdb2msa_get;

	#%get2want = transfer_pos(\%pdb2msa_want,\%msa2pdb_get);
	%get2want = transfer_pos(\%pdb2msa_get,\%msa2pdb_want);
	#print "$mer: GET2WANT ------------------------------------------------------\n";

	# Output the reformatted PDB Atoms
	#$out = outputPDBAtom(\%pdb_atom_want);
	#$out = outputPDBAtom(\%pdb_atom_get);

	$out = outputPDBAtom_resnum_conv(\%pdb_atom_get,\%get2want);
	print $out;
}


# ---- Subs ----

sub atom2seq {
	
	my $pdb_atom_ref = shift;

	my $index;
	my $seq;
	my $res;
	my $aa;
	my $atom;
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
		# Skip if not CA atom
		$atom = ${$pdb_atom_ref}{$index}{'atom'};
		next if $atom ne 'CA';
		
		# Residue number
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

sub align_pair {
	# Use MUSCLE for Pairwise Sequence Alignment
	
	my $seq1 = shift;
	my $seq2 = shift;
	my $matrix = shift;
		
	my $muscle_exe = "/home/dla/bin/muscle";
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

sub pdbSeqNum {
	
	my $pdb_atom_ref = shift;

	my $index;
	my $res_num;
	my $res_num_rem;
	my $count;
	my %aa_num;
	
	for $index (sort {$a <=> $b} keys %{$pdb_atom_ref}) {
		$res_num = ${$pdb_atom_ref}{$index}{'res_num'};
		if ($res_num_rem ne $res_num) {
			$count++;
			$aa_num{$res_num} = $count;
			#print "$res_num $count\n";
		}
		$res_num_rem = $res_num;
	}

	return %aa_num;
	
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
					#print "$seq_count $msa_count\n";
				}
			}
			last;
		}
	}
	
	die "Cannot find PDB sequence ($seq) in MSA\n" if !$msa_count;

	return %msa_num;
	
}

sub transfer_pos {
	my $posA_ref = shift;
	my $posB_ref = shift;
	
	my %posA = %{$posA_ref};
	my %posB = %{$posB_ref};
	
	my %convert_pos;
	my ($posA_num,$posB_num);
	for $posA_num (sort {$a <=> $b} keys %posA) {
		#$convert_pos{ $posB{ $posA{$posA_num} } } = $posA_num;
		$convert_pos{ $posA_num } = $posB{ $posA{$posA_num}};
		#print "$posA_num => $posB{ $posA{$posA_num} }\n";
	}
	
	return %convert_pos;
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

sub getChains {
	my $pdb_file = shift;

	open(PDB,$pdb_file) or die "Cannot open PDB file: $pdb_file\n";

	my $mer;
	my %chains;

	while (<PDB>) {
		next unless /^ATOM/;

		$mer = substr($_,20,2);
		$mer =~ s/ //g;

		$chains{$mer}++;
	}

	close(PDB);

	return %chains;
}

sub readPDBAtom {
	
	my $pdb_file = shift;
	my $pdb_chain = shift;
	
	# -- Read PDB File

	my ($index,$atom,$aa,$mer,$res_num);
	my ($x,$y,$z);
	my ($occ,$bfact);
	my %pdb_atom;

	my $mer_rem;
	my $res_num_rem;
#my $renum;

	open(PDB,$pdb_file) or die "Cannot open PDB file: $pdb_file\n";
	while (<PDB>) {
		
		# Change MSE to MET, when appropriate
		if ($_ =~ /^HETATM/ && substr($_,17,3) eq "MSE") {
			# Hetero-Atom to Atom
			s/^HETATM/ATOM  /;
			# Selenium to Sulfur
			substr($_,12,2,"SD") if substr($_,12,2) eq "SE";
			# Senenomethionine to Methionine Residue
			substr($_,17,3,"MET");
		}
		
		# Change TYI to TYR, when appropriate
		if ($_ =~ /^HETATM/ && substr($_,17,3) eq "TYI") {
			# Hetero-Atom to Atom
			s/^HETATM/ATOM  /;
			# Skip Iodines
			next if substr($_,12,2) eq "I1";
			next if substr($_,12,2) eq "I2";
			# TYI to TYR Residue
			substr($_,17,3,"TYR");
		}

		# Change SEC to CYS, when appropriate
		if ($_ =~ /^HETATM/ && substr($_,17,3) eq "SEC") {
			# Hetero-Atom to Atom
			s/^HETATM/ATOM  /;
			# Selenium to Sulfur
			substr($_,12,2,"SG") if substr($_,12,2) eq "SE";
			# SEC to CYS Residue
			substr($_,17,3,"CYS");
		}
		
		# Skip non-atom lines
		next unless /^ATOM/;
		
		$index = substr($_,6,5);
		$index =~ s/ //g;
		
		# Parse necessary PDB data
		$atom = substr($_,11,5);
		$atom =~ s/ //g;
		
		$aa = substr($_,17,3);
		$aa =~ s/ //g;
		
		$mer = substr($_,20,2);
		$mer =~ s/ //g;
		
		# Skip chain if not needed unless it's not available (NA)
		unless ($pdb_chain eq 'NA') {
			next if $mer ne $pdb_chain;
		}
		
		$res_num = substr($_,22,5);
		$res_num =~ s/ //g;

		#---- Renumber ----
		#$renum = 0 if $mer ne $mer_rem;
		#$renum++ if $res_num ne $res_num_rem;

		#$mer_rem = $mer;
		#$res_num_rem = $res_num;

		#$res_num = $renum;
		#------------------
		
		$x = substr($_,30,8);
		$x =~ s/ //g;
		
		$y = substr($_,38,8);
		$y =~ s/ //g;
		
		$z = substr($_,46,8);
		$z =~ s/ //g;
		
		$occ = substr($_,54,6);
		$occ =~ s/ //g;
		
		$bfact = substr($_,60,6);
		$bfact =~ s/ //g;
		
		# Store PDB Coordinates
		$pdb_atom{$index}{'atom'} = $atom;
		$pdb_atom{$index}{'aa'} = $aa;
		$pdb_atom{$index}{'mer'} = $mer;
		$pdb_atom{$index}{'res_num'} = $res_num;
		$pdb_atom{$index}{'x'} = $x;
		$pdb_atom{$index}{'y'} = $y;
		$pdb_atom{$index}{'z'} = $z;
		$pdb_atom{$index}{'occ'} = $occ;
		$pdb_atom{$index}{'bfact'} = $bfact;
	
	}
	close(PDB);
	
	return %pdb_atom;

}
sub outputPDBAtom_resnum_conv {
	my $pdb_atom_ref = shift;
	my $res_num_conv_ref = shift;
	
	my %pdb_atom = %{$pdb_atom_ref};
	my %res_num_conv = %{$res_num_conv_ref};
	
	my ($index,$atom,$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
	
	my $out;
	
	for my $index (sort {$a <=> $b} keys %pdb_atom) {
		
		$atom = $pdb_atom{$index}{'atom'};
		$aa = $pdb_atom{$index}{'aa'};
		$mer = $pdb_atom{$index}{'mer'};

		# Convert Residue Number
		$res_num = $res_num_conv{ $pdb_atom{$index}{'res_num'} };

		$x = $pdb_atom{$index}{'x'};
		$y = $pdb_atom{$index}{'y'};
		$z = $pdb_atom{$index}{'z'};
		$occ = $pdb_atom{$index}{'occ'};
		$bfact = $pdb_atom{$index}{'bfact'};
		
		# NOTE: Sequence Trunction: Do not output if residue does not exist in the reference structure
		# We expect these possible trunctions to be located in the N' or C' terminal ends and/or other missing densities
		# in the xtal structure.
		next if $res_num == 0;

		$out .= sprintf("ATOM %6i ",$index);
		
		# Special conditional formatting for atom names
		if ($atom =~ /^\d+\S+\d+$/) {
			$out .= sprintf("%-3s",$atom);
		}
		elsif ($atom =~ /^\d+/) {
			$out .= sprintf("%-3s ",$atom);
		}
		elsif ($atom =~ /\d\d$/) {
			$out .= sprintf("%-3s",$atom);
		}
		elsif ($atom =~ /\w\w\d$/) {
			$out .= sprintf(" %-3s",$atom);
		}
		elsif ($atom =~ /\d$/) {
			$out .= sprintf("%-3s ",$atom);
		}
		else {
			$out .= sprintf(" %-3s",$atom);
		}

		if ($res_num =~ /[A-Za-z]$/) {
			$out .= sprintf("%4s %1s %4s%11.3f%8.3f%8.3f%6.2f%6.2f\n",$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
		}
		else {
			$out .= sprintf("%4s %1s%4i %11.3f%8.3f%8.3f%6.2f%6.2f\n",$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
		}
		
	}
	
	$out .= "TER\n";
	
	return $out;

}

sub outputPDBAtom {
	my $pdb_atom_ref = shift;
	
	my %pdb_atom = %{$pdb_atom_ref};
	
	my ($index,$atom,$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
	
	my $out;
	
	for my $index (sort {$a <=> $b} keys %pdb_atom) {
		
		$atom = $pdb_atom{$index}{'atom'};
		$aa = $pdb_atom{$index}{'aa'};
		$mer = $pdb_atom{$index}{'mer'};
		$res_num = $pdb_atom{$index}{'res_num'};
		$x = $pdb_atom{$index}{'x'};
		$y = $pdb_atom{$index}{'y'};
		$z = $pdb_atom{$index}{'z'};
		$occ = $pdb_atom{$index}{'occ'};
		$bfact = $pdb_atom{$index}{'bfact'};
		
		$out .= sprintf("ATOM %6i ",$index);
		
		# Special conditional formatting for atom names
		if ($atom =~ /^\d+\S+\d+$/) {
			$out .= sprintf("%-3s",$atom);
		}
		elsif ($atom =~ /^\d+/) {
			$out .= sprintf("%-3s ",$atom);
		}
		elsif ($atom =~ /\d\d$/) {
			$out .= sprintf("%-3s",$atom);
		}
		elsif ($atom =~ /\w\w\d$/) {
			$out .= sprintf(" %-3s",$atom);
		}
		elsif ($atom =~ /\d$/) {
			$out .= sprintf("%-3s ",$atom);
		}
		else {
			$out .= sprintf(" %-3s",$atom);
		}
		
		$out .= sprintf("%4s %1s%4i %11.3f%8.3f%8.3f%6.2f%6.2f\n",$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
		
	}
	
	$out .= "TER\n";
	
	return $out;

}

sub getChainRange {
	my $pdb_file = shift;

	open(FILE,"$pdb_file") or die "Cannot open file \n";

	my %chain_range;
	my $mer;
	my $res_num;
	my $mer_rem;
	my $res_num_rem;

	while (<FILE>) {

		chomp;

		next unless /^ATOM/;

		$res_num = substr($_,22,5);
		$res_num =~ s/ //g;

		$mer = substr($_,20,2);
		$mer =~ s/ //g;

		if ($mer_rem ne $mer) {
			$chain_range{$mer}{'begin'} = $res_num;
			$chain_range{$mer_rem}{'end'} = $res_num_rem if $res_num_rem;
		}

		$mer_rem = $mer;
		$res_num_rem = $res_num;
	}
	close(FILE);

	# Get the last end residue number for the final chain
	$chain_range{$mer_rem}{'end'} = $res_num_rem;

	# Output all of the change ranges
	#for my $chain (sort {$a cmp $b} keys %chain_range) {
		#print "$chain $chain_range{$chain}{'begin'}-$chain_range{$chain}{'end'}\n"
	#}

	return %chain_range;

}

