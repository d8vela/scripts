#!/usr/bin/perl

# Written by David La
# Updated Wed May 21 01:49:29 PDT 2014

# Description:
# This script will take an input SSM data (*.tab or PSSM-style) and generate
# an resfile that represents a desired combinatorial library

my $usage = "Usage: ssm2res.pl -ssm <ssm_file.tab> -pdb <pdb_file> -chain <pdb_chain>
		-pssm [PSSM-style file, an alternative to the -ssm flag]
		-ratio_th [enrichment ratio threshold (default: 1)]
		-interface [protein-protein interface distance threshold in Angstroms (default: < 8Å)]
		-surface [percent relative area of solvent accessibility threshold (default: > .10, meaning 10%)]
		-skip_pdb_pos [PDB residue number and mutation (e.g. 421,422,425) ]
		-skip_pdb_pos_mut [PDB residue number and mutation (e.g. 421Y,421M)]
		-skip_ssm_pos [SSM residue number (e.g. 21,22,25)]
		-skip_ssm_pos_mut [SSM residue number and mutation (e.g. 92R,92H)]
		-add_pdb_pos_mut [SSM residue number and mutation (e.g. 83R,108H)]
		-add_ssm_pos_mut [SSM residue number and mutation (e.g. 34K,34E)]
		-expand_codon_pos [SSM residue number dot expand number of desired codons (e.g. 92.2,93.6,94.all)]
		-db [path to aa_calc.txt degenerate codon database (default: /work/davela/data/codons/aa-calc.txt)]\n";

use strict;
use Getopt::Long;

# ---- Special Arguments ----
my %args=();
GetOptions(	"ssm:s"=>\$args{ssm},
			"pssm:s"=>\$args{pssm},
			"pdb:s"=>\$args{pdb},
			"chain:s"=>\$args{chain},
			"ratio_th:f"=>\$args{ratio_th},
			"interface:f"=>\$args{interface},
			"surface:f"=>\$args{surface},
			"skip_pdb_pos:s"=>\$args{skip_pdb_pos},
			"skip_pdb_pos_mut:s"=>\$args{skip_pdb_pos_mut},
			"skip_ssm_pos:s"=>\$args{skip_ssm_pos},
			"skip_ssm_pos_mut:s"=>\$args{skip_ssm_pos_mut},
			"add_pdb_pos_mut:s"=>\$args{add_pdb_pos_mut},
			"add_ssm_pos_mut:s"=>\$args{add_ssm_pos_mut},
			"expand_codon_pos:s"=>\$args{expand_codon_pos},
			"db:s"=>\$args{db});

# ---- Define Arguments ----
my $ssm_file = $args{ssm}; # SSM File (The Tab File)
my $pssm_file = $args{pssm}; # PSSM File (The PSSM File)
die $usage unless ($ssm_file or $pssm_file);

my $pdb_file = $args{pdb} or die $usage; # PDB File
my $pdb_chain = $args{chain} || 'B'; # PDB Chain ID
my $ratio_th = $args{ratio_th} || 1; # Enrichment Ratio Threshold (Either Logged or Not Logged)
my $interface_dist = $args{interface}; # Only Consider Residues of at a Distance Defined Protein Interface
my $skip_rsasa = $args{surface}; # Skip Relative Solvent Accessible Surface Area Equal to or Below a Percent Threshold
my $skip_resnum_list = $args{skip_pdb_pos}; # Skip PDB Residue Position (Comma Delimited List)
my $skip_mut_list = $args{skip_pdb_pos_mut}; # Skip Mutation at PDB Position (Comma Delimited List) [e.g. 421Y,421M]
my $skip_ssm_pos_list = $args{skip_ssm_pos}; # Skip SSM Position (Comma Delimited List)
my $skip_ssm_mut_list = $args{skip_ssm_pos_mut}; # Skip Mutation at SSM Position (Comma Delimited List) [e.g. 92R,92H]
my $add_pdb_mut_list = $args{add_pdb_pos_mut}; # Skip Mutation at SSM Position (Comma Delimited List) [e.g. 92R,92H]
my $add_ssm_mut_list = $args{add_ssm_pos_mut}; # Skip Mutation at SSM Position (Comma Delimited List) [e.g. 92R,92H]
my $expand_codon_pos_list = $args{expand_codon_pos}; # Expand Codon Choices at SSM Position dot Desired Number of Codons (Comma Delimited List) [e.g. 92.2,93.3,94.all]
my $db = $args{db} || '/work/davela/data/codons/aa-calc.txt'; # Pre-computed Degenerate Codon Database

# ---- Read PDB File ----
my %pdb_atom = readPDBAtom($pdb_file,$pdb_chain);

# ---- Get Protein Interface ----
my %bind = get_interface($pdb_file,$interface_dist) if $interface_dist;

# ---- Get Protein Surface ----
my %surf = dssp_surf($pdb_file,$pdb_chain,$skip_rsasa) if $skip_rsasa;

# ---- Get PDB Residue Postions to Skip ----
my %skip_resnum;
$skip_resnum{$_}++ for (split /\,/, $skip_resnum_list);

# ---- Get Mutations at PDB Positions to Skip ----
my %skip_mut;
my ($resnum,$mut_res);
for (split /\,/, $skip_mut_list) {
	# PDB Position and Mutation
	($resnum,$mut_res) = $_ =~ /(\d+)(\w+)/;
	$mut_res =~ tr/a-z/A-Z/;
	$skip_mut{$resnum}{$mut_res}++;
}

# ---- Get SSM Positions to Skip ----
my %skip_ssm_pos;
$skip_ssm_pos{$_}++ for (split /\,/, $skip_ssm_pos_list);

# ---- Get Mutations at PDB Positions to Skip ----
my %skip_ssm_mut;
my ($ssm_pos,$ssm_res);
for (split /\,/, $skip_ssm_mut_list) {
	# PDB Position and Mutation
	($ssm_pos,$ssm_res) = $_ =~ /(\d+)(\w+)/;
	$ssm_res =~ tr/a-z/A-Z/;
	$skip_ssm_mut{$ssm_pos}{$ssm_res}++;
}

# ---- Get Positions to Expand Codons Choices ----
my ($expand_pos,$expand_num);
my %expand_codon_pos;
for my $expand_data (split /\,/, $expand_codon_pos_list) {
	($expand_pos,$expand_num) = split /\./, $expand_data;
	$expand_codon_pos{$expand_pos} = $expand_num;
}

# ---- Get the Design PDB Sequence ----
my $design_pdb_seq = atom2seq(\%pdb_atom);

# ---- Check SSM or PSSM Data ----
my $design_ssm_seq;
my $design_pssm_seq;

my $pdb_len;
my $pssm_len;

my $design_pdb_seq_aln;
my $design_ssm_seq_aln;
my $design_pssm_seq_aln;
my %design_aln;

if ($ssm_file) {
	# ---- Get the Design SSM Sequence (possibly with residues of missing densities) ----
	$design_ssm_seq = ssm2seq($ssm_file);

	# ---- Align PDB to SSM WT Sequence ----
	($design_pdb_seq_aln,$design_ssm_seq_aln) = align_pair($design_pdb_seq,$design_ssm_seq);

	%design_aln = ('PDB' => $design_pdb_seq_aln,
		       'SSM' => $design_ssm_seq_aln);
}
elsif ($pssm_file) {
	
	# ---- Check that PSSM Length is the same as the PDB Length ----
	
	# PDB Length
	my $pdb_len = length($design_pdb_seq);
	
	# PSSM Length
	my $pssm_len = column_count($pssm_file);
	
	# Check PDB Length is the same as PSSM Length
	die "\nYour PDB Length ($pdb_len) is not the Same as the PSSM Length ($pssm_len)!\n" unless $pdb_len == $pssm_len;
	
	# ---- Assume PDB and SSM sequence are the same for PSSM Inputs ----
	$design_ssm_seq = $design_pdb_seq; # Assume the same as PDB
	
	$design_pdb_seq_aln = $design_pdb_seq; # Assume the same as PDB
	$design_ssm_seq_aln = $design_pdb_seq; # Assume the same as PDB
	
	%design_aln = ('PDB' => $design_pdb_seq_aln,
		       'SSM' => $design_pssm_seq_aln);
}
else {
	die "Error: No input SSM or PSSM file specified!\n";
}

# ---- Reference PDB and Alignment Positions to Real Residue Positions ----

# Hash for Referencing PDB Residue Number to Real Residue Number
my %aa_num = pdbSeqNum(\%pdb_atom);

# Hash for Referencing Real Residue Number to Alignment Position Number
my %aln_num = alnSeqNum(\%design_aln,$design_pdb_seq);

# Convert Alignment Positions to PDB Positions
my %aln2aa = reverse %aln_num;
my %aa2pdb = reverse %aa_num;

# ---- Mutations Added by Choice ----

my %add_mut;

# Add PDB Position Mutations
%add_mut = add_mut_list($add_pdb_mut_list);
# Convert PDB Position to SSM Positions
%add_mut = convert_pos_mut(\%add_mut,\%aa_num);

# Add SSM Position Mutations
%add_mut = add_mut_list($add_ssm_mut_list,\%add_mut);

# ---- Fetch Enriched SSM Mutations ----

my %enrich;
if ($ssm_file) {
	%enrich = ssm_enrich($ssm_file,$ratio_th,\%add_mut);
}
elsif ($pssm_file) {
	%enrich = pssm_enrich($pssm_file,$ratio_th,\%add_mut);
}
else {
	die "Error: No input SSM or PSSM file specified!\n";
}

# ---- Output Resfile of the Combinatorial Library ----
my $count;
my $aa_pos;
my $pdb_pos;
my $pos;
my $mut_aa;
my $new_mut_aa;
my $codon;
my $codon_count;
my $codon_comb = 1;
my $mut_len;
my $new_mut_len;
my $mut_comb = 1;
my $new_mut_comb = 1;
my $wt_aa;

print "NATRO\n";
print "start\n";
for $wt_aa (split //, $design_ssm_seq) {
	# Increment residue position
	$pos++;
	$aa_pos = $aln2aa{$pos}; # Get Real AA position
	$pdb_pos = $aa2pdb{$aa_pos}; # Get PDB AA position
	
	# Skip Non-interface Residues
	if ($interface_dist) {
		next unless $bind{$pdb_chain}{$pdb_pos};
	}
	
	# Skip Core Residues
	if ($skip_rsasa) {
		next unless $surf{$pdb_pos};
	}
	
	# Skip defined PDB residue positions
	next if $skip_resnum{$pdb_pos};
	
	# Get List of Enriched Mutations
	$mut_aa = $enrich{$pos};
	
	# Skip Defined Mutations at PDB Positions
	for my $aa (split //, $mut_aa) {
		# Remove residue from the mutations list
		if ($skip_mut{$pdb_pos}{$aa}) {
			$mut_aa =~ s/$aa//g; # "g" incase there are duplicates
		}
	}
	
	# Skip Defined SSM Positions
	next if $skip_ssm_pos{$pos};
	
	# Skip Defined Mutations at SSM Positions
	for my $aa (split //, $mut_aa) {
		# Remove residue from the mutations list
		if ($skip_ssm_mut{$pos}{$aa}) {
			$mut_aa =~ s/$aa//g; # "g" incase there are duplicates
		}
	}
	
	# Place WT in the Front of the Mutation List
	$mut_aa =~ s/$wt_aa//; # Temperarily, remove WT residue (if it's there!)
	$mut_aa =~ s/^/$wt_aa/; # Now, place WT residue in the front
	
	# Skip Position if it has No Mutations
	$mut_len = length($mut_aa);
	next if $mut_len == 1;
	
	# Get Degenerate Codon for the Mutation List
	if ($expand_codon_pos{$pos}) {
		$expand_num = $expand_codon_pos{$pos};
		($codon,$new_mut_aa,$codon_count) = expand_codon($mut_aa,$expand_num,$db);
	}
	else {
		($codon,$new_mut_aa,$codon_count) = aa_calc($mut_aa,$db);
	}
	
	# Calculate Number of Mutation Combinations
	$mut_comb *= $mut_len;
	
	$new_mut_len = length($new_mut_aa);
	$new_mut_comb *= $new_mut_len;
	
	$codon_comb *= $codon_count;
	
	# Place WT in the Front (Again)
	$mut_aa =~ s/$wt_aa//;
	$mut_aa =~ s/^/\[$wt_aa\]/;
	$new_mut_aa =~ s/$wt_aa//;
	$new_mut_aa =~ s/^/$wt_aa/;
	
	# Output Mutations
	$count++;
	printf("%s\t%s\t%s\t%-20s\t%s\n",$pdb_pos,$pdb_chain,"PIKAA",$new_mut_aa,"# CODON=$codon POS=$pos ENRICH=$mut_aa");
	
}

# Reformat Numbers with Commas
$mut_comb = num_commify($mut_comb);
$new_mut_comb = num_commify($new_mut_comb);
$codon_comb = num_commify($codon_comb);

# Output Stats
print "\n";
print "# Number of Residue Positions: $count\n";
print "# Enrichment Library Diversity: $mut_comb\n";
print "# Protein Library Diversity: $new_mut_comb\n";
print "# DNA Library Diversity: $codon_comb\n\n";

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

sub get_interface {
	
	my $pdb_file = shift;
	my $dist_th = shift || 8.0;
	
	my $contato = "/work/davela/src/contact1.2/contato";
	
	my ($chain,$aa_num);
	my %bind;

	for (`$contato -pdb $pdb_file -distance $dist_th 2> /dev/null`) {
		next unless /^\w\d+\s+\w+/;
		chomp;
		($chain,$aa_num) = $_ =~ /^(\w)(\d+)\s\w+/;

		$bind{$chain}{$aa_num}++;
	}
	
	return %bind;
}

sub dssp_surf {
	
	my $pdb_file = shift;
	my $pdb_chain = shift;
	my $prcnt_asa_th = shift || .10;
	my $temp = shift || "/tmp/DSSP_$$";
	
	#my $dssp = "$ENV{'HOME'}/bin/dsspcmbi";
	my $dssp = "/work/davela/bin/dssp";
	
	my %total_asa = (
		'A' => 113,
		'R' => 241,
		'N' => 158,
		'D' => 151,
		'C' => 140,
		'Q' => 189,
		'E' => 183,
		'G' => 85,
		'H' => 194,
		'I' => 182,
		'L' => 180,
		'K' => 211,
		'M' => 204,
		'F' => 218,
		'P' => 143,
		'S' => 122,
		'T' => 146,
		'W' => 259,
		'Y' => 229,
		'V' => 160
	);

	# Run DSSP
	`mkdir $temp` unless -d $temp;
	
	# Change PDB file to PDB file with only chain of interest
	$pdb_file = parseChain($pdb_file,$pdb_chain,$temp);
	
	# Find atoms on the protein surface
	my @line;
	my %surf;
	my ($res_num,$res,$asa,$prcnt_asa);
	
	my $surf_count;
	my $all_count;

	my $start;
	for (`$dssp $pdb_file 2> /dev/null`) {
		if ($start) {
			$res_num = substr($_,5,6);
			$res_num =~ s/ //g;
			
			$res = substr($_,13,1);
			
			$asa = substr($_,35,4);
			$asa =~ s/ //g;
			
			next unless $total_asa{$res};
			$prcnt_asa = $asa / $total_asa{$res};
			
			if ($prcnt_asa >= $prcnt_asa_th) {
				$surf{$res_num}++;
				$surf_count++;
			}
			
			$all_count++;
		}
		elsif (/\s+#/) {
			$start++;
		}
	}
	
	# Clean TMP
	`rm -r $temp/*`;
	
	return(%surf);
	
}

sub parseChain {
	
	my $pdb_file = shift;
	my $chain_id = shift;
	my $out_dir = shift;
	
	my $mer;
	
	open(PDB,$pdb_file);
	open(CHAIN,">$out_dir/$chain_id.pdb");
	while (<PDB>) {
		next unless /^ATOM/;
		$mer = substr($_,21,1);
		unless ($chain_id eq 'NA') {
			if ($mer ne $chain_id) {
				next;
			}
		}
		print CHAIN;
	}
	close(CHAIN);
	close(PDB);
	
	return "$out_dir/$chain_id.pdb";

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

sub convert_pos_mut {
	# Convert Position Mutation Positions Numbering
	my $pos_mut_ref = shift;
	my $pos_convert_ref = shift;

	my %pos_mut = %{$pos_mut_ref};
	my %pos_convert = %{$pos_convert_ref};
	
	my %pos_mut_convert;
	for my $pos (sort {$a <=> $b} keys %pos_mut) {
		for my $aa (sort {$a <=> $b} keys %{$pos_mut{$pos}}) {
			$pos_mut_convert{$pos_convert{$pos}}{$aa} = $pos_mut{$pos}{$aa}
		}
	}

	return %pos_mut_convert;
}

sub ssm_enrich {
	
	my $ssm_file = shift;
	my $ratio_th = shift || 1;
	my $add_mut_ref = shift;

	my %add_mut = %{$add_mut_ref};
	
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

	# For First Position to the Last Position of the Protein Sequence
	my %comb;
	for $pos (sort {$a <=> $b} keys %ratio_aa) {
		# For Each of the Standard 20 Amino Acids + 1 Stop Codon (Sorted Largest to Smallest Enrichment Value)
		for $aa (sort { $ratio_aa{$pos}{$b} <=> $ratio_aa{$pos}{$a} } keys %{$ratio_aa{$pos}}) {
			if ($add_mut{$pos}{$aa}) {
				$comb{$pos} .= $aa;
			}
			elsif ($ratio_aa{$pos}{$aa} <= $ratio_th) {
				next;
			}
			else {
				$comb{$pos} .= $aa;
			}
		}
		
	}
	
	return %comb;
}

sub pssm_enrich {
	
	my $pssm_file = shift;
	my $ratio_th = shift || 1;
	my $add_mut_ref = shift;

	my %add_mut = %{$add_mut_ref};
	
	# All 20 standard amino acids + '*' stop codon
	my @aa_all = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*');
	
	open(FILE,"$pssm_file") or die "Cannot open file $pssm_file\n";

	my %ratio_aa;
	
	my $aa_type = 0;
	my $pos;
	my $aa;
	my $enrich;

	while (<FILE>) {
		chomp;
		
		for $enrich (split /\s+/, $_) {
			$pos++;
			$aa = $aa_all[$aa_type];
			$ratio_aa{$pos}{$aa} = $enrich;
		}
		
		# Reset Position
		$pos = 0;
		
		# Increment to next AA type
		$aa_type++;
		
	}
	close(FILE);
	
	# For First Position to the Last Position of the Protein Sequence
	my %comb;
	for $pos (sort {$a <=> $b} keys %ratio_aa) {
		# For Each of the Standard 20 Amino Acids + 1 Stop Codon (Sorted Largest to Smallest Enrichment Value)
		for $aa (sort { $ratio_aa{$pos}{$b} <=> $ratio_aa{$a} } keys %{$ratio_aa{$pos}}) {
			if ($add_mut{$pos}{$aa}) {
				$comb{$pos} .= $aa;
			}
			elsif ($ratio_aa{$pos}{$aa} <= $ratio_th) {
				next;
			}
			else {
				$comb{$pos} .= $aa;
			}

		}
		
	}
	
	return %comb;
}

sub add_mut_list {
	my $add_mut_list = shift;
	my $append_mut_ref = shift;

	my %append_mut;
	%append_mut = %{$append_mut_ref} if $append_mut_ref;

	# Process Mutation List to Add
	my %add_mut;
	%add_mut = %append_mut if $append_mut_ref;
	my ($pos,$aa);
	for my $mutation (split /\,/, $add_mut_list) {
		($pos,$aa) = $mutation =~ /(\d+)(\w+)/;
		next unless $mutation =~ /[ACDEFGHIKLMNPQRSTVWY\*]/;
		$aa =~ tr/a-z/A-Z/;
		$add_mut{$pos}{$aa}++;
	}

	return %add_mut;
}

sub aa_calc {
	my $mut_list = shift; # List of amino acids with no space (e.g: YMWH)
	my $db = shift || '/work/davela/data/codons/aa-calc.txt'; # Pre-computed database of degenerate codons!
	
	# Check if the pre-computed degenerate codon database exists
	die "Cannot find the precomputed database of degenerate codons!\n" unless -f $db;
	
	# For converting single amino acid letters to three amino acid letter codes!
	my %aa3code = (
		'A' => 'Ala',
		'C' => 'Cys',
		'D' => 'Asp',
		'E' => 'Glu',
		'F' => 'Phe',
		'G' => 'Gly',
		'H' => 'His',
		'I' => 'Ile',
		'K' => 'Lys',
		'L' => 'Leu',
		'M' => 'Met',
		'N' => 'Asn',
		'P' => 'Pro',
		'Q' => 'Gln',
		'R' => 'Arg',
		'S' => 'Ser',
		'T' => 'Thr',
		'V' => 'Val',
		'W' => 'Trp',
		'Y' => 'Tyr',
		'*' => 'stop'
	);
	
	# For converting three amino acid letter codes to single amino acid letters
	my %aa_code = reverse %aa3code;
	
	# Find the Best Degenerate Codon (with Minimal Number of Codons and Minimal Additional Amino Acids)
	my @query;
	
	for my $aa (split //, $mut_list) {
		push @query, $aa3code{$aa};
	}
	
	# Simply Grep-out your codon of interest from the database and take the first hit.
	my $mut_num;
	my $cmd;
	for my $res (@query) {
		$mut_num++;
		if ($mut_num == 1) {
			$cmd = "grep $res $db";
		}
		else {
			$cmd .= " | grep $res";
		}
	}
	$cmd .= " | sort -nk 2 | head -n 1";
	
	# Run series of grep commands
	my $hit = `$cmd`;
	chomp $hit;
	
	my ($codon,$codon_total,$aa_total,$aa_list);
	($codon,$codon_total,$aa_total,$aa_list) = $hit =~ /(\w+)\s+(\d+)\s+(\d+)\s+(.*)/;
	
	# Convert three letter amino acid codes to single letter amino acid codes.
	my ($three_letter_aa,$codon_count);
	my $new_mut_list;
	for my $aa_data (split /\s/, $aa_list) {
		($three_letter_aa,$codon_count) = split /\:/, $aa_data;
		$new_mut_list .= $aa_code{$three_letter_aa};
	}
	
	return ($codon,$new_mut_list,$codon_total);
	
}

sub expand_codon {
	my $mut_list = shift;
	my $expand_num = shift || 2;
	my $db = shift || '/work/davela/data/codons/aa-calc.txt'; # Pre-computed database of degenerate codons!
	
	# Residue Classes
	my @hydrophobics = ('A','V','L','I','M','C');
	my @aromatics = ('F','Y','W','H');
	my @polars = ('S','T','N','Q');
	my @acids = ('D','E');
	my @bases = ('K','R');
	my @special = ('G','P');
	
	my %class;
	my $count_all;
	
	my $aa;
	if ($expand_num == 2) {
		
		# Approximate Classes
		
		# Class 1: Hydrophobics
		for $aa (@hydrophobics) {
			$class{$aa} = 1;
		}
		for $aa (@aromatics) {
			$class{$aa} = 1;
		}
		for $aa (@special) {
			$class{$aa} = 1;
		}
		
		# Class 2: Hydrophilics
		for $aa (@polars) {
			$class{$aa} = 2;
		}
		for $aa (@acids) {
			$class{$aa} = 2;
		}
		for $aa (@bases) {
			$class{$aa} = 2;
		}
	}
	
	elsif ($expand_num == 3) {
		
		# Approximate Classes
		
		# Class 1: Hydrophobics
		for $aa (@hydrophobics) {
			$class{$aa} = 1;
		}
		for $aa (@special) {
			$class{$aa} = 1;
		}
		
		# Class 2: Aromatics
		for $aa (@aromatics) {
			$class{$aa} = 2;
		}
		
		# Class 3: Hydrophilics
		for $aa (@polars) {
			$class{$aa} = 3;
		}
		for $aa (@acids) {
			$class{$aa} = 3;
		}
		for $aa (@bases) {
			$class{$aa} = 3;
		}
	}
	
	elsif ($expand_num == 4) {
		
		# Approximate Classes
		
		# Class 1: Hydrophobics
		for $aa (@hydrophobics) {
			$class{$aa} = 1;
		}
		for $aa (@special) {
			$class{$aa} = 1;
		}
		
		# Class 2: Aromatics
		for $aa (@aromatics) {
			$class{$aa} = 2;
		}
		
		# Class 3: Polars
		for $aa (@polars) {
			$class{$aa} = 3;
		}
		
		# Class 4: Charged
		for $aa (@acids) {
			$class{$aa} = 4;
		}
		for $aa (@bases) {
			$class{$aa} = 4;
		}
	}
	
	elsif ($expand_num == 5) {
		
		# Approximate Classes
		
		# Class 1: Hydrophobics
		for $aa (@hydrophobics) {
			$class{$aa} = 1;
		}
		for $aa (@special) {
			$class{$aa} = 1;
		}
		
		# Class 2: Aromatics
		for $aa (@aromatics) {
			$class{$aa} = 2;
		}
		
		# Class 3: Polars
		for $aa (@polars) {
			$class{$aa} = 3;
		}
		
		# Class 4: Acids
		for $aa (@acids) {
			$class{$aa} = 4;
		}
		
		# Class 5: Bases
		for $aa (@bases) {
			$class{$aa} = 5;
		}
		
	}
	
	elsif ($expand_num == 6) {
		
		# Approximate Classes
		
		# Class 1: Hydrophobics
		for $aa (@hydrophobics) {
			$class{$aa} = 1;
		}
		
		# Class 2: Aromatics
		for $aa (@aromatics) {
			$class{$aa} = 2;
		}
		
		# Class 3: Polars
		for $aa (@polars) {
			$class{$aa} = 3;
		}
		
		# Class 4: Acids
		for $aa (@acids) {
			$class{$aa} = 4;
		}
		
		# Class 5: Bases
		for $aa (@bases) {
			$class{$aa} = 5;
		}
		
		# Class 6: Special
		for $aa (@special) {
			$class{$aa} = 6;
		}
	}
	
	elsif ($expand_num eq 'all') {
		
		# Classify All Residues as their Own Group
		
		for $aa (@hydrophobics,@aromatics,@polars,@acids,@bases,@special) {
			$count_all++;
			$class{$aa} = $count_all;
		}
	}
	
	# Group Mutations List
	my $group_num;
	my %muts;
	for $aa (split //, $mut_list) {
		$group_num = $class{$aa};
		$muts{$group_num} .= $aa;
	}
	
	# Convert Mutations into Codons
	my $degen_codons;
	my $new_codon;
	my ($new_mut_aa,$codon_count);
	my ($all_mut_aa,$all_codon_count);
	my $group_mut_list;
	for $group_num (keys %muts) {
		$group_mut_list = $muts{$group_num};
		($new_codon,$new_mut_aa,$codon_count) = aa_calc($group_mut_list,$db);
		$all_mut_aa .= $new_mut_aa;
		$all_codon_count += $codon_count;
		$degen_codons .= "$new_codon,";
	}
	
	# Remove Last Comma
	chop $degen_codons;
	
	return ($degen_codons,$all_mut_aa,$all_codon_count);

}

sub num_commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

sub column_count {
	my $file = shift;
	
	open(FILE,"$file") or die "Cannot open file $file\n";
	
	my $col_count;
	while (<FILE>) {
		chomp;
		for (split /\s+/) {
			$col_count++;
		}
		last;
	}
	close(FILE);
	
	return $col_count;
}

sub output_hash {
	my $hash_ref = shift;
	
	for (sort {$a cmp $b} keys %{$hash_ref}) {
		print "$_\t${$hash_ref}{$_}\n";
	}
	
}

