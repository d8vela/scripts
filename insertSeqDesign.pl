#!/usr/bin/perl

# Written by David La
# Updated Wed Mar 5 16:58:20 PST 2014


# Description: 
# This script will take a PDB file of your design and structurally align it with the 
# corresponding full length structure (with missing densities), then inserts the residues 
# of missing densities (gaps in the alignment) into the design sequence and outputs it in 
# the FASTA format.

# Note for AUTO DETECTION of the native PDB ID and Chain ID: 
# The first 5 characters of the file name WILL MATTER! It must start with a 4 letter PDB ID, 
# followed by the chain ID (no space).
# For example: The PDB ID "7TIM" and its corresponding Chain ID "A" should use the case 
# insensitive ID that corresponds to the file name: "7TIMA.pdb","7TIMA.design.pdb", or "7tima.blah.pdb".

my $usage = <<END;

\033[1m USAGE: \033[0m
	\033[0m insertSeqDesign.pl \033[0m[OPTIONS] <PDB_FILE> ...

\033[1m OPTIONS: \033[0m
	 -chain \e[4;37mCHAIN ID\e[4;0m			Chain ID of the input PDB file (default chain: B)
	 -nat_id \e[4;37mPDB ID\e[4;0m 			PDB ID of the native structure (default auto-detect from input PDB file name)
	 -nat_chain \e[4;37mCHAIN ID\e[4;0m			Chain ID of the native structure (default auto-detect from input PDB file name)
	 -insert_N \e[4;37m[ TRUE | FALSE ]\e[4;0m		Insert residues of missing densities at the N-terminal (default: FALSE)
	 -insert_C \e[4;37m[ TRUE | FALSE ]\e[4;0m		Insert residues of missing densities at the C-terminal (default: FALSE)
	 -output \e[4;37mFILE NAME\e[4;0m			Custom name of the output FASTA file (default file output name: designs.fasta)

\033[1m EXAMPLES: \033[0m
	\033[1m insertSeqDesign.pl 7tima.design.pdb \033[0m
	 
	    The above will use default chain B of the input "7tima.design.pdb" file, auto-detect native PDB ID has "7TIM" and 
	    chain ID is "A" and output results to "designs.fasta" as set by default
	 
	\033[1m insertSeqDesign.pl design.pdb -chain A -nat_id 7TIM -nat_chain A -output 7tim.fasta \033[0m
	 
	    The above will use chain ID "A" of the input "design.pdb" file, native PDB ID "7TIM", and native chain "A" and output to 
	    custom name file "7tim.fasta"
	
	\033[1m insertSeqDesign.pl *.pdb \033[0m

	    The above will use chain "B" of the input "design.pdb" file, auto-detect the PDB ID, native chain ID, and will output 
	    all the results for all PDB files to "designs.fasta" as set by default
	 
	 
END


use strict;
use lib '/work/davela/perl5/lib/perl5';
use Getopt::Long;
use File::Basename;


# -- Special Arguments
my %args=();
GetOptions(	"chain:s"=>\$args{chain},
		"nat_id:s"=>\$args{nat_id},
		"nat_chain:s"=>\$args{nat_chain},
		"insert_N:s"=>\$args{insert_N},
		"insert_C:s"=>\$args{insert_C},
		"output:s"=>\$args{output});

# For protein-protein interface design, the structure of the design is usually chain B
my $chain = $args{chain} || 'B';


# Check chain information
if ($chain) {
	die "ERROR: The -nat_chain value is suppose to a single character\n" unless $chain =~ /^\S$/;
}

# Defined native information
my $nat_id_def = $args{nat_id};
my $nat_chain_def = $args{nat_chain};

# Check native information
if ($nat_id_def or $nat_chain_def) {
	die "ERROR: Cannot use -nat_id and -nat_chain with more than one PDB input!\n" if $#ARGV > 0;
	die "ERROR: The -nat_id value is not in the proper 4 character format\n" unless $nat_id_def =~ /^\d\w{3}$/;
	die "ERROR: The -nat_chain value is suppose to a single character\n" unless $nat_chain_def =~ /^\S$/;
}

# Define if N-terminal and C-terminal missing densities should be inserted (default: No!)
my $insert_N = $args{insert_N} || 'FALSE';
my $insert_C = $args{insert_C} || 'FALSE';

# Check if flags to the optional values are in the correct format
if ($insert_N !~ /^TRUE|FALSE$/) {
	die "ERROR: The -insert_N value should be either 'TRUE' or 'FALSE'\n";
}
if ($insert_C !~ /^TRUE|FALSE$/) {
	die "ERROR: The -insert_C value should be either 'TRUE' or 'FALSE'\n";
}

# Output file name
my $output_file = $args{output} || "designs.fasta";


# -- Other Arguments
my @pdb_files = @ARGV or die $usage;


# -- Default variables
my $matrix = "/work/davela/data/matrix/IDENTITY";
my $tmp_dir = "/tmp/insertSeqDesign\_$$\_$ENV{'HOSTNAME'}";
#my $tmp_dir = "insertSeqDesign\_$$\_$ENV{'HOSTNAME'}";


# ---- Main ----

my ($design_pdb,$native_pdb);
my ($nat_id,$nat_chain);
my ($design_tmp,$native_tmp);
my ($design_aln,$native_aln1);
my ($native_full_aln,$native_aln2);
my ($native_seq,$native_full_seq);
my ($native_len,$native_full_len);
my ($seq_test1,$seq_test2);
my ($missing_res_count1,$missing_res_count2);
my $seq_ident;

my (%pdb_atom);
my $pdb_out;
my ($seq1,$seq2);
my $input_pdb;
my $input_file;
my $design_seq;
my $design_file;
my $design_full_seq;
my $missing_count;
my $warning_ref;
my $pos_des;

# Make temp directory
mkdir("$tmp_dir");

# Create FASTA file
open(FASTA,">$output_file") or die "Cannot create file: $output_file\n";

# STDOUT
print "\n------------ DESCRIPTION ------------\n\n";
print "The insertSeqDesign.pl script will insert residues of missing densities into your design sequences and output them in the FASTA format.\n\n";

print "------------ OPTIONS ------------\n\n";
print "Input PDB chain ID: $chain\n";
if ($nat_id_def and $nat_chain_def) {
	print "Native PDB ID: $nat_id_def\n";
	print "Native PDB chain ID: $nat_chain_def\n";
}
else {
	print "Native PDB ID: AUTO DETECT\n";
	print "Native PDB chain ID: AUTO DETECT\n";
}
print "Output file name: $output_file\n";
print "\n";

print "------------ STARTING ------------\n";

# Go through all input PDB files!
for $input_pdb (@pdb_files) {
	
	# Check if PDB file exist
	die "\nERROR: The PDB file \"$input_pdb\" cannot be found!\n\n" unless -f $input_pdb;
	
	print "\nProcessing PDB: $input_pdb\n";
	print "\n";
	
	# Try to automatically determine PDB ID and PDB Chain from the file name
	unless ($nat_id_def and $nat_chain_def) {
		$input_file = basename($input_pdb);
		($nat_id,$nat_chain) = $input_file =~ /^(\w{4})(\w)/;
		
		# Check PDB ID and chain ID format
		die "ERROR: The PDB ID could not be automatically detected from your input PDB file name\n" unless $nat_id =~ /^\d\w{3}$/;
		die "ERROR: The chain ID could not be automatically detected from your input PDB file name\n" unless $nat_chain =~ /^\S$/;
	}
	else {
		$nat_id = $nat_id_def;
		$nat_chain = $nat_chain_def;
	}
	
	# Reformat PDB ID and PDB Chain
	$nat_id =~ tr/a-z/A-Z/;
	$nat_chain =~ tr/a-z/A-Z/ unless $nat_chain_def; # Override Auto Cap of Chain ID if it is Pre-Defined!
	
	#print "PDB: $nat_id $nat_chain\n";
	
	# Get Specific Chain of the Design PDB Structure (Input PDB)
	$design_tmp = "$tmp_dir/design";
	mkdir("$design_tmp");
	$design_pdb = parseChain($input_pdb,$chain,"$tmp_dir/design");
	
	# Get Specific Chain of the Native PDB Structure
	$native_tmp = "$tmp_dir/native";
	mkdir("$native_tmp");
	$native_pdb = fetchPDB($nat_id,$nat_chain,"$tmp_dir/native");
	
	#print "PDB1: $design_pdb\n";
	#print "PDB2: $native_pdb\n";
	
	# Structure alignment between design and native PDBs (output as sequence alignment)
	($design_aln,$native_aln1) = TMalign($design_pdb,$native_pdb);
	$seq_ident = id_seq($design_aln,$native_aln1);
	$seq_ident = sprintf("%4.3f",$seq_ident);
	
	# Get Specific Chain of the Native Sequence (without the possible missing densities)
	%pdb_atom = readPDBAtom($native_pdb,$nat_chain);
	$native_seq = atom2seq(\%pdb_atom);
	
	# Get Specific Chain of the FULL Native Sequence (with the possible missing densities)
	$native_full_seq = pdbSeq_RESTful($nat_id,$nat_chain);
	
	# Check the assumption that native_seq is shorter or equal to native_full
	$native_len = length($native_seq);
	$native_full_len = length($native_full_seq);
	die "ERROR: Native sequence length is shorter than the full native sequence length (with missing densities)!\n" unless $native_len <= $native_full_len;
	
	# Sequence alignment between native and native_full using the IDENTITY matrix
	($native_aln2,$native_full_aln) = align_pair($native_seq,$native_full_seq);
	
	# Check that the above aligned native sequence vs design and full native sequence are identical
	$seq_test1 = $native_aln1;
	$seq_test1 =~ s/-//g;
	$seq_test2 = $native_aln2;
	$seq_test2 =~ s/-//g;
	
	#print "TEST1: $seq_test1\n";
	#print "TEST2: $seq_test2\n";
	die "ERROR: Native sequences previously generated for structure and sequence alignments are not identical!\n" unless $seq_test1 eq $seq_test2;
	
	# FASTA Output
	$missing_res_count1 = $native_aln2 =~ tr/-/-/;
	$missing_res_count2 = $native_full_aln =~ tr/-/-/;
	
	# STDOUT
	print "SEQUENCE ALIGNMENT: Full Native to Original Native:\n";
	print "NATIVE:       $native_aln2\n";
	print "FULL NATIVE:  $native_full_aln\n";
	print "\n";
	
	print "STRUCTURAL ALIGNMENT: Original Native to Design:\n";
	print "DESIGN:       $design_aln\n";
	print "NATIVE:       $native_aln1\n";
	print "COMMENT:      Sequence identity between design and native: $seq_ident\n";
	print "\n";

	if ($missing_res_count1 == 0) {
		# FASTA Sequence ID and Etc Information
		$design_file = basename($input_pdb);
		$design_seq = $design_aln;
		$design_seq =~ s/-//g;
		
		print FASTA ">$design_file\n";
		print FASTA "$design_seq\n";
		
		print "FULL DESIGN:  $design_seq\n";
		print "COMMENT:      No residues of missing densities detected.\n";
	}
	else {
		# FASTA Sequence ID and Etc Information
		$design_file = basename($input_pdb);
		print FASTA ">$design_file\n";
		
		# Insert missing residue densities into the design sequence!
		print "TRANSFERING RESIDUES OF MISSING DENSITIES TO DESIGN";
		
		# Output whether to insert or insert residues of missing densities to the N-terminal and/or C-terminal
		if ($insert_N eq 'TRUE' and $insert_C eq 'FALSE') {
			print " (Inserting residues of missing densities ONLY at the N-terminal)";
		}
		elsif ($insert_N eq 'FALSE' and $insert_C eq 'TRUE') {
			print " (Inserting residues of missing densities ONLY at the C-terminal)";
		}
		elsif ($insert_N eq 'TRUE' and $insert_C eq 'TRUE') {
			print " (Inserting residues of missing densities at the N-terminal and C-terminal)";
		}
		else {
			print " (Ignoring residues of missing densities at the N-terminal and C-terminal)";
		}
		print ":\n";
		($design_full_seq,$warning_ref) = insertMissingResDensity($design_aln,$native_aln1,$native_aln2,$native_full_aln,$insert_N,$insert_C);
				
		# New Design Sequence with Missing Residue Densities.
		print FASTA "$design_full_seq\n";
		
		print "FULL DESIGN:  $design_full_seq\n";
		
		$missing_count = $native_aln2 =~ tr/-/-/;
		print "COMMENT:      A total of $missing_count residues with missing densities detected and inserted into the design sequence.\n";
		
		# STDOUT Warnings
		if (scalar(keys %{$warning_ref}) > 0) {
			print "WARNING:      Residues of missing densities will be inserted in neighboring design sequence position:";
			for $pos_des (sort {$a <=> $b} keys %{$warning_ref}) {
				print " ", $pos_des + 1, " ";
			}
			print "\n";
		}
	}
	print "\n";
	
	# Clean sub temp directories
	print `rm -r $design_tmp`;
	print `rm -r $native_tmp`;
	
	# End of PDB
	print "-----------------------------------------------------------------------------\n";
	
}

# Clean temp directory
print `rm -r $tmp_dir`;

# Close FASTA file
close(FASTA);

# STDOUT
print "\n------------ COMPLETE! ------------\n";
print "\nFASTA formatted file of the full design sequence(s) outputted to: $output_file\n\n\n";


# ---- Subs ----

sub insertMissingResDensity {
	
	my $design_aln = shift;
	my $native_aln1 = shift;
	my $native_aln2 = shift;
	my $native_full_aln = shift;
	my $insert_N = shift || 'FALSE';
	my $insert_C = shift || 'FALSE';
	
	# Remove N' and/or C' missing densities (optional)
	($native_aln2,$native_full_aln) = remove_indel_N_term($native_aln2,$native_full_aln) if $insert_N eq 'FALSE';
	($native_aln2,$native_full_aln) = remove_indel_C_term($native_aln2,$native_full_aln) if $insert_C eq 'FALSE';
	
	my ($aa_full_nat,$aa_og_nat,$aa_nat,$aa_des,$aa_final);
	my ($pos_des,$pos_nat);
	my ($nat_id,$nat_chain);
	my %design_pos;
	my %warning;
	
	my ($aa1,$aa2,$pos,$pos_degap);
	my ($design_aln_degap,$native_aln1_degap);
	my %insert;
	
	my $design_full_seq;
	
	# Ensure that the native sequence is the same length as the design sequence
	# by temporarily removing extra design amino acids (caused by grafting)
	
	$pos = 0;
	for $aa1 (split //, $design_aln) {
		
		$aa2 = substr($native_aln1,$pos,1);
		
		# If we find a gapped position in the native sequence, temporarily remove and store it.
		# NOTE: In addition, we assume that we only keep gaps in the design (if it is shorter) to 
		# maintain equal length with the native sequence
		if ($aa2 eq '-') {
			$insert{'design'}{$pos_degap} .= $aa1;
			# print "\nTEST START: $insert{'design'}{$pos_degap} $pos_degap\n";
		}
		# Else, add all other residues to generate a degapped alignment
		else {
			$design_aln_degap .= $aa1;
			$native_aln1_degap .= $aa2;
			$pos_degap++;
		}
		
		# Keep track of design regions
		$design_pos{$pos}++ if $aa1 ne $aa2;
		
		# Increment the design position of the alignment
		$pos++;
		
	}
	
	# Check degapped alignments are equal in length
	die "ERROR: Degapped native and design alignments are not equal in length!" unless length($design_aln_degap) == length($native_aln1_degap);
	
	$pos_des = 0;
	$pos_nat = 0;
	
	# With respect to the full native sequence (with missing residue densities)
	for $aa_full_nat (split //, $native_full_aln) {
		
		# Get corresponding residue from the original native sequence (without missing residue densities)
		$aa_og_nat = substr($native_aln2,$pos_nat,1);
		
		# If missing residue densities are found for the current amino acid (gap in the native sequence)
		if ($aa_og_nat eq '-') {
			
			# Other formatting
			$aa_nat = '-';
			$aa_des = '-';
			
			# Set native residue of missing residue density to the final design sequence
			$aa_final = $aa_full_nat;
			
			# WARNING if insertion is next to a site of design (at the current position or left and/or right to it)
			if ($design_pos{$pos_des-1} or $design_pos{$pos_des} or $design_pos{$pos_des+1}) {
				$warning{$pos_des}++;
			}
			
		}
		else {
			# If no missing residues are found, correspond native sequence and design sequence with
			# their original identity (naturally, $aa_full_nat and $aa_nat1 remain the same)
			$aa_nat = substr($native_aln1_degap,$pos_des,1);
			$aa_des = substr($design_aln_degap,$pos_des,1);
			
			# Keep the identity of the design residue
			$aa_final = $aa_des;
			
			# Increment position of the design sequence if no residue of missing densities are found
			if ($aa_og_nat eq $aa_nat) {
				$pos_des++;
			}
			
		}
		
		# Increment position of the full native sequence
		$pos_nat++;
		
		# Skip if it's a insert of the original native residue (design sequence is shorter due to grafting)
		next if $insert{'native'}{$pos_des};
		# Build final design sequence with the inserted residue of missing densities (and insert gaps)
		# print "$aa_full_nat\t$aa_og_nat\t$aa_nat\t$aa_des\t->\t$aa_final\n";
		$design_full_seq .= $aa_final unless $aa_final eq '-';
		
		# Additional reinserting of residues added from grafting
		if ($insert{'design'}{$pos_des}) {
			for $aa_des (split //, $insert{'design'}{$pos_des}) {
				$aa_final = $aa_des;
				
				# print "-\t-\t-\t$aa_des\t->\t$aa_final\n";
				$design_full_seq .= $aa_final;
			}
		}
		
	}
	
	return ($design_full_seq,\%warning);
}

sub seqTransferInsert {
	# This subroutine will transfer missing residue inserts from sequence 1 to sequence 2
	# with the sequence alignment 1 and sequence alignment 2 respectively
	# Note: this assumes design sequence is shorter or equal than the actual PDB sequence
	my $seq_design_aln = shift;
	my $seq_actual_aln = shift;
	
	# Check if seq_design is shorter than seq_actual
	my $seq_design = $seq_design_aln;
	my $seq_design =~ s/-//;
	
	my $seq_actual = $seq_actual_aln;
	my $seq_actual =~ s/-//;
	
	my $seq_design_len = length($seq_design);
	my $seq_actual_len = length($seq_actual);
	
	die "Sequence 1 is not shorter or equal to Sequence 2\n" unless $seq_design_len <= $seq_actual_len;
	
	# Main
	my $pos_num = 0;
	my ($aa1,$aa2);
	my $seq;
	for $aa1 (split //, $seq_design_aln) {
		$aa2 = substr($seq_actual_aln,$pos_num,1);
		if ($aa1 eq '-') {
			# Add missing residue from actual sequence
			$seq .= $aa2;
		}
		else {
			# Add residue from the design sequence
			$seq .= $aa1;
		}
		
		$pos_num++;
	}
	
	return $seq;
}

sub remove_indel_N_term {
	my $seq_aln1 = shift;
	my $seq_aln2 = shift;
	
	# We assume that the first sequence has the N-temrinal gap(s).  If not ...
	# Reverse sequence variables when the first sequence alignment does not have N-terminal gap(s)
	my $seq_aln1_temp;
	my $seq_aln2_temp;
	
	if ($seq_aln1 !~ /^\-/ and $seq_aln2 =~ /^\-/) {
		$seq_aln1_temp = $seq_aln1;
		$seq_aln2_temp = $seq_aln2;
		
		$seq_aln1 = $seq_aln2_temp;
		$seq_aln2 = $seq_aln1_temp;
	}
	
	my $new_seq_aln1;
	my $new_seq_aln2;
	
	my $gaps_N_term;
	($gaps_N_term) = $seq_aln1 =~ /^(\-+)/;
	
	my $num_insert;
	
	# Evaluation the alignments for N-temrinal indels!
	if (!$gaps_N_term) {
		# Return the exact input alignment if no N-terminal gaps are found!
		return($seq_aln1,$seq_aln2);
	}
	else {
		# Remove Deletions (Gaps)
		$new_seq_aln1 = $seq_aln1;
		$new_seq_aln1 =~ s/^$gaps_N_term//;
		
		# Remove Insertions
		$num_insert = length($gaps_N_term);
		$new_seq_aln2 = $seq_aln2;
		$new_seq_aln2 =~ s/^\w{$num_insert}//;
		
		# Return alignment with no N-terminal indels!
		return($new_seq_aln1,$new_seq_aln2);
	}
	
}

sub remove_indel_C_term {
	my $seq_aln1 = shift;
	my $seq_aln2 = shift;
	
	# We assume that the first sequence has the C-temrinal gap(s).  If not ...
	# Reverse sequence variables when the first sequence alignment does not have C-terminal gap(s)
	my $seq_aln1_temp;
	my $seq_aln2_temp;
	
	if ($seq_aln1 !~ /\-$/ and $seq_aln2 =~ /\-$/) {
		$seq_aln1_temp = $seq_aln1;
		$seq_aln2_temp = $seq_aln2;
		
		$seq_aln1 = $seq_aln2_temp;
		$seq_aln2 = $seq_aln1_temp;
	}
	
	my $new_seq_aln1;
	my $new_seq_aln2;
	
	my $gaps_C_term;
	($gaps_C_term) = $seq_aln1 =~ /(\-+)$/;
	
	my $num_insert;
	
	# Evaluation the alignments for C-temrinal indels!
	if (!$gaps_C_term) {
		# Return the exact input alignment if no C-terminal gaps are found!
		return($seq_aln1,$seq_aln2);
	}
	else {
		# Remove Deletions (Gaps)
		$new_seq_aln1 = $seq_aln1;
		$new_seq_aln1 =~ s/$gaps_C_term$//;
		
		# Remove Insertions
		$num_insert = length($gaps_C_term);
		$new_seq_aln2 = $seq_aln2;
		$new_seq_aln2 =~ s/\w{$num_insert}$//;
		
		# Return alignment with no C-terminal indels!
		return($new_seq_aln1,$new_seq_aln2);
	}
	
}

sub pdbSeq_RESTful {
	
	use XML::Simple;
	use Data::Dumper;
	
	my $xml = new XML::Simple;
	
	my $pdb_id = shift;
	my $pdb_chain = shift;
	
	my $curl = '/usr/bin/curl';
	
	# Get PDB info for entire pdb entry
	my $data = `$curl http://www.pdb.org/pdb/rest/das/pdbchainfeatures/sequence?segment=$pdb_id.$pdb_chain 2> /dev/null`;
	my $info = $xml->XMLin($data);
	
	my $seq = $info->{'SEQUENCE'}{'content'};
	$seq =~ tr/\s\n//d;
	
	#print Dumper($info);
	
	return $seq;
	
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

sub readPDBAtom {
	
	my $pdb_file = shift;
	my $pdb_chain = shift;
	
	# -- Read PDB File

	my ($index,$atom,$aa,$mer,$res_num);
	my ($x,$y,$z);
	my ($occ,$bfact);
	my %pdb_atom;

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
			# Senenomethionine to Methionine Residue
			substr($_,17,3,"TYR");
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

sub parseChain {
	
	my $pdb_file = shift;
	my $chain_id = shift;
	my $out_dir = shift;
	
	my $mer;
	
	open(PDB,$pdb_file);
	open(CHAIN,">$out_dir/$chain_id.pdb");
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
			# Senenomethionine to Methionine Residue
			substr($_,17,3,"TYR");
		}
		
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

sub seq2atom {
	
	# This subroutine converts a sequence to a hash table of arbitrary atoms and coordinates
	
	my $pdb_seq = shift;
	my $mer = shift || 'A'; # Pre-defined chain ID
	my $index = shift || 1; # The start index
	my $res_num = shift || 1; # The start residue number
	
	my %aa_code = (
		 'A' => 'ALA',
		 'C' => 'CYS',
		 'D' => 'ASP',
		 'E' => 'GLU',
		 'F' => 'PHE',
		 'G' => 'GLY',
		 'H' => 'HIS',
		 'I' => 'ILE',
		 'K' => 'LYS',
		 'L' => 'LEU',
		 'M' => 'MET',
		 'N' => 'ASN',
		 'P' => 'PRO',
		 'Q' => 'GLN',
		 'R' => 'ARG',
		 'S' => 'SER',
		 'T' => 'THR',
		 'V' => 'VAL',
		 'W' => 'TRP',
		 'Y' => 'TYR'
	);
	
	my $aa;
	my %pdb_atom;
	my $atom;
	
	for $aa (split //, $pdb_seq) {
		
		# Note: only arbitrary backbone atoms per residue
		for $atom ('N','CA','C','O') {
			$pdb_atom{$index}{'atom'} = $atom;
			$pdb_atom{$index}{'aa'} = $aa_code{$aa};
			$pdb_atom{$index}{'mer'} = $mer;
			$pdb_atom{$index}{'res_num'} = $res_num;
			$pdb_atom{$index}{'x'} = 0;
			$pdb_atom{$index}{'y'} = 0;
			$pdb_atom{$index}{'z'} = 0;
			$pdb_atom{$index}{'occ'} = 1;
			$pdb_atom{$index}{'bfact'} = 0;
			
			$index++;
		}
		
		$res_num++;
	}
	
	return %pdb_atom;
	
}

sub cutPDBAtoms {
	my $pdb_atom_ref = shift;
	my $start_res_cut = shift;
	my $stop_res_cut = shift;
	
	my %pdb_atom = %{$pdb_atom_ref};
	
	my ($index,$atom,$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
	
	my $out;
	my $res_count++;
	
	my %pdb_atom_shorter = %pdb_atom;
	
	for my $index (sort {$a <=> $b} keys %pdb_atom) {
		
		$res_num = $pdb_atom{$index}{'res_num'};
		
		if ($res_num >= $start_res_cut and $res_num <= $stop_res_cut) {
			delete $pdb_atom_shorter{$index};
		}
		
	}
	
	return %pdb_atom_shorter;
}

sub pastePDBAtoms {
	my $pdb_atom_ref = shift; # Original PDB Atoms
	my $pdb_atom_paste_ref = shift; # PDB Atoms to paste
	my $start_res_paste = shift; # Residue number to start pasting
	
	# Residue position before the sequence paste
	my $start_res_paste_before = $start_res_paste - 1;
	
	my %pdb_atom = %{$pdb_atom_ref};
	my %pdb_atom_paste = %{$pdb_atom_paste_ref};
	
	my ($index,$atom,$aa,$mer,$res_num,$x,$y,$z,$occ,$bfact);
	my $res_num_rem;
	my $aa_rem;
	
	my %pdb_atom_longer;
	
	my $index_paste;
	my $index_new;
	my $paste_rem;
	
	for my $index (sort {$a <=> $b} keys %pdb_atom) {
		
		$res_num = $pdb_atom{$index}{'res_num'};
		
		if ($paste_rem) {
			# Use incremented index after the residue paste!
			$index_new++;
		}
		else {
			# Use regular index
			$index_new = $index;
		}
		
		if ($res_num > $start_res_paste_before and !$paste_rem) {
			
			$res_num = $res_num_rem;
			
			for my $index_paste (sort {$a <=> $b} keys %pdb_atom_paste) {
				$index_new++;
				
				$res_num++ if $aa_rem ne $pdb_atom_paste{$index_paste}{'aa'};
				
				$pdb_atom_longer{$index_new}{'atom'} = $pdb_atom_paste{$index_paste}{'atom'};
				$pdb_atom_longer{$index_new}{'aa'} = $pdb_atom_paste{$index_paste}{'aa'};
				$pdb_atom_longer{$index_new}{'mer'} = $pdb_atom_paste{$index_paste}{'mer'};
				$pdb_atom_longer{$index_new}{'res_num'} = $res_num;
				$pdb_atom_longer{$index_new}{'x'} = $pdb_atom_paste{$index_paste}{'x'};
				$pdb_atom_longer{$index_new}{'y'} = $pdb_atom_paste{$index_paste}{'y'};
				$pdb_atom_longer{$index_new}{'z'} = $pdb_atom_paste{$index_paste}{'z'};
				$pdb_atom_longer{$index_new}{'occ'} = $pdb_atom_paste{$index_paste}{'occ'};
				$pdb_atom_longer{$index_new}{'bfact'} = $pdb_atom_paste{$index_paste}{'bfact'};
				
				$aa_rem = $pdb_atom_paste{$index_paste}{'aa'};
			}
			
			$paste_rem++;
		}
		else {
			$pdb_atom_longer{$index_new}{'atom'} = $pdb_atom{$index}{'atom'};
			$pdb_atom_longer{$index_new}{'aa'} = $pdb_atom{$index}{'aa'};
			$pdb_atom_longer{$index_new}{'mer'} = $pdb_atom{$index}{'mer'};
			$pdb_atom_longer{$index_new}{'res_num'} = $pdb_atom{$index}{'res_num'};
			$pdb_atom_longer{$index_new}{'x'} = $pdb_atom{$index}{'x'};
			$pdb_atom_longer{$index_new}{'y'} = $pdb_atom{$index}{'y'};
			$pdb_atom_longer{$index_new}{'z'} = $pdb_atom{$index}{'z'};
			$pdb_atom_longer{$index_new}{'occ'} = $pdb_atom{$index}{'occ'};
			$pdb_atom_longer{$index_new}{'bfact'} = $pdb_atom{$index}{'bfact'};
		}
		
		$res_num_rem = $res_num;
		
	}
	
	return %pdb_atom_longer;
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

sub align_all {
	# Use MUSCLE for Multiple Sequence Alignment
	
	my $seqs_ref = shift;
	my $matrix = shift;
	
	my %seqs = %{$seqs_ref};
		
	my $muscle_exe = "$ENV{'HOME'}/bin/muscle";
	my $temp_dir = "/tmp/muscle_$$";
	
	# Make temp dir
	mkdir($temp_dir);
	
	my $in_file = "$temp_dir/in.txt";
	my $out_file = "$temp_dir/out.txt";
	
	# Generate Temperary FASTA formatted file for MUSCLE alignment
	open(FILE,">$in_file") or die "Cannot open file $in_file\n";
	# Print out all sequences
	print FILE ">$_\n$seqs{$_}\n" for (keys %seqs);
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
	my %seqs_aln = fasta2hash($out_file);
	
	# Clean temp dir
	rmdir($temp_dir);
	
	return (%seqs_aln);
	
}

sub TMalign {
	# Use TMalign for structural alignment
	
	my $pdb1 = shift;
	my $pdb2 = shift;
	
	my $TMalign_path = "/work/davela/bin/TMalign";
	
	my $line;
	my $start;
	
	my ($seq_aln1,$seq_aln2);
	
	for $line (`$TMalign_path $pdb1 $pdb2`) {
		$start++ if $line =~ /denotes aligned residue pairs/;
		
		if ($start) {
			next if $line =~ /\:/;
			
			if (!$seq_aln1) {
				$seq_aln1 = $line;
			}
			else {
				$seq_aln2 = $line;
				$start = 0;
			}
		}
		
	}
	
	chomp $seq_aln1;
	chomp $seq_aln2;
	
	return($seq_aln1,$seq_aln2);
	
}

sub id_seq {
        # Find sequence identity
        my $seq1 = shift;
        my $seq2 = shift;

        die "Alignments not equal in length:\nSEQ1: $seq1\n\nSEQ2: $seq2\n" unless length($seq1) == length($seq2);

        my @seq1 = split //, $seq1;
        my @seq2 = split //, $seq2;

        my $id;
        my $total;

        my $i;
        for $i (0..$#seq1) {

                # Skip gaps
                next if $seq1[$i] eq '-';
                next if $seq2[$i] eq '-';

                # Count identities
                $id++ if $seq1[$i] eq $seq2[$i];
                $total++;

        }

        my $ident;
        if ($total == 0) {
                # For all gaps and nothing compared?
                $ident = 0;
        }
        else {
                $ident = $id / $total;
        }

        return $ident;
}

sub fetchPDB {
	# Fetch PDB File
	my $pdb_id = shift;
	my $pdb_chain = shift;
	my $tmp_dir = shift;
	
	# Download PDB File
	my $pdb_out = `curl http://www.rcsb.org/pdb/files/$pdb_id.pdb 2> /dev/null`;
	die "Cannot find retreive PDB: $pdb_id\n" unless $pdb_out =~ /ATOM/;
	
	# Place PDB info in temp directory
	my $pdb_path = "$tmp_dir/$pdb_id$pdb_chain.pdb";
	write_file($pdb_out,$pdb_path);
	
	# Parse chain if it was specified
	if ($pdb_chain) {
		$pdb_path = parseChain($pdb_path,$pdb_chain,$tmp_dir);
	}
	
	return $pdb_path;
} 

sub write_file {
	my $text = shift;
	my $out_file = shift;
	
	open(DATA,">$out_file") or die "Cannot open file $out_file\n";
	print DATA "$text";
	close(DATA);
	
}
