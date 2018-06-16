#!/usr/bin/perl

# Written by David La
# Updated Sat Sep 20 18:37:33 PDT 2014

# Description:
# Converts resfile from ssm2res.pl to combinatorial libraries that can be 
# ordered as IDT overlaping ultramer fragments (max 200 nucleotides each)
# if you are using the default PCR-based method!

# Note: 
# Adds 40 nucleotide pETCON 5' (with ATG in front of your gene) and 3' flanks
# to each side of your gene!
#

my $usage = "Usage: res2lib.pl -resfile <resfile> -dna <wt_dna_seq> 
		  -method [PCR|kunkel (default: PCR)]
		  -max_fragment_length [max_fragment_length for the PCR Method (default: 200)]
		  -overlap_fragment_length [overlap_fragment_length for the PCR Method (default: 20)]
		  -primer_length [primer_length for the Kunkel Method (default: 60)]
		  -primer_extend [extend_length for the Kunkel Method (default: 10)]
		  -offset [number_of_positions_to_shift]\n";

use strict;
use Getopt::Long;

# ---- Special Arguments ----

my %args=();
GetOptions(	"resfile:s"=>\$args{resfile},
		"dna:s"=>\$args{dna},
		"method:s"=>\$args{method},
		"max_fragment_length:s"=>\$args{max_fragment_length},
		"overlap_fragment_length:s"=>\$args{overlap_fragment_length},
		"primer_length:s"=>\$args{primer_length},
		"primer_extend:s"=>\$args{primer_extend},
		"offset:n"=>\$args{offset});

# ---- Define Arguments ----

my $resfile = $args{resfile} or die $usage; # Resfile (needs the commented output from ssm2res.pl)
my $wt_dna = $args{dna} or die $usage; # WT DNA Sequence
my $method = $args{method} || 'PCR'; # Experimental Method
my $max_frag_len = $args{max_fragment_length} || 200; # Maximum Fragment Length for PCR Method
my $overlap_frag_len = $args{overlap_frag_length} || 20; # Overlap Fragment Length
my $primer_len = $args{primer_length} || 60; # Primer Length for Kunkel Mutagensis Method
my $primer_extend = $args{primer_extend} || 10; # Primer Extend Length for Kunkel Mutagensis Method
my $shift_pos = $args{offset}; # Increment or Decrement Position (if you lengthen or shorten the gene a priori)

# ---- Main ----

# Get Degenerate Codons from Resfile
my %degen_codon = res2codon($resfile);

# Make Sure DNA Sequence is All Uppercase!
$wt_dna = uc($wt_dna);

# Remove possible ATG at the beginning of the WT DNA Sequence, because ATG is added later as part of the 5' pETCON flank!
$wt_dna =~ s/^ATG// if $wt_dna =~ /^ATG/;

# Translate WT DNA sequence to protein
my $wt_protein = dna2protein($wt_dna);

# Select Experimental Method
my $degen_dna;
my $degen_dna_pETCON;
my $wt_dna_pETCON;
my %frag;
my %frag_method;


if ($method =~ /pc?r?/i) {
	
	# Generate Combinatorial Library (Degenerate DNA)
	$degen_dna = dna2lib($wt_dna,\%degen_codon,$shift_pos);

	# Add 5' and 3' pETCON Flanks
	$degen_dna_pETCON = add_flank_pETCON($degen_dna);
	
	# Generate Overlaping Fragments (Fragment Max Length=200, Fragment Overlap Length=20)
	%frag = dna2frag($degen_dna_pETCON,$max_frag_len,$overlap_frag_len);
	
	# Check if there are More Fragments to be Made for Additional Codons
	if ($degen_dna =~ /\+\+\+/) {
		# If there are mutliple codons at one position marked by "+++"
	%frag = frag_expand(\%frag,\%degen_codon);
	}

	# PCR Fragment Assembly: Generate Alternating Forward and Reverse Fragments for Overlap PCR
	%frag_method = gen_fwd_rev(\%frag);

}
elsif ($method =~ /ku?n?k?e?l?/i) {
	
	# Generate Short Fragment Primers
	%frag = dna2primers($wt_dna,\%degen_codon,$primer_len,$primer_extend,$shift_pos);
	
	# Expand Multiple Codons
	#%frag = frag_expand(\%frag,\%degen_codon);
	
	# Kunkel Mutagenesis: Generate Reverse Fragment Primers
	#%frag_method = %frag;
	%frag_method = gen_rev(\%frag);
}
else {
	
	# Exit if the Defined Method is Not Recognized!
	die "\nMethod not recognized: $method\n";
}

# Output in a Convenient IDT Cut-and-Paste Ready Format
print "\n";
out_IDT(\%frag_method);
print "\n";

# Output Number of Fragments
my $frag_num;
$frag_num++ for keys %frag_method;
print STDERR "Number of Fragments: $frag_num\n\n";


# ---- Subs ----

sub res2codon {
	my $resfile = shift;
	
	# Read SSM2RES Codon Output
	open(FILE,"$resfile") or die "Cannot open file $resfile\n";

	my ($codon,$pos);
	my %degen_codon;
	while (<FILE>) {
		chomp;
		next unless /^[\d\t]/;
		($codon,$pos) = $_ =~ /CODON=([\w\,]+)\s+POS=(\d+)/;
		$degen_codon{$pos} = $codon;
	}
	close(FILE);
	
	return %degen_codon;
}

sub dna2protein {
	my $dna_seq = shift;
	
	my %genetic_code = (
	 'TCA'=>'S', # Serine
	 'TCC'=>'S', # Serine
	 'TCG'=>'S', # Serine
	 'TCT'=>'S', # Serine
	 'TTC'=>'F', # Phenylalanine
	 'TTT'=>'F', # Phenylalanine
	 'TTA'=>'L', # Leucine
	 'TTG'=>'L', # Leucine
	 'TAC'=>'Y', # Tyrosine
	 'TAT'=>'Y', # Tyrosine
	 'TAA'=>'*', # Stop
	 'TAG'=>'*', # Stop
	 'TGC'=>'C', # Cysteine
	 'TGT'=>'C', # Cysteine
	 'TGA'=>'*', # Stop
	 'TGG'=>'W', # Tryptophan
	 'CTA'=>'L', # Leucine
	 'CTC'=>'L', # Leucine
	 'CTG'=>'L', # Leucine
	 'CTT'=>'L', # Leucine
	 'CCA'=>'P', # Proline
	 'CAT'=>'H', # Histidine
	 'CAA'=>'Q', # Glutamine
	 'CAG'=>'Q', # Glutamine
	 'CGA'=>'R', # Arginine
	 'CGC'=>'R', # Arginine
	 'CGG'=>'R', # Arginine
	 'CGT'=>'R', # Arginine
	 'ATA'=>'I', # Isoleucine
	 'ATC'=>'I', # Isoleucine
	 'ATT'=>'I', # Isoleucine
	 'ATG'=>'M', # Methionine
	 'ACA'=>'T', # Threonine
	 'ACC'=>'T', # Threonine
	 'ACG'=>'T', # Threonine
	 'ACT'=>'T', # Threonine
	 'AAC'=>'N', # Asparagine
	 'AAT'=>'N', # Asparagine
	 'AAA'=>'K', # Lysine
	 'AAG'=>'K', # Lysine
	 'AGC'=>'S', # Serine
	 'AGT'=>'S', # Serine
	 'AGA'=>'R', # Arginine
	 'AGG'=>'R', # Arginine
	 'CCC'=>'P', # Proline
	 'CCG'=>'P', # Proline
	 'CCT'=>'P', # Proline
	 'CAC'=>'H', # Histidine
	 'GTA'=>'V', # Valine
	 'GTC'=>'V', # Valine
	 'GTG'=>'V', # Valine
	 'GTT'=>'V', # Valine
	 'GCA'=>'A', # Alanine
	 'GCC'=>'A', # Alanine
	 'GCG'=>'A', # Alanine
	 'GCT'=>'A', # Alanine
	 'GAC'=>'D', # Aspartic Acid
	 'GAT'=>'D', # Aspartic Acid
	 'GAA'=>'E', # Glutamic Acid
	 'GAG'=>'E', # Glutamic Acid
	 'GGA'=>'G', # Glycine
	 'GGC'=>'G', # Glycine
	 'GGG'=>'G', # Glycine
	 'GGT'=>'G', # Glycine
	 );
	
	my $protein_seq;
	my $codon;
	
	# Split sequence by 3s (codons)
	for $codon (split /(\w{3})/, $dna_seq) {
		# Skip empty strings
		next if !$codon;
		# Concatentate converted protein amino acids
		$protein_seq .= $genetic_code{$codon};
	}
	
	return $protein_seq;
}

sub dna_complement {
	my $dna_seq = shift;
	
	my %dna_complem = (
		'A' => 'T',
		'T' => 'A',
		'G' => 'C',
		'C' => 'G',
		'R' => 'Y',
		'Y' => 'R',
		'M' => 'K',
		'K' => 'M',
		'S' => 'S',
		'W' => 'W',
		'H' => 'D',
		'D' => 'H',
		'V' => 'B',
		'B' => 'V',
		'N' => 'N',
		'+' => '+'
	);
	
	my $nuc;
	my $nuc_complem;
	my $dna_complem;
	for $nuc (split //, $dna_seq) {
		$nuc_complem = $dna_complem{$nuc};
		$dna_complem .= $nuc_complem;
	}
	
	return $dna_complem;
}

sub dna2lib {
	# Generate Library for Positions with Single Degenerate Codons
	# and Marks positions with multiple codons using "+++"
	
	my $wt_dna = shift;
	my $degen_codon_ref = shift;
	my $shift_pos = shift;
	
	my %degen_codon = %{$degen_codon_ref};
	
	# Initial WT Sequence
	my $degen_dna = $wt_dna;
	my $codon;
	my $first_codon;
	my $aa_pos;
	my $dna_pos;
	
	for my $pos (keys %degen_codon) {
		# Get Degenerate Codon
		$codon = $degen_codon{$pos};
		
		# Mark Codon if there are Multiple
		if ($codon =~ /\,/) {
			$codon = "+++"; # Temporary Marker
		}
		
		# Shift Position by Defined Number
		$aa_pos = $pos + $shift_pos;
		
		# Find DNA position
		$dna_pos = $aa_pos * 3; # Codon per AA
		$dna_pos = $dna_pos - 2; # First Position of Codon
		$dna_pos = $dna_pos - 1; # Offset for Internal Numbering
		
		# Replace WT Sequence Codon with Degenerate Codon
		#print "$codon $aa_pos $dna_pos\n"; # TEST OUTPUT
		substr($degen_dna,$dna_pos,3,$codon);
		
	}
	
	return $degen_dna;
}

sub dna2mark {
	# Marks positions with all codons with "+++"
	
	my $wt_dna = shift;
	my $degen_codon_ref = shift;
	my $shift_pos = shift;
	
	my %degen_codon = %{$degen_codon_ref};
	
	# Initial WT Sequence
	my $degen_dna = $wt_dna;
	my $codon;
	my $first_codon;
	my $aa_pos;
	my $dna_pos;
	
	for my $pos (keys %degen_codon) {
		# Get Degenerate Codon
		$codon = $degen_codon{$pos};
		
		# Mark Codon if there are Multiple
		if ($codon =~ /\,/) {
			$codon = "+++"; # Multiple Codons
		}
		else {
			$codon = '???' # Single Codon
		}
		
		# Shift Position by Defined Number
		$aa_pos = $pos + $shift_pos;
		
		# Find DNA position
		$dna_pos = $aa_pos * 3; # Codon per AA
		$dna_pos = $dna_pos - 2; # First Position of Codon
		
		# Replace WT Sequence Codon with Degenerate Codon
		substr($degen_dna,$dna_pos,3,$codon);
		
	}
	
	return $degen_dna;
}

sub count_mark {
	my $mark_seq = shift;
	
	my $mark_count;
	$mark_count += () = $mark_seq =~ /\?\?\?/g; # Count number of '???' single codons
	$mark_count += () = $mark_seq =~ /\+\+\+/g; # Count number of '+++' multiple codons
	
	return $mark_count;
}

sub determine_extend {
	# Determine how much to extend the 5' and/or 3' ends!
	
	my $dna_seq = shift; # We expected marked codons!
	my $extend_len = shift || 10; # 10 nucleotides between mutations and ends
	
	my $five_prime_dna;
	my $three_prime_dna;
	
	my $five_prime_len;
	my $three_prime_len;
	
	my $five_prime_extend = 0;
	my $three_prime_extend = 0;
	
	if ($dna_seq =~ /[\?\+]/) {
		($five_prime_dna) = $dna_seq =~ /^(\w+)[\?\+]/;
		($three_prime_dna) = $dna_seq =~ /[\?\+](\w+)$/;
		
		$five_prime_len = length($five_prime_dna);
		$three_prime_len = length($three_prime_dna);
		
		# Evaluate 5' End
		if ($five_prime_len < 10 ) {
			$five_prime_extend = $extend_len - $five_prime_len;
		}
		else {
			$five_prime_extend = 0;
		}
		
		# Evaluate 3' End
		if ($three_prime_len < 10) {
			$three_prime_extend = $extend_len - $three_prime_len;
			
		}
		else {
			$three_prime_extend = 0;
		}
	}
	
	return ($five_prime_extend,$three_prime_extend);
}

sub dna2frag {
	my $dna = shift;
	my $max_len = shift || 200; # Max for IDT Ultramer
	my $overlap_len = shift || 20; # Overlap Between Fragments
	
	# For Each Nucleotide in the Fragment
	my $len;
	my %frag;
	my $num = 1;
	my $overlap_dna;
	my $overlap_no_degen;
	my $overlap_remain;
	for my $nuc (split //, $dna) {
		# Keep track of DNA Fragment Length
		$len++;
		
		# If the Max Fragment Length is Reached
		if ($len > $max_len) {
			# Overlap for Next Fragment (Tail Segment of the Fragment)
			$overlap_dna = substr($frag{$num},-1*$overlap_len,$overlap_len);
			$num++;
			
			# Check if Overlap DNA has Degenerate Nucleotide Positions
			if ($overlap_dna =~ /[RYMKSWHDVBN\+]/) {
				
				# Slide Overlap Window Backwards to Find Overlapping Regions without Degenerate Nucleotide Positions
				for my $pos (1..$max_len-$overlap_len) {
					# Check Previous Fragment
					$overlap_no_degen = substr($frag{$num-1},-1*$overlap_len-$pos,$overlap_len);
					
					# Found Overlapping Window without Degenerate Nucleotides
					if ($overlap_no_degen !~ /[RYMKSWHDVBN\+]/) {
						
						# Replace Overlapping Region of the Previous Fragment with new Non-Degenerate Overlapping Region
						$frag{$num-1} =~ s/($overlap_no_degen.*)/$overlap_no_degen/;
						
						# Remember Non-Degenerate Overlap + Degenerate Overlap
						$overlap_remain = $1;
						last;
					}
				}
				
			}
			
			# If a DNA Overlap Degeneracy was Evaluated and Another Overlap of Non-Degeneracy was Found
			if ($overlap_remain) {
				# Start New Fragment with Non-Degenerate Overlap + Degenerate Overlap
				$frag{$num} = $overlap_remain;
				$len = length($overlap_remain) + 1;
				$overlap_remain = ''; # Reset if used for Evaluating Next the Fragment Overlap
			}
			# Else a DNA Overlap Degeneracy was NOT Evaluated or was Evaluated but Overlap of Non-Degeneracy was NOT Found
			else {
				# Incase DNA Degeneracy was Evaluated and No Overlapping Regions of Non-Degeneracy was Found
				print STDERR "\nNote: Could not prevent DNA overlapping Fragment " . ($num-1) . " and Fragment " . $num . " from containing degenerate codon(s)!  Try making the overlap sequence length between these fragments shorter!\n" if $overlap_dna =~ /[RYMKSWHDVBN]/;
				
				# Start New Fragment with Overlap as Regularly
				$frag{$num} = $overlap_dna;
				$len = $overlap_len + 1;
			}
			
		}
		
		# Add Nucleotide to Fragment
		$frag{$num} .= $nuc;
	}
	
	return %frag;
}

sub dna2primers {
	# Generate Primer Fragment Library for Positions with Single Degenerate Codons
	
	my $wt_dna = shift;
	my $degen_codon_ref = shift;
	my $primer_base_len = shift || 60; # The standard length of each primer
	my $primer_extend_len = shift || 10; # How much to extend if ends has mutations
	my $shift_pos = shift;
	
	my %degen_codon = %{$degen_codon_ref};
	
	# pETCON Flanks
	my $five_prime = "TGGAGGCGGTAGCGGAGGCGGAGGGTCGGCTAGCCATATG"; # Note, this adds Methionine to the 5'
	my $three_prime = "CTCGAGGGAGGCGGATCCGAACAAAAGCTTATTTCTGAAG";
	
	# Initial WT Sequence
	my $degen_dna = $wt_dna;
	my $degen_dna_full;
	my $degen_dna_mark_full;
	my $primer_pos_full;
	my $codon;
	my $codon_single;
	my $first_codon;
	my $aa_pos;
	my $dna_pos;
	my $primer_pos;
	my $primer_count;
	my $sub_primer_count;
	my %primer_lib;
	my $skip_rem;
	my $primer_seq;
	my $mark_count;
	my $mark_seq;
	my $five_ext;
	my $three_ext;
	my $primer_pos_rem;
	
	# Add Entire WT DNA Sequence (For Reference)
	$primer_lib{'WT'} = $wt_dna;
	
	# Determine DNA Length
	my $dna_len = length($degen_dna);
	
	# Determine Flanks Lengths
	my $flank_len = ($primer_len - 3) / 2;
	
	# Round-off to Nearest Integer
	$flank_len = sprintf("%.0f",$flank_len);
	
	# Generate Mark-up Degenerate DNA ('???' for single codons and '+++' for multliple codons)
	my $degen_dna_mark = dna2mark($wt_dna,$degen_codon_ref,$shift_pos);
	$degen_dna = $degen_dna_mark;
	
	#print "\nTEST: $degen_dna\n\n"; # TEST
	
	# All Defined Codon Positions
	for my $pos (sort {$a <=> $b} keys %degen_codon) {
		
		# Get Degenerate Codon
		$codon = $degen_codon{$pos};
		
		# Shift Position by Defined Number
		$aa_pos = $pos + $shift_pos;
		
		# Find DNA position
		$dna_pos = $aa_pos * 3; # Codon per AA
		$dna_pos = $dna_pos - 2; # First Position of Codon
		
		# Mark Codon if there are Multiple
		if ($codon =~ /\,/) {
			$codon = "+++"; # Temporary Marker
		}
		
		# Replace WT Sequence Codon with Degenerate Codon
		substr($degen_dna,$dna_pos,3,$codon);
		
		# Skip future codons already observed to avoid primer sequence overlaps
		if ($skip_rem > 0) {
			# Keep Multiple '+++' Mark Codon and Skip
			if ($codon =~ /\+\+\+/) {
				# Decrement Skips Number
				$skip_rem--;
				
				# Skip!
				next;
			}
			# Replace Single Codon '???' Mark with Actual Codon and Skip
			else {
				# Edit same primer
				$primer_lib{"Primer\.$primer_count\_$mark_count"} =~ s/[\?][\?][\?]/$codon/;
				
				# Remember Last DNA position of Last codon
				$skip_rem--;
				
				# Skip!
				next;
			}
			
		}
		
		# Primer Counter
		$primer_count++;
		
		# Reformat Numbers
		$primer_count = sprintf("%03d",$primer_count);
		
		# Position of DNA to Start Parsing the Primer Out from Degenerate DNA
		$primer_pos = $dna_pos - $flank_len;
		
		# Check if the Previous Position will Overlap with the Current Position
		if ($primer_pos_rem+$primer_base_len >= $primer_pos && $primer_pos_rem) {
			# Start primer from left to right, instead from the middle!
			$primer_seq = substr($degen_dna,$dna_pos,$primer_base_len);
			$mark_seq = substr($degen_dna_mark,$dna_pos,$primer_base_len);
		}
		else {
			# Build the Library
			if ($dna_pos <= $flank_len) {
				# Add 5' pETCON flank if mutations are close to 5'
				$degen_dna_full = "$five_prime$degen_dna";
				$degen_dna_mark_full = "$five_prime$degen_dna_mark";
				$primer_pos_full = $primer_pos + 40;
			}
			elsif ($dna_pos <= $dna_len+$flank_len) {
				# Add 3' pETCON flank if mutations are close to 3'
				$degen_dna_full = "$degen_dna$three_prime";
				$degen_dna_mark_full = "$degen_dna_mark$three_prime";
				$primer_pos_full = $primer_pos;
			}
			else {
				# Determine Primer Normally
				$degen_dna_full = $degen_dna;
				$degen_dna_mark_full = $degen_dna_mark;
				$primer_pos_full = $primer_pos;
			}
			
			# Get the Primers
			$primer_seq = substr($degen_dna_full,$primer_pos_full,$primer_base_len);
			$mark_seq = substr($degen_dna_mark_full,$primer_pos_full,$primer_base_len);
			
			# Extend when necessary (when mutations are near the 5' or 3' ends)
			#($five_ext,$three_ext) = determine_extend($mark_seq,$primer_extend_len);
			#if ($five_ext > 0 or $three_ext > 0) {
			#	$primer_seq = substr("$five_prime$degen_dna",$primer_pos+40+$five_ext,$primer_base_len+$three_ext);
			#}
			
		}
		
		# Count number of mutation via codons marks
		$mark_count = count_mark($mark_seq);
		
		# Number of times to skip based on number of marks observed
		$skip_rem = $mark_count - 1;
		
		# Store Primer
		$primer_lib{"Primer\.$primer_count\_$mark_count"} = $primer_seq;
		
		# Remember Previous DNA Position of Last Codon
		$primer_pos_rem = $primer_pos;
		
	}
	
	#print "\nTEST: $degen_dna\n\n"; # TEST
	
	return %primer_lib;
}

sub frag_expand {
	
	my $frag_ref = shift;
	my $degen_codon_ref = shift;
	
	my %frag = %{$frag_ref};
	my %degen_codon = %{$degen_codon_ref};
	
	my $dna_frag;
	my $og_id;
	my $codon_list;
	my $codon;
	my %expand;
	my %expand_new;
	my $format_id;
	my $sub_frag_count = 0;
	my $start;
	my $num_multi_codon;
	my $pos_count;
	my %pos_rem;
	
	for my $id (sort {$a <=> $b} keys %frag) {

		#print "FRAG: $id $frag{$id}\n";
		
		# DNA Fragment Sequence
		$dna_frag = $frag{$id};
		
		# New Fragment Codon Position
		$start = 1;

		# Count number of positions with multiple codons
		$num_multi_codon = 0;
		$num_multi_codon++ for split /\+\+\+/, $dna_frag;
		$num_multi_codon--;
		$pos_count = 0;
		#print "POS COUNT: $num_multi_codon\n";
		
		# If We Find "+++" Positions
		if ($dna_frag =~ /\+\+\+/) {
			
			# Replace "+++" with List of Codons
			for my $pos (sort {$a <=> $b} keys %degen_codon) {

				# Codons (Comma Delimited)
				$codon_list = $degen_codon{$pos};
				
				# Skip position if it only has a single codon
				next unless $codon_list =~ /\,/;
				#print "NOW POS: $pos\n";

				# End loop if no more '+++' (multi codon) positions to evalauate
				$pos_count++ unless $pos_rem{$pos};
				#print "POS COUNTER: $pos_count\n";
				last if $pos_count > $num_multi_codon;

				# Skip previous positions already evaluated
				#print "BEFORE CURRENT POS: $pos\n";
				next if $pos_rem{$pos};
				##print "SKIP: $pos\n" && next if $pos_rem{$pos};
				#print "AFTER CURRENT POS: $pos\n";

				# Remember previous positions
				#print "POS REM: $pos\n";
				$pos_rem{$pos}++;

				# Evaluate Each Codon
				for $codon (split /\,/, $codon_list) {
					
					# First Codon Position for this Fragment
					unless (!$start) {
						
						# Increment Fragment Count for Next Combination Fragment
						$sub_frag_count++;
						$sub_frag_count = sprintf("%03d",$sub_frag_count);
						
						# Get Initial Fragment
						$expand_new{"$id.$sub_frag_count"} = $dna_frag;
						#print qq(FIRST TEST1: $expand_new{"$id.$sub_frag_count"}\n);
						
						# Search and Replace Only First Match
						$expand_new{"$id.$sub_frag_count"} =~ s/\+\+\+/$codon/;
						#print qq(FIRST TEST2: $expand_new{"$id.$sub_frag_count"}\n);
						
						#print "FIRST: $codon ($codon_list)\t$pos\t$id.$sub_frag_count |\n";
					}

					# Remaining Codon Positions for this Fragment
					else {
						
						for my $expand_id (sort {$a <=> $b} keys %expand) {
							
							# Skip Other Fragments (Processed Previously)
							($og_id) = $expand_id =~ /(\w+)\./;
							next unless $og_id eq $id;
							
							# Increment Fragment Count for Next Combination Fragment
							$sub_frag_count++;
							$sub_frag_count = sprintf("%03d",$sub_frag_count);
							
							# Get Previous Fragment (if it exists)
							$dna_frag = $expand{$expand_id};
							$expand_new{"$og_id.$sub_frag_count"} = $dna_frag;
							#print qq(OTHER TEST1: $expand_new{"$og_id.$sub_frag_count"}\n);
							
							# Search and Replace Only First Match
							$expand_new{"$og_id.$sub_frag_count"} =~ s/\+\+\+/$codon/;
							#print qq(OTHER TEST2: $expand_new{"$og_id.$sub_frag_count"}\n);
							
							#print "OTHER: $codon ($codon_list)\t$pos\t$id.$sub_frag_count\n";
						}
					}
					
				}

				# Add or update %expand_new to the main %expand
				for my $new_id (keys %expand_new) {
					$expand{$new_id} = $expand_new{$new_id};
				}
				
				# Reset %expand_new
				%expand_new = ();

				# Not the First Codon Position for this Fragment After this Point!
				$start = 0;

				# Reset Fragment Counter
				$sub_frag_count = 0;
				
			}

		}
		# Add all Other Fragments
		else {
			$expand{$id} = $dna_frag;
		}
		
		# Reset New Expand Fragments
		%expand_new = ();

	}
	
	return %expand;
}

sub gen_fwd_rev {
	
	# Generate Alternating Forward (FWD) and Reverse (REV) Fragments for Overlap PCR
	
	my $frag_ref = shift;
	
	my %frag = %{$frag_ref};
	
	# Generate Forward and Reverse Fragments
	my %frag_pcr;
	my $count;
	my $sub_id;
	for my $id (sort {$a <=> $b} keys %frag) {
		if ($id =~ /\d+\.\d+/) {
			# Prevent Counting Expanded Fragments
			($sub_id) = $id =~ /\d+\.(\d+)/;
			$count++ unless $sub_id > 1;
		}
		else {
			$count++;
		}
		# Reverse Complement of Even Numbered Fragments
		if ($count % 2 == 0) {
			$frag_pcr{"Fragment\_$id\_REV"} = reverse dna_complement($frag{$id});
		}
		else {
			$frag_pcr{"Fragment\_$id\_FWD"} = $frag{$id};
		}
	}
	
	return %frag_pcr;
	
}

sub gen_rev {
	
	# Generate ONLY Reverse (REV) Complement Fragments
	
	my $frag_ref = shift;
	
	my %frag = %{$frag_ref};
	
	# Generate Reverse Complement DNA Fragments
	my %frag_rev;
	my $count;
	my $sub_id;
	for my $id (sort {$a <=> $b} keys %frag) {
		if ($id =~ /\d+\.\d+/) {
			# Prevent Counting Expanded Fragments
			($sub_id) = $id =~ /\d+\.(\d+)/;
			$count++ unless $sub_id > 1;
		}
		else {
			$count++;
		}
		
		# Include WT
		if ($id =~ /WT/) {
			$frag_rev{"WT\_FWD"} = $frag{$id};
		}
		# Add all other Fragments
		else {
			$frag_rev{"Fragment\_$id\_REV"} = reverse dna_complement($frag{$id});
		}
	}
	
	return %frag_rev;
	
}

sub add_flank_pETCON {
	my $dna = shift;
	
	my $five_prime = "TGGAGGCGGTAGCGGAGGCGGAGGGTCGGCTAGCCATATG"; # Note, this adds Methionine to the 5'
	my $three_prime = "CTCGAGGGAGGCGGATCCGAACAAAAGCTTATTTCTGAAG";
	
	my $dna_pETCON = $five_prime . $dna . $three_prime;
	
	return $dna_pETCON;
}

sub out_IDT {
	my $seq_hash_ref = shift;
	
	my %seq_hash = %{$seq_hash_ref};
	
	print "ID\tSequence\n";
	for my $id (sort {$a cmp $b} keys %seq_hash) {
		if ($id =~ /WT/) {
			print "\n$id\t$seq_hash{$id}\n";
		}
		else {
			print "$id\t$seq_hash{$id}\n";
		}
	}
}
