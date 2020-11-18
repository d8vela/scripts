#!/usr/bin/perl

# Written by David La
# Updated: Fri Aug 17 12:38:15 PDT 2018

# Description:
# This script will renumber a PDB file from first residue number of the original chain to the last residue number for all chains
# Note: Used for renumbering chothia pdb numbering formats for RMSD calcuations (input to rmsd_ab.pl).

my $usage = "Usage: renumPDB.pl <pdb_file>\n";

use strict;

my $pdb_file = $ARGV[0] or die $usage;

my %chains = getChains($pdb_file);
my %chain_range = getChainRange($pdb_file);

my %pdb_atom;
my $out;

for my $mer (sort {$a cmp $b} keys %chains) {
	# Read PDB File (this also renumbers)
	%pdb_atom = readPDBAtom($pdb_file,$mer,\%chain_range);
	
	# Output the reformatted PDB Atoms
	$out = outputPDBAtom(\%pdb_atom);
	print $out;
}


# ---- Subs ----

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
	my $chain_range_ref = shift;

	%chain_range = %{$chain_range_ref};
	
	# -- Read PDB File

	my ($index,$atom,$aa,$mer,$res_num);
	my ($x,$y,$z);
	my ($occ,$bfact);
	my %pdb_atom;

	my $mer_rem;
	my $res_num_rem;
	my $renum;

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
		$renum = $chain_range{$mer}{'begin'} - 1 if $mer ne $mer_rem;
		$renum++ if $res_num ne $res_num_rem;

		$mer_rem = $mer;
		$res_num_rem = $res_num;

		$res_num = $renum;
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

