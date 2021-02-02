#!/usr/bin/perl

use strict;
use File::Basename;

my $usage = "pdb2seq.pl <pdb_file> <pdb_chain>\n";

my $file = $ARGV[0] or die $usage;
my $chain = $ARGV[1] || 'NA';


# ---- Main ----

my $base = basename($file);
my %pdb = readPDBAtom($file,$chain);
my $seq = atom2seq(\%pdb);


print ">$base | $chain\n";
print "$seq";
print "\n";

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
		
		# End loop if MODEL 2 found (For other NMR structures)
		# We use only the first NMR model by default
		last if /^ENDMDL/;


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

