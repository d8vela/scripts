#!/usr/bin/perl

# Written by David La
# Updated: Fri Jun 27 04:34:14 PDT 2014

# Description:
# This script will get information about your PDB structure 
# by getting the data directly from the PDB!

my $usage = "pdbInfo.pl <pdb_id> <pdb_chain>\n";

use strict;
use lib '/work/davela/perl5/lib/perl5';

my $pdb_id = $ARGV[0] or die $usage;
my $pdb_chain = $ARGV[1];

# Check if the pdb_id argument was combined with pdb_chain
if (!$pdb_chain) {
	if (length($pdb_id) != 5) {
		die $usage;
	}
	else {
		# Get the last char
		$pdb_chain = chop $pdb_id;
	}
}

# Prep chain
$pdb_id =~ tr/a-z/A-Z/;
$pdb_chain =~ tr/a-z/A-Z/;

# Get the info
my ($info1,$info2) = pdbInfo("$pdb_id","$pdb_chain");
die "PDB ID: $pdb_id-$pdb_chain not found!\n" unless $$info1{'PDB'};

# The info about the entire PDB entry
my $nr_atoms = $$info1{'PDB'}{'nr_atoms'};
my $status = $$info1{'PDB'}{'status'};
my $citation_authors = $$info1{'PDB'}{'citation_authors'};
my $revision_date = $$info1{'PDB'}{'revision_date'};
my $nr_residues = $$info1{'PDB'}{'nr_residues'};
my $keywords = $$info1{'PDB'}{'keywords'};
my $expMethod = $$info1{'PDB'}{'expMethod'};
my $publish_date = $$info1{'PDB'}{'publish_date'};
my $structure_authors = $$info1{'PDB'}{'structure_authors'};
my $structureId = $$info1{'PDB'}{'structureId'};
my $resolution = $$info1{'PDB'}{'resolution'};
my $title = $$info1{'PDB'}{'title'};
my $nr_entities = $$info1{'PDB'}{'nr_entities'};

# The info about the specific PDB chain
my $chainId = $$info2{'structureId'}{'chainId'};
my $Taxonomy_name = $$info2{'structureId'}{'polymer'}{'Taxonomy'}{'name'};
my $Taxonomy_id = $$info2{'structureId'}{'polymer'}{'Taxonomy'}{'id'};

my $macroMolecule_name = $$info2{'structureId'}{'polymer'}{'macroMolecule'}{'name'};
my $macroMolecule_accession_id = $$info2{'structureId'}{'polymer'}{'macroMolecule'}{'accession'}{'id'};

my $length = $$info2{'structureId'}{'polymer'}{'length'};
my $entityNr = $$info2{'structureId'}{'polymer'}{'entityNr'};
my $weight = $$info2{'structureId'}{'polymer'}{'weight'};
my $polymerDescription = $$info2{'structureId'}{'polymer'}{'polymerDescription'}{'description'};

my $type = $$info2{'structureId'}{'polymer'}{'type'};
my $id = $$info2{'structureId'}{'id'};


# Output
print "$id\t$chainId\t$type\t$polymerDescription\t$length\t$resolution\t$expMethod\t$Taxonomy_name\t$title\t$keywords\n";


# ---- Subs ----

sub pdbInfo {
	
	use XML::Simple;
	#use Data::Dumper;
	
	my $xml = new XML::Simple;
	
	my $pdb_id = shift;
	my $pdb_chain = shift;
	
	my $curl = '/usr/bin/curl';
	
	# Get PDB info for entire pdb entry
	my $data1 = `$curl http://www.pdb.org/pdb/rest/describePDB?structureId=$pdb_id 2> /dev/null`;
	my $info1 = $xml->XMLin($data1);
	
	# Get PDB info for specific chain
	my $data2 = `$curl http://www.pdb.org/pdb/rest/describeMol?structureId=$pdb_id.$pdb_chain 2> /dev/null`;
	my $info2 = $xml->XMLin($data2);
	
	
	#print Dumper($info1);
	#print Dumper($info2);
	
	
	return ($info1,$info2);
	
}


