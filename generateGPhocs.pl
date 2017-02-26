#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::Genotypes;
use Classes::Nexus;
use Classes::Misc;

# written by Alexander Nater, October 2016
if (@ARGV < 6){
	die "Wrong number of arguments!\nUsage: ./generateGPhocs.pl [output file] [comma separated list of chromosome names or 'all' if working over whole genome] [windowsize] [stepsize] [comma separated list of individual ids] [verbose output 1=yes/0=no]\n";
	}
print STDERR $_+1, ". Argument: ", $ARGV[$_], "\n" foreach (0..$#ARGV);

# Set your options here:
my $workdir="/home/ubuntu/output/GPhocs";
my $fasta_file="/data/data/refgenom/PonAbe_7x_callable_mostconservedmasked.fasta";
my $fasta_index=$fasta_file . ".fai";
my $genes_file="/data/data/refgenom/masking/pongo_genic.bed";
my $allsitesfolder="/data/data/vcf/allSites";
my $allsiteslist="$workdir/allsiteslist_5groups.txt";
my $vcffolder="/data/data/vcf/SNPs";
my $vcflist="$workdir/vcflist_5groups.txt";
my $npopulations=5;
my @poplabels=("PP", "ST", "NT", "HS", "PT");
my @nindividuals=(20, 2, 14, 2, 2);
my @minind=(1, 1, 1, 1, 1);
my @popgroups=( [0], [1], [2], [3], [4] );
my @grouplabels=("PP", "ST", "NT", "HS", "PT");
my $maxsize=10000000;	# combined maximum size of regions to process together.
my $mincoverage=5;
my $minmapingquality=20;
my $minvariantquality=0;
my $mingenotypequality=0;	# this setting biases towards lower diversity if > 0!
my $minprop_valid=0.7;	# minimum proportion of valid sites per locus.

my $outfile=shift @ARGV;
my $regions_file=shift @ARGV;
my @ids=Misc::parseArglist(shift @ARGV);
my $windowsize=shift @ARGV;
my $stepsize=shift @ARGV;
my $subindstring=shift @ARGV;
my $verbose=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------


# obtain array of individual indices and array of subpop labels:
my ($ref_indnames, $ref_files, $ref_samples)=Misc::getIndnames($allsitesfolder, $allsiteslist);
unless (@$ref_files==@nindividuals){ die "Wrong number of populations in popsizes array!\n" }
my ($ref_subsamples, $ref_subinds, @subpoplabels);
if ($subindstring eq 'all'){ $ref_subsamples=$ref_samples; $ref_subinds=$ref_indnames; @subpoplabels=@poplabels }
else {
	($ref_subsamples, $ref_subinds)=Misc::getIndicesfromIndnames($ref_indnames, $ref_samples, [ split(',', $subindstring) ] );
	for my $pop (0..$#{ $ref_subsamples }){
		push @subpoplabels, $poplabels[$pop] if (@{ $ref_subsamples->[$pop] });
		$nindividuals[$pop]=@{ $ref_subsamples->[$pop] };
		}
	}
print STDERR join(',', @subpoplabels), "\n";

# prepare new arrays for population settings for subset of populations with selected individuals:
my @subpopindices=Misc::setSubArrays(\@subpoplabels, \@poplabels, \@nindividuals, \@minind, \@popgroups, \@grouplabels);
print "Selected population indices: ", join(',', @subpopindices), "\nNumber of individuals per population: ", join(',', @nindividuals), "\nMinimum number of individuals per population: ", join(',', @minind), "\n";
print "$grouplabels[$_]: ", join(',', @{ $popgroups[$_] }), "\n" foreach (0..$#popgroups);
print STDERR join(',', @nindividuals), "\n";

# prepare reference fasta file:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_index, $fasta_file, 0);

# prepare windows:
my $regions;
if (defined $regions_file && $regions_file ne 'None'){
	$regions=Regions->new();
	$regions->readBED($regions_file, 'bed', 0);
	}
else { $regions=$fasta->locifromIndex() }
my $windows=$regions->windowingLoci($windowsize, $windowsize, $stepsize, 0, @ids);

# exclude windows near genic regions:
my $genes=Regions->new();
$genes->readBED($genes_file, 'bed', 0);
$windows=$windows->excludeRegions($genes, 1000);
print STDERR $windows->getSize($_), " loci remaining after filtering genic regions on chromosome $_\n" foreach (@ids);

# prepare nexus object:
my $nexus=Nexus->new();


SUBSET: while (1){
	# get a subset of loci with a combined size smaller than maxsize:
	my $subset=$windows->subsetLoci($maxsize);
	unless ( $subset->getSize() ){ last SUBSET }
	$subset->printLoci("-");

	# acquire missing data:
	my $missing=Missing->new($allsitesfolder, $allsiteslist);
	$missing->setPops(\@subpopindices, $ref_subsamples);
	$missing->readAllsites($subset, $minmapingquality, $minvariantquality, \@nindividuals);

	# acquire genotype data:
	my $genotypes=Genotypes->new($vcffolder, $vcflist);
	$genotypes->setPops(\@subpopindices, $ref_subsamples);
	$genotypes->readvcfList($subset, $mingenotypequality, 0, 1, 0, 0, 0, $verbose);

	# generate GPhocs input file:
	$nexus->genotypestoGPhocs($outfile, $subset, $fasta, $genotypes, $missing, \@minind, $mincoverage, $minprop_valid, 0);
	}


