#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::Genotypes;

# written by Alexander Nater, October 2015

# Set your options here:
my $workdir="/home/ubuntu/output/sumstats";
my $fasta_file="/data/data/refgenom/PonAbe_7x_callable.fasta";
my $fasta_index=$fasta_file . ".fai";
my $allsitesfolder="/data/data/vcf/allSites";
my $allsiteslist="$workdir/allsiteslist_7pops.txt";
my $vcffolder="/data/data/vcf/SNPs";
my $vcflist="$workdir/vcflist_7pops.txt";
my $npopulations=7;
my @poplabel=("NKSK", "EK", "WKCK", "SAR", "BT", "NALK", "WA");
my @nindividuals=(4, 3, 9, 4, 2, 8, 6);
my @minind=(2, 2, 4, 2, 2, 4, 4);
my @popgroups=( [0, 1, 2, 3], [4, 5, 6] );
my $maxsize=1000000000;
my $minlength=100000;	# minimum length of chromosome/scaffold to be considered.
my $minmapingquality=20;
my $minvariantquality=0;
my $mingenotypequality=0;	# this setting is questionable, as it biases towards lower variability!
my $mincoverage=5;
my $MAFfilter=0;
my $excludeRefN=1;	# FALSE=0/TRUE=1 - exclude sites with N in reference sequence.
my $excludenonbiallelic=1;	# FALSE=0/TRUE=1 - exclude sites with more than two alleles in population pair.
my $haplotize=0;	# FALSE=0/TRUE=1 - only take one random allele per individual.
my $minpropcovered=0.3;	# minimum proportion of valid sites per window.
my $outfolder=shift @ARGV;
my $outprefix=shift @ARGV;
my $chromosome=shift @ARGV;
my $windowsize=shift @ARGV;
my $stepsize=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

$outfolder=~s{/\z}{};	# remove trailing slash from folder path.
mkdir $outfolder unless (-d "$outfolder");
my $outfile="$outprefix" . "_" . "$chromosome" . ".bed";

# prepare reference fasta file:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_index, $fasta_file);
my $regions=$fasta->locifromIndex();
my $windows=$regions->windowingLoci($minlength, $windowsize, $stepsize, 0, $chromosome);

# prepare bed files of regions to be excluded:
# my $genes=Regions->new();
# $genes->readBED($genes_file, 'onebased', 0);

# regions from file:
# my $windows=Regions->new();
# $windows->readBED($regions_file, 'onebased', 1);


SUBSET: while (1){
	my $subset=$windows->subsetLoci($maxsize);
	$subset->printLoci("-");
	unless ( $subset->getSize() ){ last SUBSET }
	my $concat=$subset->concatRegions();
	$concat->addRefseq($fasta);
	$concat->printLoci("-");
	my $missing=Missing->new($allsitesfolder, $allsiteslist);
	$missing->readAllsites($concat, $minmapingquality, $minvariantquality, \@nindividuals, $excludeRefN);
	my $genotypes=Genotypes->new($vcffolder, $vcflist);
	$genotypes->readvcfList($concat, $mingenotypequality, 0, 1);
# 	$genotypes->readOutgroups($outgroupfolder, $outgroupvcflist, $concat);
	$subset->calculateSSWithinPops($genotypes, $missing, \@popgroups, $mincoverage, $minpropcovered, \@minind, $haplotize, $MAFfilter, $excludenonbiallelic);
#	$windows->storeSFS($concat, $genotypes, $missing, \@popgroups, $mincoverage, \@minind, $haplotize);
	$subset->printLociSSWithinPops("$outfolder/$outfile", 1);
	}



