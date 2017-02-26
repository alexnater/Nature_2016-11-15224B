#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::Genotypes;

# written by Alexander Nater, September 2016

# Set your options here:
my $workdir="/home/ubuntu/output/genomescans";
my $fasta_file="/data/data/refgenom/PonAbe_7x_callable.fasta";
my $fasta_index=$fasta_file . ".fai";
my $allsitesfolder="/data/data/vcf/allSites";
my $allsiteslist="$workdir/allsiteslist_3groups.txt";
my $vcffolder="/data/data/vcf/SNPs";
my $vcflist="$workdir/vcflist_3groups.txt";
my $npopulations=5;
my @poplabels=("PP", "ST", "NT", "HS", "PT");
my @nindividuals=(20, 2, 14, 2, 2);
my @minind=(2, 2, 2, 1, 1);
my @popgroups=( [0], [1], [2], [3,4] );
my @grouplabels=("PP", "ST", "NT", "Outgroups");
my $maxsize=1000000000;
my $minlength=100000;	# minimum length of chromosome/scaffold to be considered.
my $minmapingquality=20;
my $minvariantquality=0;
my $mingenotypequality=0;	# this setting is questionable, as it biases towards lower variability!
my $mincoverage=5;
my $MAFfilter=0;
my $excludeRefN=0;	# FALSE=0/TRUE=1 - exclude sites with N in reference sequence.
my $excludenonbiallelic=1;	# FALSE=0/TRUE=1 - exclude sites with more than two alleles in population pair.
my $haplotize=0;	# FALSE=0/TRUE=1 - only take one random allele per individual.
my $minpropcovered=0.3;	# minimum proportion of valid sites per window.
my $outfolder=shift @ARGV;
my $outprefix=shift @ARGV;
my $chromosome=shift @ARGV;
my $subindstring=shift @ARGV;
my $windowsize=shift @ARGV;
my $stepsize=shift @ARGV;
my $verbose=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

$outfolder=~s{/\z}{};	# remove trailing slash from folder path.
mkdir $outfolder unless (-d "$outfolder");
my $outfile="$outprefix" . "_" . "$chromosome" . ".bed";

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
	$concat->addRefseq($fasta) if $excludeRefN;
	$concat->printLoci("-");

	# acquire missing data:
	my $missing=Missing->new($allsitesfolder, $allsiteslist);
	$missing->setPops(\@subpopindices, $ref_subsamples);
	$missing->readAllsites($concat, $minmapingquality, $minvariantquality, \@nindividuals, $excludeRefN, $verbose);

	# acquire genotype data:
	my $genotypes=Genotypes->new($vcffolder, $vcflist);
	$genotypes->setPops(\@subpopindices, $ref_subsamples);
	$genotypes->readvcfList($concat, $mingenotypequality, 0, 0, 0, 0, 0, $verbose);

	# calculate and print summary statistics:
# 	$genotypes->readOutgroups($outgroupfolder, $outgroupvcflist, $concat);
	$subset->calculateSS($genotypes, $missing, \@popgroups, $mincoverage, $minpropcovered, \@minind, $haplotize, $MAFfilter, $excludenonbiallelic);
#	$windows->storeSFS($concat, $genotypes, $missing, \@popgroups, $mincoverage, \@minind, $haplotize);
	$subset->printLociSS("$outfolder/$outfile", 1);
	}



