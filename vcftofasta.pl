#!/usr/bin/env perl

use strict;
use warnings;
use Classes::Fasta;
use Classes::Regions;
use Classes::Missing;
use Classes::Genotypes;
use Classes::Misc;

# written by Alexander Nater, January 2015

# Set your options here:
my $workdir="/home/ubuntu/output/vcftofasta";
my $fasta_file="/data/data/refgenom/PonAbe_7x_callable.fasta";
my $fasta_index=$fasta_file . ".fai";
my $allsitesfolder="/data/data/vcf/allSites";
my $allsiteslist="$workdir/allsiteslist.txt";
my $vcffolder="/data/data/vcf/SNPs";
my $vcflist="$workdir/vcflist.txt";
my $npopulations=7;
my @poplabels=("NKSK", "EK", "WKCK", "SAR", "BT", "NALK", "WA");
my @nindividuals=(4, 3, 9, 4, 2, 8, 7);
my @minind=(0, 0, 0, 0, 0, 0, 0);
my @popgroups=( [0, 1, 2, 3], [4, 5, 6] );
my @grouplabels=("Borneo", "Sumatra");
my @mingroupind=(8, 8);
my $maxsize=2000000;
my $rowlength=100;
my $minmapingquality=20;
my $minvariantquality=5;
my $mingenotypequality=20;
my $outfolder=shift @ARGV;
my $outprefix=shift @ARGV;
my $mincoverage=shift @ARGV;
my $excluderefN=shift @ARGV;	# FALSE=0/TRUE=1
my $haplotize=shift @ARGV;	# FALSE=0/TRUE=1
my $idlist=shift @ARGV;
my @subpoplabels=split(',', shift @ARGV);
my $subindstring=shift @ARGV;

#-----------------------------------------------------------------------------------------------------------------------------

$outfolder=~s{/\z}{};	# remove trailing slash from folder path.
mkdir "$outfolder" unless (-d "$outfolder");

# prepare new arrays for population settings if subset of populations is selected:
my @subpopindices;
if (@subpoplabels && $subpoplabels[0] ne 'all'){ @subpopindices=Misc::setSubArrays(\@subpoplabels, \@poplabels, \@nindividuals, \@minind, \@popgroups, \@grouplabels, \@mingroupind) }
else { @subpopindices=(0..$npopulations-1) }
print "Selected population indices: ", join(',', @subpopindices), "\nNumber of individuals per population: ", join(',', @nindividuals), "\nMinimum number of individuals per population: ", join(',', @minind), "\n";
print "$grouplabels[$_]: ", join(',', @{ $popgroups[$_] }), "\n" foreach (0..$#popgroups);

# obtain array of individual indices to be printed:
my ($ref_indnames, $ref_files, $ref_samples)=Misc::getIndnames($allsitesfolder, $allsiteslist);
my %subindlabels;
%subindlabels=map { $_=>1 } split(',', $subindstring) if (defined $subindstring);
my (@subindividuals, @indnames);
my $totsize=0;
for my $pop (@subpopindices){
	my @tempind;
	if (defined $subindstring){
		for my $index (0..$#{ $ref_indnames->[$pop] }){
			push @tempind, $index if ($subindlabels{ $ref_indnames->[$pop][$index] });	# test if list of subpoplabels contain individual ids.
			}
		}
	else { @tempind=( 0..$#{ $ref_indnames->[$pop] } ) }
	push @subindividuals, \@tempind;
	push @indnames, map { $ref_indnames->[$pop][$_] } @tempind;
	$totsize+=@tempind;
	}
print STDERR join(',', @indnames), "\n";
die "Wrong number of individual names!\n" unless (@indnames==$totsize);

# prepare list of chromosomes/scaffolds to be processed:
my @scaffolds=Misc::parseArglist($idlist);

# prepare reference fasta index:
my $fasta=Fasta->new();
$fasta->readIndex($fasta_index, $fasta_file);

# prepare windows to process:
my $regions=$fasta->locifromIndex();
my $windows=$regions->windowingLoci(0, $maxsize, $maxsize);

# prepare fasta object for each individual to be printed:
my @fasta_objects;
push @fasta_objects, Fasta->new() foreach (0..$totsize-1);


SUBSET: while (1){
	my $subset=$windows->subsetLoci($maxsize, @scaffolds);
	$subset->printLoci("-");
	unless ( $subset->getSize() ){ last SUBSET }
	my $missing=Missing->new($allsitesfolder, $allsiteslist);
	$missing->setPops(\@subpopindices);
	$missing->readAllsites($subset, $minmapingquality, $minvariantquality, \@nindividuals);
	my $genotypes=Genotypes->new($vcffolder, $vcflist);
	$genotypes->setPops(\@subpopindices);
	$genotypes->readvcfList($subset, $mingenotypequality, 0, 1, 0);
	$subset->generateFasta($fasta, \@fasta_objects, $genotypes, $missing, \@minind, \@popgroups, \@mingroupind, $mincoverage, $excluderefN, $haplotize, \@subindividuals);
	}

for my $ind (0..$totsize-1){
	my $indfolder=$indnames[$ind];
	mkdir "/$outfolder/$indfolder" unless (-d "/$outfolder/$indfolder");
	my $outname="$outfolder" . "/" . "$indfolder" . "/" . "$outprefix" . "_" . "$indnames[$ind]" . "_" . $scaffolds[0] . "-" . $scaffolds[-1] . ".fasta";
	$fasta_objects[$ind]->printSeqs($outname, $rowlength);
	}


