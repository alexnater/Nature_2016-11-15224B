package Outgroups;

use strict;
use warnings;
use Classes::Misc;



sub new{
	my $class=shift;
	my $frequencies=shift;
	my $alleles=shift;
	my $self={};
	$self->{_Frequencies}=defined $frequencies ? $frequencies : {};
	$self->{_Alleles}=defined $alleles ? $alleles : {};
	bless $self, $class;
	return $self;
	}


sub getValue{
	my ($self, $chrom, $pos, $pop)=@_;
	if (@_==3){ return($self->{_Frequencies}{$chrom}{$pos}) }
	elsif (@_==4){ return($self->{_Frequencies}{$chrom}{$pos}[$pop]) }
	else { die "Invalid number of arguments!\n" }
	}


sub setValue{
	unless (@_==5){die "Wrong number of arguments!\n"}
	my ($self, $chrom, $pos, $pop, $value)=@_;
	$self->{_Frequencies}{$chrom}{$pos}[$pop]=$value;
	return($value);
	}


sub getAncestral{
	unless (@_==3){ die "Wrong number of arguments!\n" }
	my ($self, $id, $pos)=@_;
	return($self->{_Ancestral}{$id}{$pos});
	}


sub getBase{
	my ($self, $chrom, $pos, $i, $j)=@_;
	if (@_==4){
		if ($i eq '.'){ return("N") } 
		elsif ($i>=0 && $i<2){ return( substr($self->{_Bases}{$chrom}{$pos}, $i, 1) ) }
		elsif ($i>=3 && $i<4){ return("N") }
		else { return }
		}
	else { die "Invalid number of arguments!\n" }
	}


sub readMAFs{
	my ($self, $mafsfolder, $mafslist, $ref_loci, $minfreq, $recordbases, $storagelength)=@_;

	my ($mafsfiles, $samples)=Misc::openvcflist($mafsfolder, $mafslist);

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			for my $pop (0..$#{$mafsfiles}){
				open my $mafsfile, "-|", "tabix $mafsfiles->[$pop] $chrom:$start-$end" or die "Could not open vcf file $mafsfiles->[$pop]!\n$!\n";
				print "Aquiring frequencies for $chrom:$start-$end, population $pop.\n";
				my ($scaffold, $prevscaff, $startpos);
				LINE: while (<$mafsfile>){
					if(/^\w+\s/){
						my @line=split(/\s/, $_);
						$scaffold=$line[0];
						my $pos=$line[1]-1;	# The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
						my $nsamples=defined $line[6] ? $line[6] : 0;
						$nsamples=255 if $nsamples>255;
						my $offset=$pos % $storagelength;
#						print "$pop - $scaffold:$pos - $offset\n";
						if (!defined $prevscaff || $scaffold ne $prevscaff || $pos-$startpos>=$storagelength){
							$startpos=$pos-$offset;
							$self->{_NSamples}{$chrom}{$startpos}[$pop]=(chr(0) x $storagelength);
							}
						if ($recordbases){
							warn "Major base not the same as ancestral base at $scaffold:$pos!\n" if $line[2] ne $line[4];
							$self->{_Bases}{$scaffold}{$pos}="$line[2]" . "$line[3]";
							$self->{_Ancestral}{$scaffold}{$pos}=$line[4];
							}
						substr( $self->{_NSamples}{$chrom}{$startpos}[$pop], $offset, 1, chr($nsamples) );
						$self->{_Frequencies}{$scaffold}{$pos}[$pop]=$line[5] unless ($line[5]<$minfreq);
						}
					} continue { $prevscaff=$scaffold }
				close $mafsfile;
				}
			}
		}
	return;
	}



sub readOutgroups{	# under construction!
	my ($self, $ref_loci, $allsitesfolder, $allsiteslist, $minmq, $minvq, $mingtq, $suballeles, $recordbases)=@_;

	my %missing;	# $missing{chrom/scaff}[locus][allind]="char string";

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my $indsum=0;
			my $iterator;
			if ($suballeles){ $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus) }
			for my $pop (0..$#{$files}){
				open my $vcffile, "-|", "tabix $files->[$pop] $chrom:$start-$end" or die "Could not open vcf file $files->[$pop]!\n$!\n";
				print "Aquiring coverage and genotype data for $chrom:$start-$end, population $pop.\n";
				my $prevpos=$start-2;
				my ($pos, $scaffold, $snp, $dp);
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1]-1;	# vcf files are 1-based, internal storage is 0-based.
						my $vq=$line[5];

						my $gapsize=$pos-$prevpos-1;
						if ($gapsize){
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=(chr(0) x $gapsize) }
							$mq{$scaffold}[$locus].=(chr(0) x $gapsize) if ($recordmq);
							}
						$prevpos=$pos;

						my ($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
						if ($recordmq){
							if (defined $mq && $mq<255){ $mq{$scaffold}[$locus].=chr( int($mq) ) }
							elsif (defined $mq){ $mq{$scaffold}[$locus].=chr(255) }
							else { $mq{$scaffold}[$locus].=chr(0) }
							}
						unless ($vq ne '.' && $vq>=$minvq){
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=chr(0) }
							next LINE;
							}
						unless (defined $mq && $mq>=$minmq){
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=chr(0) }
							next LINE
							}

						if ($line[4] eq '.'){
							$snp=0; $dp=1;
							}
						else {
							$snp=1; $dp=2;
							$genotypes->{_Alleles}{$scaffold}{$pos}=0 unless defined $genotypes->{_Alleles}{$scaffold}{$pos};
							unless ( $pop && exists $genotypes->{_Genotypes}{$scaffold}{$pos} ){
								for my $i (0..$npops-1){ $genotypes->{_Genotypes}{$scaffold}{$pos}[$i]='00' x $popsize->[$i] }
								}
							$genotypes->{_Genotypes}{$scaffold}{$pos}[$pop]='';
							if ($recordbases){ $genotypes->{_Bases}{$scaffold}{$pos}="$line[3]" . join('', split(',' , $line[4]) ) }
							}

						my $allind=$indsum;
						my $index=0;
						for my $ind (@{$samples->[$pop]}){
							my @data=split(":", $line[$ind+8]);
							if (@data>1 && $data[$dp]<255){ $missing{$scaffold}[$locus][$allind].=chr($data[$dp]) }
							elsif (@data>1){ $missing{$scaffold}[$locus][$allind].=chr(255) }
							else { $missing{$scaffold}[$locus][$allind].=chr(0) }
							if ($snp){
								if ($data[0]=~/([0-3])\/([0-3])/){
									if ( $data[3]>=$mingtq){ $genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=$1 . $2 }
									else { $genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=".." }
									if ($suballeles && defined $iterator->[$pop][$index]){	# only records allele states for subsampled individuals.
										if (2*$ind==$iterator->[$pop][$index]){$genotypes->setAlleles($scaffold, $pos, $1); ++$index}
										if (2*$ind+1==$iterator->[$pop][$index]){$genotypes->setAlleles($scaffold, $pos, $2); ++$index}
										}
									else { $genotypes->setAlleles($scaffold, $pos, $1, $2) }
									}
								else {$genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=".."}
								}
							++$allind;
							}
						}
					}
				close $vcffile;
				my $gapsize=$end-$pos-1;
				if ($recordmq && $gapsize){ $mq{$scaffold}[$locus].=(chr(0) x $gapsize) }
				for my $ind (0..$popsize->[$pop]-1){
					if ($gapsize){ $missing{$chrom}[$locus][$indsum+$ind].=(chr(0) x $gapsize) }
					my $stringlength=length($missing{$chrom}[$locus][$indsum+$ind]);
					if ($stringlength != $locuslength){
						print "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				$indsum+=$popsize->[$pop];
				}
			}
		}
	$self->{_Missing}=\%missing;
	$self->{_Pops}=$popsize;
	$self->{_Indnames}=\@indnames;
	$self->{_MQ}=\%mq;
	return($genotypes);
	}


sub extractRegions{
	my ($self, $regions)=@_;

	my %reduced;

	CHROM: for my $chrom (sort $regions->allKeys()){
		if (defined $self->{_Frequencies}{$chrom}){
			REG: for my $reg ($regions->allValues($chrom)){
				POSITION: for my $pos (sort {$a<=>$b} keys %{$self->{_Frequencies}{$chrom}}){
					if ($pos>=$reg->{'start'} && $pos<=$reg->{'end'}){$reduced{$chrom}{$pos}=$self->{_Frequencies}{$chrom}{$pos}}
					elsif ($pos>$reg->{'end'}){next REG}
					}
				}
			}
		else {print "No frequencies available for chromosome $chrom!\n"}
		}

	my $newfrequencies=Frequencies->new(\%reduced);
	return($newfrequencies);
	}


sub translateScafftoChrom{
	my ($self, $chrommap)=@_;	# bed file format with chromosome locations

	my %tr_frequencies;
	my %tr_biallelic;

	SCAFF: while (my ($scaff, $ref_scaff)=each %{$self->{_Frequencies}}){
		CHROM: for my $chrom ($chrommap->allKeys()){
			if (defined $chrommap->getValue($chrom, $scaff)){
				my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
				if ($dir eq '-'){
					my $offset=$chrommap->getValue($chrom, $scaff, 'end');
					POSITION: while (my ($pos, $ref_pos)=each %{$ref_scaff}){
						$tr_frequencies{$chrom}{$offset-$pos}=$ref_pos;
						$tr_biallelic{$chrom}{$offset-$pos}=$self->{_Biallelic}{$scaff}{$pos};
						}
					}
				else {
					my $offset=$chrommap->getValue($chrom, $scaff, 'start');
					POSITION: while (my ($pos, $ref_pos)=each %{$ref_scaff}){
						$tr_frequencies{$chrom}{$pos+$offset}=$ref_pos;
						$tr_biallelic{$chrom}{$pos+$offset}=$self->{_Biallelic}{$scaff}{$pos};
						}
					}
				if ($dir ne '+' && $dir ne '-') {print "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n"}
				next SCAFF;
				}
			}
		print "Warning: scaffold $scaff has not been located on a chromosome! Loci/Positions on this scaffold are excluded in the output!\n";
		}

	my $tr_loci=Frequencies->new(\%tr_frequencies, \%tr_biallelic, $self->{_Header}, $self->{_Samples});
	return($tr_loci);
	}


sub translateChromtoScaff{
	my ($self, $chrommap)=@_;	# bed file format with chromosome locations

	my %tr_frequencies;
	my %tr_biallelic;

	CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_Frequencies}}){
		my ($lastscaff, $scaffstart, $scaffend);
		POSITION: for my $pos (sort {$a<=>$b} keys %{$ref_chrom}){
			unless (defined $lastscaff && $pos<=$scaffend){
				SCAFF: for my $scaff ($chrommap->allKeys($chrom)){
					if ($pos>=$chrommap->getValue($chrom, $scaff, 'start') && $pos<=$chrommap->getValue($chrom, $scaff, 'end')){
						$lastscaff=$scaff;
						$scaffstart=$chrommap->getValue($chrom, $lastscaff, 'start');
						$scaffend=$chrommap->getValue($chrom, $lastscaff, 'end');
						last SCAFF;
						}
					}
				unless (defined $lastscaff && $pos<=$scaffend){
					print "Warning: Chromosome $chrom position $pos has not been located on a scaffold! This position is not included in the output!\n";
					next POSITION;
					}
				my $dir=$chrommap->getValue($chrom, $lastscaff, 'dir');
				if ($dir eq '-'){
					$tr_frequencies{$lastscaff}{$scaffend-$pos}=$ref_chrom->{$pos};
					$tr_biallelic{$lastscaff}{$scaffend-$pos}=$self->{_Biallelic}{$chrom}{$pos};
					}
				else {
					$tr_frequencies{$lastscaff}{$pos-$scaffstart}=$ref_chrom->{$pos};
					$tr_biallelic{$lastscaff}{$pos-$scaffstart}=$self->{_Biallelic}{$chrom}{$pos};
					}
				if ($dir ne '+' && $dir ne '-') {print "Warning: no information about direction for scaffold $chrom: $lastscaff, assuming forward direction!\n"}
				}
			}
		}

	my $tr_loci=Frequencies->new(\%tr_frequencies, \%tr_biallelic, $self->{_Header}, $self->{_Samples});
	return($tr_loci);
	}


1;

