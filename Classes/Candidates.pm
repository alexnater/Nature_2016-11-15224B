package Candidates;

use strict;
use warnings;
#use functions;


sub new{
	my ($class, $candidates)=@_;
	my $self={};
	$self->{_Cand}=defined $candidates ? $candidates : {};
	bless $self, $class;
	return $self;
	}


sub getValue{
	my $self=shift;
	my ($locus, $key)=@_;

	if (defined $key){return($self->{_Cand}[$locus]{$key})}
	elsif (defined $locus){return($self->{_Cand}[$locus])}
	else {return($self->{_Cand})}
	}


sub setValue{
	my $self=shift;
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($locus, $key, $value)=@_;

	$self->{_Cand}[$locus]{$key}=$value;
	return($self->{_Cand}[$locus]{$key});
	}


sub allKeys{
	my $self=shift;
	my $locus=shift;

	if (defined $locus){return(keys %{$self->{_Loci}[$locus]})}
	else {return(0..$#{$self->{_Cand}})}
	}


sub allValues{
	my $self=shift;
	my $locus=shift;

	if (defined $locus){return(values %{$self->{_Loci}[$locus]})}
	else {return(@{$self->{_Cand}})}
	}


sub translateScafftoChrom{
	my ($self, $chrommap, $include_NA)=@_;

	my @translated;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $scaffold=$candidate->{'id'};
		my $start=$candidate->{'start'};
		my $end=$candidate->{'end'};

		CHROM: for my $chrom ($chrommap->allKeys()){
			if (defined $chrommap->getValue($chrom, $scaffold)){
				my $dir=$chrommap->getValue($chrom, $scaffold, 'dir');
				if ($dir eq '+'){
					my $offset=$chrommap->getValue($chrom, $scaffold, 'start');
					push @translated, {'id'=>$chrom, 'start'=>$start+$offset, 'end'=>$end+$offset };
					}
				elsif ($dir eq '-'){
					my $offset=$chrommap->getValue($chrom, $scaffold, 'end');
					push @translated, {'id'=>$chrom, 'start'=>$offset-$end, 'end'=>$offset-$start };
					}
				else {
					print "Warning: no information about direction for scaffold $chrom: $scaffold, assuming forward direction!\n";
					my $offset=$chrommap->getValue($chrom, $scaffold, 'start');
					push @translated, {'id'=>$chrom, 'start'=>$start+$offset, 'end'=>$end+$offset };
					}
				print "Scaffold $scaffold: $start-$end translated to chromosome $chrom: $translated[-1]{'start'}-$translated[-1]{'end'}.\n";
				next CAND;
				}
			}
		print "Warning: scaffold $scaffold has not been located on a chromosome!\n";
		if ($include_NA){push @translated, {'id'=>"NA", 'start'=>0, 'end'=>0 } }
		}

	my $tr_candidates=Candidates->new(\@translated);
	return($tr_candidates);
	}


sub translateChromtoScaff{
	my ($self, $chrommap, $include_NA)=@_;

	my @translated;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'};
		my $end=$candidate->{'end'};

		SCAFF: for my $scaff ($chrommap->allKeys($chrom)){
			if ($start>=$chrommap->getValue($chrom, $scaff, 'start') && $end<=$chrommap->getValue($chrom, $scaff, 'end')){
				my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
				if ($dir eq '+'){
					my $offset=$chrommap->getValue($chrom, $scaff, 'start');
					push @translated, {'id'=>$scaff, 'start'=>$start-$offset, 'end'=>$end-$offset };
					}
				elsif ($dir eq '-'){
					my $offset=$chrommap->getValue($chrom, $scaff, 'end');
					push @translated, {'id'=>$scaff, 'start'=>$offset-$end, 'end'=>$offset-$start };
					}
				else {
					print "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n";
					my $offset=$chrommap->getValue($chrom, $scaff, 'start');
					push @translated, {'id'=>$scaff, 'start'=>$start-$offset, 'end'=>$end-$offset };
					}
				print "Chromosome $chrom: $start-$end translated to scaffold $scaff: $translated[-1][1]-$translated[-1][2].\n";
				next CAND;
				}
			}
		print "Warning: Chromosome $chrom: $start-$end has not been located on a scaffold!\n";
		if ($include_NA){push @translated, {'id'=>"NA", 'start'=>0, 'end'=>0 } }
		}

	my $tr_candidates=Candidates->new(\@translated);
	return($tr_candidates);
	}


sub excludeRegions{
	my ($self, $regions, $min_dist)=@_;

	my @reduced;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'};
		my $end=$candidate->{'end'};

		if (defined $regions->getValue($chrom)){
			for my $reg ($regions->allValues($chrom)){
				if (($reg->{'end'}+$min_dist)>=$start && ($reg->{'start'}-$min_dist)<=$end){
					print "Overlap detected between $chrom $reg->{'start'}-$reg->{'end'} and $start-$end\n";
					next CAND;
					}
				}
			}
		push @reduced, $candidate;
		}

	my $newcand=Candidates->new(\@reduced);
	return($newcand);
	}


sub testbyRegions{
	my ($self, $regions, $min_dist)=@_;

	my @accepted;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'};
		my $end=$candidate->{'end'};

		if (defined $regions->getValue($chrom)){
			for my $reg ($regions->allValues($chrom)){
				if (($reg->{'end'}+$min_dist)>=$start && ($reg->{'start'}-$min_dist)<=$end){
					print "Overlap detected between $chrom $reg->{'start'}-$reg->{'end'} and $start-$end\n";
					push @accepted, 0;
					next CAND;
					}
				}
			}
		push @accepted, 1;
		}

	return(@accepted);
	}


sub testBiallelicMasterlist{
	my ($self, $masterlist, $maxnotbiallelic)=@_;

	my @accepted;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'}+1;
		my $end=$candidate->{'end'}+1;
		my $isbiallelic=1;

		open my $vcffile, "-|", "tabix $masterlist $chrom:$start-$end" or die "Could not open vcf file $masterlist!\n$!\n";
		print "Checking $chrom:$start-$end for non-biallelic variants...\n";
		LINE: while (<$vcffile>){
			if(/^\w+\s/){
				my @line=split("\t", $_);
				if ($line[4] eq '.'){ next LINE }
				my @temp=split(',', $line[4]);
				if (@temp>1){
					print "Position $line[1] not biallelic!\n";
					if (!$maxnotbiallelic--){ $isbiallelic=0; last LINE }
					}
				}
			}
		close $vcffile;
		push @accepted, $isbiallelic;
		}

	return(@accepted);
	}


sub testCoverage{
	my ($self, $vcffiles, $samples, $minind, $minvq, $minmq, $mincov, $minpropcov)=@_;

	my @accepted;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'}+1;
		my $end=$candidate->{'end'}+1;
		my $locuslength=$end-$start+1;
		print "Checking $chrom:$start-$end for sufficient coverage...\n";

		for my $pop (0..$#{$vcffiles}){
			open my $vcffile, "-|", "tabix $vcffiles->[$pop] $chrom:$start-$end" or die "Could not open vcf file $vcffiles->[$pop]!\n$!\n";
			my @coveredind=(0) x scalar( @{$samples->[$pop]} );
			my ($scaffold, $pos, $vq, $mq);
			my $refseq=$candidate->{'refseq'};
			my $prevpos=$start;
			LINE: while (<$vcffile>){
				if(/^\w+\s/){
					my @line=split("\t", $_);
					$scaffold=$line[0];
					$pos=$line[1];
					substr($refseq, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
					$prevpos=$pos;
					if ( uc( substr($refseq, 0, 1) ) eq 'N' ){ next LINE }	# discard if base in reference sequence is hardmasked.
					$vq=$line[5];
					unless ($vq ne '.' && $vq>=$minvq){ next LINE }
					($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
					unless (defined $mq && $mq>=$minmq){ next LINE }
					my $dp=$line[4] eq '.' ? 1 : 2;
					my $index=0;
					for my $ind (@{$samples->[$pop]}){
						my @data=split(":", $line[$ind+8]);
						if (@data>1 && $data[$dp]>=$mincov){ ++$coveredind[$index] }
						++$index;
						}
					}
				}
			close $vcffile;
			my $validind=0;
			for my $indcov (@coveredind){
				my $propcov=($indcov/$locuslength);
#				print "$pop: $indcov - $propcov\n";
				if ($propcov>=$minpropcov){++$validind}
				}
#			print "Population $pop - $chrom:$start-$end: $validind\n";
			unless ($validind>=$minind->[$pop]){push @accepted, 0; print "$chrom:$start-$end has not enough individuals sufficiently covered for population $pop!\n"; next CAND}
			}
		print "$chrom:$start-$end has enough individuals sufficiently covered.\n";
		push @accepted, 1;
		}
	return(@accepted);
	}


sub testBiallelic{
	my ($self, $vcffiles, $samples, $minvq, $mincov, $maxnotbiallelic)=@_;

	my @accepted;

	CAND: for my $candidate (@{$self->{_Cand}}){
		my $chrom=$candidate->{'id'};
		my $start=$candidate->{'start'}+1;
		my $end=$candidate->{'end'}+1;
		my $locuslength=$end-$start+1;
		my %alleles;
		my %genocount;
		my $nbacounter=$maxnotbiallelic;
		print "Checking $chrom:$start-$end for non-biallelic variants...\n";

		for my $pop (0..$#{$vcffiles}){
			open my $vcffile, "-|", "tabix $vcffiles->[$pop] $chrom:$start-$end" or die "Could not open vcf file $vcffiles->[$pop]!\n$!\n";
			my ($scaffold, $pos, $vq);
			LINE: while (<$vcffile>){
				if(/^\w+\s/){
					my @line=split("\t", $_);
					$scaffold=$line[0];
					$pos=$line[1]-1; # The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
					$vq=$line[5];
					unless ($vq ne '.' && $vq>=$minvq){next LINE}
					$genocount{$pos}=[0, 0];
					for my $ind (@{$samples->[$pop]}){
						my @data=split(":", $line[$ind+8]);
						if (@data>1 && $data[1]>=$mincov){	# this is vcf format specific!
							if ($data[0]=~/([0-3])\/([0-3])/){
								++$alleles{$pos}{$1};
								++$alleles{$pos}{$2};
								++$genocount{$pos}[0] if $1 eq $2;
								++$genocount{$pos}[1] if $1 ne $2;
								}
							}
						}
					}
				}
			close $vcffile;
			}
		for my $pos (keys %alleles){
			my $totcount=$genocount{$pos}[0]+$genocount{$pos}[1];
			if ($totcount>5 && !$genocount{$pos}[0]){
				push @accepted, 0;
				print "$chrom:$pos has only heterozygous genotypes , locus excluded! Number of covered individuals: ", $totcount, "\n";
				next CAND; 
				}
			my $nalleles=keys %{ $alleles{$pos} };
			if ($nalleles>2){
				print "$chrom:$pos is non-biallelic with $nalleles alleles!\n";
				for my $allele ( keys %{ $alleles{$pos} } ){
					if ($alleles{$pos}{$allele}==1){ print "$chrom:$pos is non-biallelic with a singleton for allele $allele.\n" }
					}
				if (!$nbacounter--){push @accepted, 0; print "$chrom:$start-$end has more than $maxnotbiallelic non-biallelic sites!\n"; next CAND}
				}
			}
		print "$chrom:$start-$end has ", $maxnotbiallelic-$nbacounter, " non-biallelic loci.\n";
		push @accepted, 1;
		}
	return(@accepted);
	}



1;

