package Frequencies;

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


sub setPops{
	my ($self, $mafsfolder, $mafslist, $npops, $ref_subpops)=@_;

	my ($files, $samples)=Misc::openvcflist($mafsfolder, $mafslist);
	unless (@{$files}==$npops){ die "Wrong number of mafs files specified (", scalar(@{$files}), " vs. $npops)!\n" }
	my ($subfiles, $subsamples);

	for my $i (@$ref_subpops){
		push @{ $subfiles }, $files->[$i];
		push @{ $subsamples }, $samples->[$i];
		}
	$self->{_Files}=$subfiles;
	$self->{_Samples}=$subsamples;
	return;
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


sub getNSamplesDirect{
	my ($self, $chrom, $pos, $pop)=@_;
	if (@_==3){ return($self->{_NSamples}{$chrom}{$pos}) }
	elsif (@_==4){ return($self->{_NSamples}{$chrom}{$pos}[$pop]) }
	else { die "Invalid number of arguments!\n" }
	}

sub setNSamples{
	unless (@_==5){die "Wrong number of arguments!\n"}
	my ($self, $chrom, $pos, $pop, $value)=@_;
	$value=255 if $value>255;
#	my $storagelength=100;	# length of stored character strings. Should be at least 20 to maximize storage efficiency.
	my $storagelength=$self->{_StorageLength};
	my $offset=$pos % $storagelength;
	my $startpos=$pos-$offset;
	unless (defined $self->{_NSamples}{$chrom}{$startpos}[$pop]){
		$self->{_NSamples}{$chrom}{$startpos}[$pop]=(chr(0) x $storagelength);
		}
	substr( $self->{_NSamples}{$chrom}{$startpos}[$pop], $offset, 1, chr($value) );
	return;
	}

sub setNSamplesString{
	unless (@_==5){die "Wrong number of arguments!\n"}
	my ($self, $chrom, $pos, $pop, $value)=@_;
#	my $storagelength=100;	# length of stored character strings. Should be at least 20 to maximize storage efficiency.
	unless ( length($value)==$self->{_StorageLength} ){ die "Invalid input string length!\n" }
	$self->{_NSamples}{$chrom}{$pos}[$pop]=$value;
	return;
	}

sub getNSamples{
	my ($self, $chrom, $pos, $pop)=@_;
#	my $storagelength=100;
	my $offset=$pos % $self->{_StorageLength};
	my $startpos=$pos-$offset;
	if (@_==3){ return($self->{_NSamples}{$chrom}{$startpos}) }
	elsif (@_==4){
		my $value=defined $self->{_NSamples}{$chrom}{$startpos}[$pop] ? ord( substr($self->{_NSamples}{$chrom}{$startpos}[$pop], $offset, 1) ) : 0;
#		unless (defined $value){ warn "Invalid access $chrom:$pos - population: $pop, starting position: $startpos, offset: $offset\n" }
		return($value);
		}
	else { die "Invalid number of arguments!\n" }
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
	my ($self, $mafsfolder, $mafslist, $npops, $ref_loci, $minfreq, $recordbases, $storagelength)=@_;

	$self->{_StorageLength}=$storagelength;
	my ($mafsfiles, $samples);
	if ( defined $self->{_Files} && defined $self->{_Samples} ){
		$mafsfiles=$self->{_Files};
		$samples=$self->{_Samples};
		}
	elsif (defined $mafsfolder && defined $mafslist){
		($mafsfiles, $samples)=Misc::openvcflist($mafsfolder, $mafslist);
		}
	else { die "List of maf files have not been provided or set before calling readMAFs!\n" }
	unless (@{$mafsfiles}==$npops){ die "Wrong number of maf files specified!\n" }

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
						$self->{_Frequencies}{$scaffold}{$pos}[$pop]=($line[5]<$minfreq) ? 0 : $line[5];
						}
					} continue { $prevscaff=$scaffold }
				close $mafsfile;
				}
			}
		}
	return;
	}


sub readDoC{
	my ($self, $infile, $ref_loci, $npops, $ref_poolsizes, $mincoverage, $storagelength)=@_;
	$self->{_StorageLength}=$storagelength;
	my ($minpoolsize, $maxpoolsize)=Misc::minmax($ref_poolsizes);
	print "Maximum of poolsize = $maxpoolsize\n";
	die "Pool sizes larger than 255 individuals currently not supported!\n" if ($maxpoolsize>255);

	CHROM: for my $chrom ($ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# DoC files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			open my $docfile, "-|", "tabix $infile $chrom:$start-$end" or die "Could not open DoC file $infile!\n$!\n";
			print "Aquiring coverage for $chrom:$start-$end.\n";
			my $startpos;
			LINE: while (my $line = <$docfile>){
				my @fields=split(/\s/, $line);
				my $pos=$fields[1]-1;	# DoC files are 1-based. All internal arrays and hashes use a 0-based indexing.
#				my $totcov=($fields[2] <= 255) ? $fields[2] : 255;
				my @coverage=map { ($fields[$_+4] < $mincoverage) ? 0 : $ref_poolsizes->[$_] } (0..$npops-1);
				my $offset=$pos % $storagelength;
				if (!defined $startpos || $pos-$startpos >= $storagelength){
					$startpos=$pos-$offset;
					$self->{_NSamples}{$chrom}{$startpos}=[ map { (chr(0) x $storagelength) } (1..$npops) ];
					}
				substr( $self->{_NSamples}{$chrom}{$startpos}[$_], $offset, 1, chr($coverage[$_]) ) foreach (0..$npops-1);
				}
			close $docfile;
			}
		}
	return;
	}

sub readPoolVcf{
	my ($self, $infile, $ref_loci, $npops, $ref_poolsizes)=@_;
	die "Wrong size of poolsizes array!\n" unless ($npops==@$ref_poolsizes);

	CHROM: for my $chrom ($ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# DoC files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			open my $vcffile, "-|", "tabix $infile $chrom:$start-$end" or die "Could not open DoC file $infile!\n$!\n";
			print "Aquiring frequencies for $chrom:$start-$end.\n";
			LINE: while (my $line = <$vcffile>){
				my @fields=split(/\s/, $line);
				my $pos=$fields[1]-1;	# vcf files are 1-based. All internal arrays and hashes use a 0-based indexing.
				if ($npops+8>$#fields){ warn "Invalid number of populations for $chrom:$pos!\n"; $self->setNSamples($chrom, $pos, $_, 0) foreach (0..$npops-1); next LINE }
				POP: for my $pop (0..$npops-1){
					my $col=$pop+9;
					my @call=split(/:/, $fields[$col]);
					unless (@call>1){ warn "$chrom:$pos - population $pop has invalid genotype field: $fields[$col]\n"; $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					if ($call[0] eq '.'){ $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					my @altalleles=split(/,/, $call[0]);
					my $mac=0;
					my $alt=0;
					for my $ac (@altalleles){ ++$alt if ($ac); $mac+=$ac }
					if ($mac>2*$ref_poolsizes->[$pop]){
						warn "$chrom:$pos - population $pop has an alternative allele count of $mac, while poolsize is only ", 2*$ref_poolsizes->[$pop], "\n"; $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					if ($alt>1){ $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					$self->{_Frequencies}{$chrom}{$pos}[$pop]=$mac/(2*$ref_poolsizes->[$pop]);
					}
				}
			close $vcffile;
			}
		}
	return;
	}

sub readPoolVcfAF{
	my ($self, $infile, $ref_loci, $npops, $ref_poolsizes)=@_;
	die "Wrong size of poolsizes array!\n" unless ($npops==@$ref_poolsizes);

	CHROM: for my $chrom ($ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# DoC files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			open my $vcffile, "-|", "tabix $infile $chrom:$start-$end" or die "Could not open DoC file $infile!\n$!\n";
			print "Aquiring frequencies for $chrom:$start-$end.\n";
			LINE: while (my $line = <$vcffile>){
				my @fields=split(/\s/, $line);
				my $pos=$fields[1]-1;	# vcf files are 1-based. All internal arrays and hashes use a 0-based indexing.
				if ($npops+8 != $#fields){ warn "Invalid number of populations for $chrom:$pos!\n"; $self->setNSamples($chrom, $pos, $_, 0) foreach (0..$npops-1); next LINE }
				if (split(/,/, $fields[4]) > 1){ $self->setNSamples($chrom, $pos, $_, 0) foreach (0..$npops-1); next LINE }
				POP: for my $pop (0..$npops-1){
					my $col=$pop+9;
					my @call=split(/:/, $fields[$col]);
					unless (@call>1){ warn "$chrom:$pos - population $pop has invalid genotype field: $fields[$col]\n"; $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					my $AF=$call[1];
					if ($AF =~ /,/){ warn "$chrom:$pos has multiple allele frequencies!\n"; $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					if ($AF < 0 || $AF > 1.00){ warn "$chrom:$pos has invalid allele frequency: $AF!\n"; $self->setNSamples($chrom, $pos, $pop, 0); next POP }
					$self->{_Frequencies}{$chrom}{$pos}[$pop]=$AF;
					}
				}
			close $vcffile;
			}
		}
	return;
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

