package Genotypes;

use strict;
use warnings;
use Classes::Misc;
use File::Basename;


sub new{
	my ($class, $vcffolder, $vcflist, $simulated, $ref_genotypes, $ref_alleles, $ref_samples, $ref_files)=@_;
	my $self={};

	if (defined $vcffolder && defined $vcflist){
		($ref_files, $ref_samples)=Misc::openvcflist($vcffolder, $vcflist);
		}
	else { warn "No list and folder of vcf files provided!\n" }

	$self->{_Genotypes}=defined $ref_genotypes ? $ref_genotypes : {};
	$self->{_Alleles}=defined $ref_alleles ? $ref_alleles : {};
	$self->{_Files}=$ref_files;
	$self->{_Samples}=$ref_samples;
	$self->{_Simulated}=$simulated ? 1 : 0;
	bless $self, $class;
	return $self;
	}


sub setPops{
	my ($self, $ref_subpops, $ref_subsamples)=@_;

	my $files=$self->{_Files};
	my $samples=$self->{_Samples};
	unless (defined $files && defined $samples){ die "List and folder of vcf files have not been defined!\n" }
	$samples=$ref_subsamples if (defined $ref_subsamples);
	my ($subfiles, $subsamples);

	for my $i (@$ref_subpops){
		push @{ $subfiles }, $files->[$i];
		push @{ $subsamples }, $samples->[$i];
		print STDERR "$i - $files->[$i]: ", join( ',', @{ $samples->[$i] } ), "\n"; 
		}
	$self->{_Files}=$subfiles;
	$self->{_Samples}=$subsamples;
	return;
	}


sub getValue{
	my ($self, $chrom, $pos, $pop, $index)=@_;
	if (@_==3){ return($self->{_Genotypes}{$chrom}{$pos}) }
	elsif (@_==4){ return($self->{_Genotypes}{$chrom}{$pos}[$pop]) }
	elsif (@_==5){
#		warn "Invalid access to $chrom:$pos, $pop - $index!\n" unless (defined $self->{_Genotypes}{$chrom}{$pos}[$pop]);
		return( substr($self->{_Genotypes}{$chrom}{$pos}[$pop], $index, 1) );
		}
	else { die "Invalid number of arguments!\n" }
	}

sub setValue{
	unless (@_==6){die "Wrong number of arguments!\n"}
	my ($self, $chrom, $pos, $pop, $index, $value)=@_;
	substr($self->{_Genotypes}{$chrom}{$pos}[$pop], $index, 1, $value);
	return($value);
	}


sub allKeys{
	unless (@_==2){die "Wrong number of arguments!\n"}
	my ($self, $id)=@_;
	return( keys %{$self->{_Genotypes}{$id}} )
	}


sub getAncestral{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	return($self->{_Ancestral}{$id}{$pos});
	}

sub setAncestral{
	unless (@_==4){ die "Wrong number of arguments!\n" }
	my ($self, $id, $pos, $value)=@_;
	$self->{_Ancestral}{$id}{$pos}=$value;
	return;
	}


sub getBase{
	my ($self, $chrom, $pos, $i, $j)=@_;
	if (@_==3){ return( split('', $self->{_Bases}{$chrom}{$pos}) ) }	# returns all bases present at site as list.
	elsif (@_==4){
		if ($i eq '.'){ return("N") } 
		elsif ($i>=0 && $i<4){ return( substr($self->{_Bases}{$chrom}{$pos}, $i, 1) ) }
		else { return }
		}
	elsif (@_==5){
		my $gt=substr($self->{_Genotypes}{$chrom}{$pos}[$i], $j, 1);
		unless (defined $gt){ return }
		if ($gt eq '.'){ return("N") }
		else { return( substr($self->{_Bases}{$chrom}{$pos}, $gt, 1) ) }
		}
	else {die "Invalid number of arguments!\n"}
	}

sub getAlleleforBase{
	unless (@_==3 || @_==4){ die "Wrong number of arguments!\n" }
	my ($self, $chrom, $pos, $base)=@_;
	my @baselist=split('', $self->{_Bases}{$chrom}{$pos});
	my %basehash=map { $baselist[$_] => $_ } (0..$#baselist);
	if (defined $base){ return( $basehash{$base} ) }
	else { return(\%basehash) }		
	}


sub isAllele{
	unless (@_==4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, $allele)=@_;
	unless ($allele>=0 && $allele<=3){die "Invalid allele specifier!\n"}
	if ($self->{_Alleles}{$id}{$pos} & (1 << $allele)){return 1}	# shift first bit to the left to position specified by $allele and test if bit is set.
	else {return 0}
	}

sub getAlleles{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	my @alleles=map { $self->{_Alleles}{$id}{$pos} & (1 << $_) ? $_ : () } (0..3);
	return (@alleles);
	}

sub getPolarizedAlleles{	# returns a list where the first element is the ancestral allele, even if not physically present in the set.
	unless (@_==3 || @_==4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, $ancestral)=@_;
	$ancestral=$self->{_Ancestral}{$id}{$pos} unless (defined $ancestral);
	my $bits=$self->{_Alleles}{$id}{$pos};	# make copy of bitstring.
	if (defined $ancestral && $ancestral>=0 && $ancestral<4){
		$bits &= ~(1 << $ancestral);	# unset bit for ancestral allele.
		}
	else { $ancestral=undef }
	my @alleles=map { $bits & (1 << $_) ? $_ : () } (0..3);
	return ($ancestral, @alleles);
	}

sub setAlleles{
	unless (@_>=4){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos, @alleles)=@_;
	for my $allele (@alleles){
		unless ($allele>=0 && $allele<=3){die "Invalid allele specifier!\n"}
		$self->{_Alleles}{$id}{$pos} |= (1 << $allele);	# shift first bit to the left to position specified by $allele and set bit to 1.
		}
	return;
	}

sub NAlleles{
	unless (@_==3){die "Wrong number of arguments!\n"}
	my ($self, $id, $pos)=@_;
	my $count1=unpack( '%32b*', chr($self->{_Alleles}{$id}{$pos}) );
	my $count2=0;
	for my $allele (0..3){ if ($self->{_Alleles}{$id}{$pos} & (1 << $allele)){ ++$count2 } }
	die "Non-matching number of alleles found for $id:$pos, $count1 vs. $count2! Bitvalue: $self->{_Alleles}{$id}{$pos}\n" unless ($count1==$count2); 
	return($count1);
	}


sub getCoverage{
	my ($self, $id, $pos, $pop, $ind)=@_;
	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==3){ return($self->{_Coverage}{$id}{$pos}) }
	elsif (@_==4){ return($self->{_Coverage}{$id}{$pos}[$pop]) }
	elsif (@_==5){
		my $substring=substr($self->{_Coverage}{$id}{$pos}[$pop], $ind, 1);
		unless (defined $substring){ warn "Invalid access: $id, $pos, $pop, $ind!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub validatevcfList{
	my ($self, $ref_loci, $vcffolder, $vcflist, $plcol)=@_;
	$plcol=4 unless (defined $plcol);
	my @samples;
	my @line;
	my $pos;
	my $iterator;

	my ($vcffiles, $samples)=Misc::openvcflist($vcffolder, $vcflist);

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		my %alleles;
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			for my $pop (0..$#{$vcffiles}){
				open my $vcffile, "-|", "tabix $vcffiles->[$pop] $chrom:$start-$end" or die "Could not open vcf file $vcffiles->[$pop]!\n$!\n";
				print "Validating genotypes for $chrom:$start-$end, population $pop.\n";
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						@line=split(/\s/, $_);
						next LINE if ($line[4] eq '.');
						my $index=0;
						for my $ind (0..$#{$samples->[$pop]}){
							my $col=$samples->[$pop][$ind]+8;
							if ($col>$#line){die "Invalid sample index specified for population $pop!\n"}
							if ($line[$col]=~/([0-1])\/([0-1])/){
								my @pl=split(',', ( split(':', $line[$col]) )[$plcol]);
								if ($1==0 && $2==0){ warn "$chrom:$line[0] $pop - $ind has wrong genotype! $1/$2 vs. ", join(',', @pl), "\n" unless ($pl[0]==0) }
								elsif ($1==1 && $2==1){ warn "$chrom:$line[0] $pop - $ind has wrong genotype! $1/$2 vs. ", join(',', @pl), "\n" unless ($pl[2]==0) }
								else { warn "$chrom:$line[0] $pop - $ind has wrong genotype! $1/$2 vs. ", join(',', @pl), "\n" unless ($pl[1]==0) }
								}
							}
						}
					}
				close $vcffile;
				}
			}
		}
	return;
	}


sub readvcfList{
	my ($self, $ref_loci, $mingtq, $suballeles, $recordbases, $recordcoverage, $onlyphased, $randomize, $verbose)=@_;

	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of vcf files have not been defined!\n" }
	my $ref_samples=$self->{_Samples};
	my (@line, $pos, $iterator, @reqphasing);
	if (defined $onlyphased && ref($onlyphased) eq 'ARRAY'){
		die "Wrong number of elements in onlyphased array!\n" unless (@$onlyphased==@{ $self->{_Files} });
		@reqphasing=@$onlyphased;
		}
	elsif (defined $onlyphased){ @reqphasing=map { $onlyphased } @{ $self->{_Files} } }
	else { @reqphasing=(0) x @{ $self->{_Files} } }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			if ($suballeles){ $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus) }
			for my $pop (0..$#{ $self->{_Files} }){
				my $phased=$reqphasing[$pop];
				open my $vcffile, "-|", "tabix $self->{_Files}[$pop] $chrom:$start-$end" or die "Could not open vcf file $self->{_Files}[$pop]!\n$!\n";
				print "Aquiring genotypes for $chrom:$start-$end, population $pop.\n" if $verbose;
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						@line=split(/\s/, $_);
						$pos=$line[1]-1;	# The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
						next LINE if ($self->{_Genotypes}{$line[0]}{$pos}[$pop]);	# skip segregating position if already read in before.
						if ($recordbases){
							my $basestring="$line[3]" . join('', split(',' , $line[4]) );
							my $current=$self->{_Bases}{$line[0]}{$pos};
							if (!defined $current){ $self->{_Bases}{$line[0]}{$pos}=$basestring }
							elsif (defined $current && $current ne $basestring){
								warn "$chrom:$pos, population $pop, bases $basestring don't match with bases from other populations $current!\n";
								my ($newbasestring, $ref_translation)=Misc::harmonizeBases($current, $basestring);
								$self->{_Bases}{$line[0]}{$pos}=$newbasestring;
								warn "Merged the two basestrings to $newbasestring.\n";
								warn "Applying the following translation of alleles to population $pop:\n";
								print STDERR "$_ to $ref_translation->{$_}\n" foreach (keys %$ref_translation);
#								print STDERR "Old line:\n", join("\t", @line), "\n";
								$line[4]=join(',', split('', substr($newbasestring, 1) ) );
								for my $call (@line[9..$#line]){
									$call=~s/([0-3])(\/|\|)([0-3])/$ref_translation->{$1}$2$ref_translation->{$3}/;
									}
#								print STDERR "New line:\n", join("\t", @line), "\n";
								}
							}
#							$self->{_Bases}{$line[0]}{$pos}="$line[3]" . join('', split(',' , $line[4]) ) }
						if ($recordcoverage){ $self->{_Coverage}{$line[0]}{$pos}[$pop]=( chr(0) x scalar(@{ $ref_samples->[$pop] }) ) }
						$self->{_Alleles}{$line[0]}{$pos}=0 unless defined $self->{_Alleles}{$line[0]}{$pos};
						my $index=0;
						for my $ind (0..$#{ $ref_samples->[$pop] }){
							my $col=$ref_samples->[$pop][$ind]+8;
							if ($col>$#line){ die "Invalid sample index specified for population $pop!\n" }
							if ($line[$col]=~/([0-3])(\/|\|)([0-3])/){
								my @call=split(':', $line[$col]);
								my $infofields=(@call==5) ? 1 : 0;	# in case sample fields only contain genotype call.
								if ($phased && $2 ne '|' && $1 ne $3){
									if ($randomize){ $self->{_Genotypes}{$line[0]}{$pos}[$pop].=( int(rand(2)) ) ? $1 . $3 : $3 . $1 }
									else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
									}
								elsif ($infofields && $call[-2]<$mingtq){ $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
								else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=$1 . $3 }
								if ($recordcoverage && $infofields){
									my $coverage=$call[-3]<255 ? $call[-3] : 255;
									substr($self->{_Coverage}{$line[0]}{$pos}[$pop], $ind, 1)=chr($coverage);
									}
								if ($suballeles && defined $iterator->[$pop][$index]){	# only records allele states for subsampled individuals.
									if (2*$ind==$iterator->[$pop][$index]){ $self->setAlleles($line[0], $pos, $1); ++$index }
									if (2*$ind+1==$iterator->[$pop][$index]){ $self->setAlleles($line[0], $pos, $3); ++$index }
									}
								else { $self->setAlleles($line[0], $pos, $1, $3) }
								}
							else { $self->{_Genotypes}{$line[0]}{$pos}[$pop].=".." }
							}
						}
					}
				close $vcffile;
				}
			}

		SITE: for my $seg ( keys %{ $self->{_Genotypes}{$chrom} } ){	# check if genotype vcf files have synchronized polymorphic positions.
			my $stringref=$self->getValue($chrom, $seg);
			for my $pop (0..$#{ $self->{_Files} }){
				my $nalleles=2*scalar( @{ $ref_samples->[$pop] } );
				if (defined $stringref->[$pop]){
					unless (length($stringref->[$pop])==$nalleles){
						warn "Genotype string for $chrom:$seg - population $pop has wrong length (", length($stringref->[$pop]), " vs. $nalleles)!",
							"Setting all genotypes in population to missing.\n";
						$stringref->[$pop]='.' x $nalleles;
						}
					}
				else {
#					warn "Genotype string for $chrom:$seg - population $pop is missing! Setting all genotypes to reference allele.\n";
					$stringref->[$pop]='0' x $nalleles;
					$self->setAlleles($chrom, $seg, '0');	# make sure that reference allele is set as present for the site.
					if ($recordcoverage){ $self->{_Coverage}{$chrom}{$seg}[$pop]=chr(0) x $nalleles }	# the use of the record coverage function in combination with non-synchronyzed vcf files is dangerous!

#					if ($reqphasing[$pop] && !$randomize){	# this is questionable, as it depends on if there is a genotype-specific filter. the problem is that we don't now if the position is covered in a population without an entry in the SNP vcf. phasing wouldn't be an issue if we are confident that all samples homozygot reference genotypes. prime example why unsynchronized vcf files are a bad idea!
#						$stringref->[$pop]='.' x $nalleles;
#						}
#					else {
#						$stringref->[$pop]='0' x $nalleles;
#						$self->setAlleles($chrom, $seg, '0');	# make sure that reference allele is set as present for the site.
#						if ($recordcoverage){ $self->{_Coverage}{$chrom}{$seg}[$pop]=chr(0) x $nalleles }	# the use of the record coverage function in combination with non-synchronyzed vcf files is dangerous!
#						}
					}
				}
			}

		}
	return;
	}


sub readOutgroups{
	my ($self, $ref_loci, $useinner)=@_;
	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of vcf files have not been defined!\n" }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		my %alleles;
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			for my $pop (0..$#{ $self->{_Files} }){
				open my $vcffile, "-|", "tabix $self->{_Files}[$pop] $chrom:$start-$end" or die "Could not open vcf file $self->{_Files}[$pop]!\n$!\n";
				print "Aquiring genotypes for $chrom:$start-$end, population $pop.\n";
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						my @call=split(/\s/, $_);
						my $pos=$call[1]-1;	# The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
#						print "Acquiring outgroup data for position $pos\n";
						$alleles{$pos}[$pop]=0 unless defined $alleles{$pos}[$pop];
						for my $ind (@{ $self->{_Samples}[$pop] }){
							my $col=$ind+8;
							if ($col>$#call){die "Invalid sample index specified for population $pop!\n"}
							if ($call[$col]=~/([0-3])\/([0-3])/){
								$alleles{$pos}[$pop] |= (1 << $1);
								$alleles{$pos}[$pop] |= (1 << $2);
								}
							}
						}
					}
				close $vcffile;
				}
			}
		POS: for my $pos (keys %alleles){
#			print "$chrom:$pos; $alleles{$pos}[0], $alleles{$pos}[1], $alleles{$pos}[2], $alleles{$pos}[3]\n";
			unless ($alleles{$pos}[0] || $alleles{$pos}[1]){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# Both outgroups are missing.
			my $missing=0;
			unless ($alleles{$pos}[0] && $alleles{$pos}[1]){ $missing=1 }
			my $all=$alleles{$pos}[0] | $alleles{$pos}[1];
			my $shared=$alleles{$pos}[0] & $alleles{$pos}[1];
#			print "$chrom:$pos; all: $all, shared: $shared, missing: $missing\n"; 
			my $outer;
			my ($nalleles, $outstate)=Misc::countbits2($all);
#			print "$chrom:$pos; nalleles: $nalleles, outstate: $outstate\n";
			if (!$missing && $nalleles==1){ $self->{_Ancestral}{$chrom}{$pos}=$outstate; next POS }
			elsif (!$missing && $nalleles==2){
				if ($alleles{$pos}[0]==$alleles{$pos}[1]){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# both outgroups polymorphic.
				elsif (!$shared){ $outer=$all }	# both outgroups fixed for different alleles.
				else { $outer=$shared }	# 1 outgroup fixed, one polymorphic.
				}
			elsif ($missing && $nalleles==1){ $outer=$all }	# 1 outgroup missing, 1 fixed.
			elsif ($missing && $nalleles>1){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# 1 outgroup missing, one polymorphic.
			else { $self->{_Ancestral}{$chrom}{$pos}=4; next POS }

			unless ($useinner){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }
			unless ($alleles{$pos}[2] & $alleles{$pos}[3]){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# One of the two inner outgroups is missing, or they are fixed for different alleles.
			my $inner=$alleles{$pos}[2] | $alleles{$pos}[3];
			($nalleles, $outstate)=Misc::countbits2($inner);
#			print "$chrom:$pos; inner: $inner\n"; 
#			print "$chrom:$pos; nalleles: $nalleles, outstate: $outstate\n";
			unless ($nalleles==1){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# One or both of the two inner outgroups is polymorphic.

			my $tempoverall=$outer & $inner;
			unless ($tempoverall){ $self->{_Ancestral}{$chrom}{$pos}=4; next POS }	# no overlap between inner and outer outgroup alleles.
			($nalleles, $outstate)=Misc::countbits2($tempoverall);
			if ($nalleles==1){ $self->{_Ancestral}{$chrom}{$pos}=$outstate }
			else { $self->{_Ancestral}{$chrom}{$pos}=4 }
			} # continue {print "Outgroup state: $self->{_Ancestral}{$chrom}{$pos}\n"}
		}
	return;
	}


sub readArlequin{
	my ($self, $infile, $chrom, $recordbases, $shuffle)=@_;

	my $pop;
	my ($samplesize, $sample);
	my (@positions, @samplesizes, @haplotypes);
	my $posflag=0;

	open my $input, "<", $infile or die "Could not open Arlequin file $infile!\n$!\n";
	LINE: while (<$input>){
		if (/polymorphic positions/){ $posflag=1 }
		elsif(/^#\d+,\s/ && $posflag){
			my $line=substr($_, 1);
			@positions=split(', ', $line);
			$posflag=0;
			}
		elsif(/SampleName=/){
			if (defined $pop){ ++$pop }
			else { $pop=0 }
			$samplesize=0; $sample=0;
			}
		elsif(/SampleSize=(\d+)/){ $samplesize=$1; push @samplesizes, $samplesize/2 }
		elsif(/^\s?(\d+)_\d+\s\d+\s/){
			if ($1 != $pop+1){ die "Wrong population specified in Arlequin file!\n" }
			unless ($samplesize){ die "Wrong sample size specified in Arlequin file!\n" }
			my @line=split(/\s+/, $_);
			unless ( @positions==length($line[2]) ){ die "Haplotype string has wrong length! ", scalar(@positions), " vs. ", length($line[2]), "\n" }
			push @{$haplotypes[$pop]}, $line[2];
			--$samplesize;
			++$sample;
			}
		}
	close $input;

	my @indices;
	for my $popindex (0..$#haplotypes){
		my @temp;
		if ($shuffle){ @temp=Misc::shuffle( 0..$#{$haplotypes[$popindex]} ) }
		else { @temp=( 0..$#{$haplotypes[$popindex]} ) }
		print STDERR "$popindex: ", join(',', @temp), "\n";
		push @indices, \@temp;
		}

	for my $offset (0..$#positions){
		my $position=$positions[$offset]-1;	# internal storage is 0-based.
		my ($refbase, $altbase);
		for my $popindex (0..$#haplotypes){
			for my $ind ( @{ $indices[$popindex] } ){
				my $state=substr($haplotypes[$popindex][$ind], 0, 1, "");	# chew up haplotype string to speed up substring extraction.
				$refbase=$state unless defined $refbase;	# first base is defined as reference base.
				if ($state eq $refbase){ $self->{_Genotypes}{$chrom}{$position}[$popindex].=0 }
				else {
					if (defined $altbase && $altbase ne $state){ warn "Found a triallelic SNP in the data at $position! Method assumes strictly biallelic data!\n" }
					elsif (!defined $altbase){ $altbase=$state }
					$self->{_Genotypes}{$chrom}{$position}[$popindex].=1;	# Attention! Assumes infinitive site model! 
					}
				}
			}
		if ($recordbases){ $self->{_Bases}{$chrom}{$position}=$refbase . $altbase }
		$self->setAlleles($chrom, $position, 0, 1);
		}
	return(@samplesizes);
	}

sub readMS{
	my ($self, $infile, $chrom_prefix, $ref_chromlength, $ref_samplesizes, $ref_bases, $shuffle, $verbose)=@_;

	my $locus=-1;
	my $locusname;
	my ($pop, $ind);
	my (@positions, @haplotypes);
	my $nsegsites;
	my $monom=0;

	open my $input, "<", $infile or die "Could not open MS file $infile!\n$!\n";
	LINE: while (<$input>){
		chomp;
		if (/^\/\//){
			++$locus;
			@positions=();
			@haplotypes=();
			if (defined($pop) && !$monom){
				unless ($ind==$ref_samplesizes->[$pop] && $pop==$#{ $ref_samplesizes } ){ die "Wrong information about sample sizes and number of populations for locus $locus!\n" }
				}
			$locusname=$chrom_prefix . ($locus + 1);
			$pop=0; $ind=0;
			$monom=0;
			}
		elsif (/segsites: 0$/){ $monom=1; warn "Locus $locus is monomorphic!\n" if $verbose }
		elsif (/segsites:\s(\d+)$/){ $nsegsites=$1 }
		elsif (/positions:\s(.*)/){
			@positions=( split(/\s/, $1) );
			die "Wrong number of segregating sites in list of positions for locus $locus!\n" unless (@positions==$nsegsites);
			}
		elsif(/^\d+$/){
			unless (length($_)==$nsegsites){ die "Haplotype string has wrong length for locus $locus! $nsegsites vs. ", length($_), "\n" }
			push @{ $haplotypes[$pop] }, $_;
			if ( $ind<($ref_samplesizes->[$pop]-1) ){ ++$ind }
			elsif ( $pop<($#{ $ref_samplesizes }) ){ $ind=0; ++$pop }
			elsif ( $ind==($ref_samplesizes->[$pop]-1) && $pop==$#{ $ref_samplesizes }){
				my @indices;
				for my $popindex (0..$#haplotypes){
					my @temp;
					if ($shuffle){ @temp=Misc::shuffle( 0..$#{$haplotypes[$popindex]} ) }
					else { @temp=( 0..$#{$haplotypes[$popindex]} ) }
					push @indices, \@temp;
					}
				$self->{_Genotypes}{$locusname}={};
				for my $offset (0..$#positions){
					my $position=int($positions[$offset]*$ref_chromlength->[$locus]);
					while ( exists $self->{_Genotypes}{$locusname}{$position} ){ warn "Position $position occurs twice in the data set for locus $locusname\n" if $verbose; ++$position }
					if ($position>=$ref_chromlength->[$locus]){ warn "Segsite ", $offset+1, " could not be placed within specified region!\n" }
					for my $popindex (0..$#haplotypes){
						for my $ind ( @{ $indices[$popindex] } ){
							my $state=substr($haplotypes[$popindex][$ind], 0, 1, "");	# chew up haplotype string to speed up substring extraction.
							if ($state==0 || $state==1){ $self->{_Genotypes}{$locusname}{$position}[$popindex].=$state }
							elsif ($state==2){ $self->{_Genotypes}{$locusname}{$position}[$popindex].="." }
							else { warn "Unknown allele state $state found for locus $locusname, at segsite ", $offset+1, "!\n" }
							}
						}
					if (defined $ref_bases){ $self->{_Bases}{$locusname}{$position}=$ref_bases->[0] . $ref_bases->[1] }
					$self->setAlleles($locusname, $position, 0, 1);
					}
				++$ind;
				}
			else { die "Wrong information about sample sizes and number of populations for locus $locus!\n" }
			}
		}
	close $input;
	return($locus+1);
	}


sub extractRegions{
	my ($self, $regions)=@_;

	my %reduced;

	CHROM: for my $chrom (sort $regions->allKeys()){
		if (defined $self->{_Genotypes}{$chrom}){
			REG: for my $reg ($regions->allValues($chrom)){
				POSITION: for my $pos (sort {$a<=>$b} keys %{$self->{_Genotypes}{$chrom}}){
					if ($pos>=$reg->{'start'} && $pos<=$reg->{'end'}){$reduced{$chrom}{$pos}=$self->{_Genotypes}{$chrom}{$pos}}
					elsif ($pos>$reg->{'end'}){next REG}
					}
				}
			}
		else {print "No genotypes available for chromosome $chrom!\n"}
		}

	my $newgenotypes=Genotypes->new(undef, undef, $self->{_Simulated}, \%reduced, undef, $self->{_Samples}, $self->{_Files});
	return($newgenotypes);
	}


sub translateScafftoChrom{	# not up-to-date!
	my ($self, $chrommap)=@_; # bed file format with chromosome locations

	my %tr_genotypes;
	my %tr_biallelic;

	SCAFF: while (my ($scaff, $ref_scaff)=each %{$self->{_Genotypes}}){
		CHROM: for my $chrom ($chrommap->allKeys()){
			if (defined $chrommap->getValue($chrom, $scaff)){
				my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
				if ($dir eq '-'){
					my $offset=$chrommap->getValue($chrom, $scaff, 'end');
					POSITION: while (my ($pos, $ref_pos)=each %{$ref_scaff}){
						$tr_genotypes{$chrom}{$offset-$pos}=$ref_pos;
						$tr_biallelic{$chrom}{$offset-$pos}=$self->{_Biallelic}{$scaff}{$pos};
						}
					}
				else {
					my $offset=$chrommap->getValue($chrom, $scaff, 'start');
					POSITION: while (my ($pos, $ref_pos)=each %{$ref_scaff}){
						$tr_genotypes{$chrom}{$pos+$offset}=$ref_pos;
						$tr_biallelic{$chrom}{$pos+$offset}=$self->{_Biallelic}{$scaff}{$pos};
						}
					}
				if ($dir ne '+' && $dir ne '-') {print "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n"}
				next SCAFF;
				}
			}
		print "Warning: scaffold $scaff has not been located on a chromosome! Loci/Positions on this scaffold are excluded in the output!\n";
		}
	my $tr_loci=Genotypes->new(undef, undef, $self->{_Simulated}, \%tr_genotypes, undef, $self->{_Samples}, $self->{_Files});
	return($tr_loci);
	}


sub translateChromtoScaff{	# not up-to-date!
	my ($self, $chrommap)=@_; # bed file format with chromosome locations

	my %tr_genotypes;
	my %tr_biallelic;

	CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_Genotypes}}){
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
					$tr_genotypes{$lastscaff}{$scaffend-$pos}=$ref_chrom->{$pos};
					$tr_biallelic{$lastscaff}{$scaffend-$pos}=$self->{_Biallelic}{$chrom}{$pos};
					}
				else {
					$tr_genotypes{$lastscaff}{$pos-$scaffstart}=$ref_chrom->{$pos};
					$tr_biallelic{$lastscaff}{$pos-$scaffstart}=$self->{_Biallelic}{$chrom}{$pos};
					}
				if ($dir ne '+' && $dir ne '-') {print "Warning: no information about direction for scaffold $chrom: $lastscaff, assuming forward direction!\n"}
				}
			}
		}

	my $tr_loci=Genotypes->new(undef, undef, $self->{_Simulated}, \%tr_genotypes, undef, $self->{_Samples}, $self->{_Files});
	return($tr_loci);
	}


sub genotypestoMS{
	my $self=shift @_;
	my $outfile=shift @_;
	my $ref_loci=shift @_; # %loci{chromosome}[locus]{start/end/individuals}([population][ind_index])
	my $ref_missing=shift @_; # %missing{chromosome}[locus][allind]
	my $mode=shift @_;
	my $haplotize=shift @_;
	my $maxnotbiallelic=shift @_;
	my $outinfo=shift @_;

	my $npops=$ref_missing->getNPop();
	my @popsizes=map { $ref_missing->getNInd($_) } (0..$npops-1);

	open my $output, ">", $outfile or die "Could not open ms file $outfile!\n";

	my $l=0;
	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $iterator;
			my @haplotypes;
			my $notbiallelic=0;
			my @segsites;
			my @positions;
			my $outgroupstates;
			SITE: for my $seg (sort {$a<=>$b} keys %{$self->{_Genotypes}{$chrom}}){
				if ($seg>=$start && $seg<=$end){
					my $stringref=$self->getValue($chrom, $seg);
					for my $pop (0..$npops-1){	# check if genotype vcf files have synchronized polymorphic positions.
						if ($stringref->[$pop]){
							if(length($stringref->[$pop])!=2*$popsizes[$pop]){
								warn "Genotype string for $chrom:$seg - population $pop has wrong length! Setting all genotypes in population to missing.\n";
								$stringref->[$pop]='.' x (2*$popsizes[$pop]);
								}
							}
						else {
#							warn "Genotype string for $chrom:$seg - population $pop is missing! Setting all genotypes to reference allele.\n";
							$stringref->[$pop]='0' x (2*$popsizes[$pop]);
							$self->setAlleles($chrom, $seg, '0');	# make sure that reference allele is set as present for the site.
							}
						}
					my $nalleles=$self->NAlleles($chrom, $seg);
					if ($nalleles>2){
						++$notbiallelic;
						if ($notbiallelic>$maxnotbiallelic){
							$ref_loci->setExcluded($chrom, $locus, 1);
							print "Locus $chrom - $locus is not biallelic at position $seg. ${notbiallelic}th occurrence, locus excluded!\n";
							next LOCUS;
							}
						else {
							for my $pop (0..$npops-1){ $stringref->[$pop]='.' x length($stringref->[$pop]) }
							print "Locus $chrom - $locus is not biallelic at position $seg. ${notbiallelic}th occurrence, so still included.\n";
							}
						}
					elsif ($nalleles<2){
						print "Position $seg of locus $chrom - $locus is monomorphic! Position excluded.\n";
						next SITE;
						}
					push(@segsites, $seg);
					push @positions, ($seg+1-$start)/$locuslength;	# the first base has position 1/locuslength, the last base 1.
					}
				last SITE if $seg>=$end;
				}

			if ($mode eq 'locus' && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus)){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				if ($haplotize){
					my @hapiterator;
					for my $pop (0..$npops-1){
						for (my $ind=0;$ind<$#{$iterator->[$pop]};$ind=$ind+2){
							push @{$hapiterator[$pop]}, $iterator->[$pop][$ind]+int(rand(2));
							}
						}
					$iterator=\@hapiterator;
					}
				print "Locus $chrom: $locus in locus mode\n";
				}
			elsif ($mode eq 'position' && defined $ref_loci->getreftoIndividuals('position', $chrom, $segsites[0])){ print "Locus $chrom: $locus in position mode\n" }
			else { warn "Information about subsampling mode not available!\n" }

			SITE: for my $seg (@segsites){
				if ($mode eq 'position'){$iterator=$ref_loci->getreftoIndividuals('position', $chrom, $seg)}
				if ($npops!=scalar(@{$iterator})){die "Wrong information about number of populations!\n"}
				my $offset=$seg-$start;
				my $outstate=defined $self->{_Ancestral}{$chrom}{$seg} ? $self->{_Ancestral}{$chrom}{$seg} : 4;
				my $rt=0;

				for my $pop (0..$npops-1){
					my $gtstring=$self->getValue($chrom, $seg, $pop);
#					print "$chrom: $seg - $gtstring\n";
					if ($self->{_Alleles}{$chrom}{$seg}>3 && $self->{_Alleles}{$chrom}{$seg}<7){
						if ($self->{_Alleles}{$chrom}{$seg}==5){
							my $N=($gtstring=~tr/2/1/);
							if ($outstate==2){ $outstate=1 }
							print "Position $seg of locus $chrom - $locus, $N occurences of Allele 2 have been recoded to Allele 1 for subsampled individuals of population $pop.\n";
							}
						elsif ($self->{_Alleles}{$chrom}{$seg}==4 || $self->{_Alleles}{$chrom}{$seg}==6){
							my $N=($gtstring=~tr/2/0/);
							if ($outstate==2){ $outstate=0 }
							print "Position $seg of locus $chrom - $locus, $N occurences of Allele 2 have been recoded to Allele 0 for subsampled individuals of population $pop.\n";
							}
						}
						
					my @newiterator;
					if ($haplotize && $mode eq 'position'){
						for (my $ind=0;$ind<$#{$iterator->[$pop]};$ind=$ind+2){
							push @newiterator, $iterator->[$pop][$ind]+int(rand(2));
							}
						}
					elsif ($mode eq 'locus'){ @newiterator=( @{ $iterator->[$pop] } ) }
					else { @newiterator=( 0..2*$popsizes[$pop]-1 ) }

					for my $ind (0..$#newiterator){
						my $index=$newiterator[$ind];
						my $allind=int($index/2)+$rt;
						my $gt=substr($gtstring, $index, 1);
#						my $gt=$self->getValue($chrom, $seg, $pop, $index);
						my $miss=$ref_missing->getValue($chrom, $start, $allind, $offset);
#						my $miss=$ref_missing->getValue3($chrom, $seg, $allind);
						if (!defined $gt){
							print "Genotypes missing for position $chrom: $seg!\n";
							$ref_loci->setExcluded($chrom, $locus, 1);
							next LOCUS;
							}
						if (!defined $miss){
							print "Missing data missing for position $chrom: $seg!\n";
							$ref_loci->setExcluded($chrom, $locus, 1);
							next LOCUS;
							}
						if ($miss==0){$haplotypes[$pop][$ind].=2}
						elsif ($gt eq '.'){
							$ref_missing->setValue($chrom, $start, $allind, $offset, 0);
#							$ref_missing->setValue3($chrom, $seg, $allind, 0);
							$haplotypes[$pop][$ind].=2;	# Missing genotypes can either be represented as 2 for msABC compatibility or 3 for additional information.
							}
						else {$haplotypes[$pop][$ind].=$gt}
						}
					$rt+=$popsizes[$pop];
					}
				$outgroupstates.="$outstate";
				}
			$ref_loci->setExcluded($chrom, $locus, 0);
			print $output "// Locus: $l - $chrom: $segsites[0]-$segsites[-1]\n";
			print $output "segsites: ", scalar(@positions), "\npositions: ", join("\t", @positions), "\n";
			if ($outinfo){print $output "ancestral states: $outgroupstates\n"}
			for my $pop (0..$#haplotypes){
				for my $haplo (@{$haplotypes[$pop]}){print $output "$haplo\n"}
				}
			print $output "\n";
			++$l;
			}
		}
	close $output;
	return;
	}


sub phaseGenotypes{
	my ($self, $id, $ref_positions, $ref_haplotypes, $inpute)=@_;
	for my $pop (0..$#{ $ref_haplotypes }){
		$self->{_isPhased}{$id}{$_-1}[$pop]=1 foreach (@$ref_positions);
		my $nind=@{ $ref_haplotypes->[$pop] }/2;
		for my $ind (0..$nind-1){
			(my $hapstring1=$ref_haplotypes->[$pop][2*$ind])=~s/\?/\./g;
			(my $hapstring2=$ref_haplotypes->[$pop][2*$ind+1])=~s/\?/\./g;
			unless (length($hapstring1)==@$ref_positions){ 
				die "Unequal length of haplotype string 1 (", length($hapstring1), ") and position array (", scalar(@$ref_positions), ")!\n";
				}
			unless (length($hapstring2)==@$ref_positions){ 
				die "Unequal length of haplotype string 2 (", length($hapstring2), ") and position array (", scalar(@$ref_positions), ")!\n";
				}
			POSITION: for my $pos (@$ref_positions){
				my $newallele1=substr($hapstring1, 0, 1, "");
				my $newallele2=substr($hapstring2, 0, 1, "");
				my $oldalleles=substr($self->{_Genotypes}{$id}{$pos-1}[$pop], 2*$ind, 2);
				my $oldmissing=($oldalleles=~/\./);
				my $newmissing=($newallele1 eq '.' || $newallele2 eq '.');
				substr($self->{_Genotypes}{$id}{$pos-1}[$pop], 2*$ind, 2, $newallele1 . $newallele2) unless ($oldmissing && !$inpute);
				unless ($oldmissing || $newmissing || $oldalleles eq $newallele1 . $newallele2 || $oldalleles eq $newallele2 . $newallele1){
					warn "Error: phased genotype doesn't match previous genotype at $id:$pos ($oldalleles vs. ", $newallele1 . $newallele2, ")!\n";
					}
				}
			}
		}
	return;
	}


sub updateVCFs{
	my ($self, $vcffolder, $vcflist, $outfolder, $outsuffix, $npops, $ref_loci, $impute, $onlyphased, $verbose)=@_;

	my ($vcffiles, $samples);
	if ( defined $self->{_Files} && defined $self->{_Samples} ){
		$vcffiles=$self->{_Files};
		$samples=$self->{_Samples};
		}
	elsif (defined $vcffolder && defined $vcflist){
		($vcffiles, $samples)=Misc::openvcflist($vcffolder, $vcflist);
		}
	else { die "List of vcf files have not been provided or set before calling updateVCFs!\n" }
	unless (@{$vcffiles}==$npops){ die "Wrong number of genotype files specified!\n" }

	my @outputfiles;
	for my $pop (0..$#{$vcffiles}){	# print header.
		my ($filename,$path,$suffix)=fileparse($vcffiles->[$pop], qr/\.vcf\.b?gz/);
		my $outputfile=$outfolder . "/" . $filename . "_" . $outsuffix . ".vcf";
#		(my $outputfile=$vcffiles->[$pop])=~s/(.*)(\.vcf\.?.*)/${1}_${outsuffix}${2}/;
		push(@outputfiles, $outputfile);

		unless (-e $outputfile){
			open my $vcffile, "-|", "tabix -H $vcffiles->[$pop]" or die "Could not open vcf file $vcffiles->[$pop]!\n$!\n";
			open my $output, ">", $outputfile or die "Could not open output vcf file $outputfile!\n$!\n";
			while (<$vcffile>){ print $output $_ if (/^#/) }
			close($vcffile);
			close($output);
			}
		}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			for my $pop (0..$#{$vcffiles}){
				open my $vcffile, "-|", "tabix $vcffiles->[$pop] $chrom:$start-$end" or die "Could not open vcf file $vcffiles->[$pop]!\n$!\n";
				open my $output, ">>", $outputfiles[$pop] or die "Could not open output vcf file $outputfiles[$pop]!\n$!\n";
				my ($nsnps, $nupdated)=(0, 0);
				print "Updating genotypes for $chrom:$start-$end, population $pop.\n" if $verbose;
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						chomp;
						my @line=split(/\s/, $_);
						my $pos=$line[1]-1;	# The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
						++$nsnps;
						if ($self->{_isPhased}{$chrom}{$pos}[$pop]){
							my $gtstring=$self->{_Genotypes}{$chrom}{$pos}[$pop];
							IND: for my $ind (@{$samples->[$pop]}){
								my $allele1=substr($gtstring, 0, 1, "");
								my $allele2=substr($gtstring, 0, 1, "");
								my $col=$ind+8;
								if ($col>$#line){ die "Invalid sample index specified for population $pop!\n" }
								$line[$col]=~/([\.0-3])\/([\.0-3])/;
								my $wasmissing=($1 eq '.' || $2 eq '.');
								my $ismissing=($allele1 eq '.' || $allele2 eq '.');
								if ($wasmissing && !$ismissing){
									warn "Missing data can be imputed at $chrom:$pos (${1}/${2} vs. ${allele1}/${allele2})!\n" if $verbose;
									if ($impute){ $line[$col]=~s/\.\/\./${allele1}|${allele2}/ }
									else { $line[$col]=~s/\.\/\./\.|\./ }
									}
								elsif (!$wasmissing && !$ismissing){
									unless (($1 eq $allele1 && $2 eq $allele2) || ($1 eq $allele2 && $2 eq $allele1) ){
										warn "Error: phased genotype doesn't match previous genotype at $chrom:$pos (${1}/${2} vs. ${allele1}/${allele2})!\n";
										next IND;
										}
									$line[$col]=~s/[\.0-3]\/[\.0-3]/${allele1}|${allele2}/;
									}
								else {
									$line[$col]=~s/[\.0-3]\/[\.0-3]/${allele1}|${allele2}/;
									}
								}
							if (length($gtstring)){ die "Invalid length of genotype string at $chrom:$pos!\n" }
							++$nupdated;
							}
						print $output join("\t", @line), "\n" unless ($onlyphased && !$self->{_isPhased}{$chrom}{$pos}[$pop]);
						}
					}
				print STDERR "$chrom:$start-$end - vcf file: $outputfiles[$pop], number of SNPs: $nsnps, number of SNPs updated: $nupdated\n";
				close $vcffile;
				close $output;
				}
			}
		}
	return;
	}


sub addOutgroupSNPs{
	my ($self, $ref_loci, $ref_fasta, $vcffile, $allsitesfile, $outfile, $verbose)=@_;

	unless (-e $outfile){
		open my $input, "-|", "tabix -H $vcffile" or die "Could not open vcf file $vcffile!\n$!\n";
		open my $output, ">", $outfile or die "Could not open output vcf file $outfile!\n$!\n";
		while (<$input>){ print $output $_ if (/^#/) }
		close($input);
		close($output);
		}

	CHROM: for my $chrom ($ref_fasta->getSortedIDs() ){
		LOCUS: for my $ref_locus (sort { $a->{'start'}<=>$b->{'start'} } $ref_loci->allValues($chrom) ){
			my $start=$ref_locus->{'start'}+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_locus->{'end'}+1;
			my %snps;
			print "Updating genotypes for $chrom:$start-$end\n" if $verbose;

			open my $input, "-|", "tabix $vcffile $chrom:$start-$end" or die "Could not open vcf file $vcffile!\n$!\n";
			while (<$input>){
				if(/^\w+\s/){
					my $pos=( split(/\s/, $_) )[1];
					$snps{$pos}=$_;
					}
				}
			close $input;

			open my $allsites, "-|", "tabix $allsitesfile $chrom:$start-$end" or die "Could not open vcf file $allsitesfile!\n$!\n";
			while (<$allsites>){
				if(/^\w+\s/){
					my @line=split(/\s/, $_);
					if ($line[4] ne '.'){ $snps{ $line[1] } //= $_ }
					}
				}
			close $allsites;

			open my $output, ">>", $outfile or die "Could not open output vcf file $outfile!\n$!\n";
			for my $pos (sort {$a<=>$b} keys %snps){
				print $output $snps{$pos}
				}
			close $output;

			}
		}
	return;
	}


sub validatePhasing{
	my ($self, $ref_loci)=@_;

	my ($vcffiles, $samples);
	if ( defined $self->{_Files} && defined $self->{_Samples} ){
		$vcffiles=$self->{_Files};
		$samples=$self->{_Samples};
		}
	else { die "List of vcf files have not been provided or set before calling validatePhasing!\n" }
	unless (@{$vcffiles}==2){ die "Wrong number of genotype files specified! Works only with exactly two vcf files to compare.\n" }
	unless (@{ $samples->[0] }==@{ $samples->[1] }){ die "Unequal number of samples in the two vcf files to compare! (", scalar @{ $samples->[0] }, " vs. ", scalar @{ $samples->[1] }, ")\n" }
	my $nsamples=@{ $samples->[0] };

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		my @counts_by_dist=[] x $nsamples;
		my @total_by_dist=[] x $nsamples;
		my @prevgt=[] x $nsamples;	
		SITE: for my $seg (sort {$a<=>$b} keys %{ $self->{_Genotypes}{$chrom} } ){
			my $gtstring1=$self->{_Genotypes}{$chrom}{$seg}[0];
			my $gtstring2=$self->{_Genotypes}{$chrom}{$seg}[1];
			IND: for my $ind (0..$#{$samples->[0]}){
				my $allele11=substr($gtstring1, 0, 1, "");
				my $allele12=substr($gtstring1, 0, 1, "");
				my $allele21=substr($gtstring2, 0, 1, "");
				my $allele22=substr($gtstring2, 0, 1, "");
				next IND if ($allele11 eq '.' || $allele12 eq '.' || $allele21 eq '.' || $allele22 eq '.');
				unless (($allele11 eq $allele21 && $allele12 eq $allele22) || ($allele11 eq $allele22 && $allele12 eq $allele21) ){
					warn "Error: phased genotypes don't match at $chrom:$seg (${allele11}/${allele12} vs. ${allele21}/${allele22})!\n";
					next IND;
					}
				if ($allele11 ne $allele12 && $allele21 ne $allele22){	# only consider sites were an individual is heterozygous in both files.
					if (@{ $prevgt[$ind] }){
						my $preval11=$prevgt[$ind][1];
						my $preval12=$prevgt[$ind][2];
						my $preval21=$prevgt[$ind][3];
						my $preval22=$prevgt[$ind][4];
						my $distcat=($seg-$prevgt[$ind][0])/10;
						unless ($preval11 . $allele11 eq $preval21 . $allele21 || $preval11 . $allele11 eq $preval22 . $allele22){
							warn "Unequal phasing between the two files at $chrom:$seg (", $preval11 . $allele11, " vs. ", $preval21 . $allele21, ")!\n";
							++$counts_by_dist[$ind][$distcat];
							$prevgt[$ind]=[$seg, $allele11, $allele12, $allele21, $allele22];
							}
						++$total_by_dist[$ind][$distcat];
						}
					}
				}
			}

		print "Chromosome $chrom, missphasing by individual and distance:\n";
		my @missperind;
		for my $ind (0..$nsamples-1){
			for my $cat (0..$#{ $counts_by_dist[$ind] }){
				if ($total_by_dist[$ind][$cat]){
					$missperind[$ind][$cat]=$counts_by_dist[$ind][$cat] ? sprintf("%.3f", $counts_by_dist[$ind][$cat]/$total_by_dist[$ind][$cat]) : "0.000";
					}
				else { $missperind[$ind][$cat]="NA" }
				}
			print "Individual $ind: ", join(',', @missperind), "\n";
			}
		}
	return;
	}



sub printGenotypes{
	my ($self, $outfile)=@_;
	my $output;

	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	SCAFF: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Genotypes}}){
		POSITION: for my $pos (sort {$a<=>$b} keys %{$self->{_Genotypes}{$chrom}}){
			print $output "$chrom\t$pos";
			POP: for my $pop (@{$self->{_Genotypes}{$chrom}{$pos}}){
				my @data;
				for (my $i=0; $i<length($$pop); $i+=2){push @data, substr($$pop, $i, 1) . "/" . substr($$pop, $i+1, 1)}
				print $output join(',', @data), "\t";
				}
			print $output "\n";					
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}



1;

