package Haplotypes;

use strict;
use warnings;
use Classes::Misc;


sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	return $self;
	}


sub genotypestoFastPhase{
	my ($self, $outfolder, $ref_loci, $refseq, $ref_genotypes, $ref_popsizes, $ref_groups, $ref_poplabels, $ref_minind, $excluderefN, $snplimit, $overlap, $allelemode, $cluster, $starts, $iter, $nhaps, $qual, $inpute, $verbose)=@_;
	unless ($snplimit>2*$overlap){ die "SNP limit has to be larger than twice the specified overlap!\n" }
	unless ($overlap){ warn "Warning: overlap need to be larger than 0! Setting to 1.\n"; $overlap=1 }

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my (@haplotypes, @phasedhaps, @positions);
			for my $pop (0..$#{ $ref_popsizes }){
				$haplotypes[$pop][$_]="" foreach ( 0..2*$ref_popsizes->[$pop]-1 );
				$phasedhaps[$pop][$_]="" foreach ( 0..2*$ref_popsizes->[$pop]-1 );
				}
			my ($windowstart, $windowend)=($start, $end);
			my $withoverlap=0;

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless( defined $ref_genotypes->getValue($chrom, $pos) );
				next POSITION if ($excluderefN && $refbase eq 'N');	# if reference has a N at the current position, skip position.
				if ($ref_genotypes->NAlleles($chrom, $pos)>2){	# fastPHASE does not accept non-biallelic SNPs.
					warn "Warning: position $chrom:", $pos+1, " is not bi-allelic, skipping position!\n" if $verbose;
					next POSITION;
					}
				my $ancestral;
				if ($allelemode eq 'ancder'){
					$ancestral=$ref_genotypes->getAncestral($chrom, $pos);
					unless (defined $ancestral && $ancestral<4){
						warn "Warning: ancestral state not defined for $chrom:", $pos+1, ", skipping position!\n" if $verbose;
						next POSITION;
						} 
					}
				my (@states, %alleles);

				POP: for my $pop (0..$#{ $ref_popsizes }){
					my $validal=0;
					my $gtstring=$ref_genotypes->getValue($chrom, $pos, $pop);
					unless ( $gtstring && length($gtstring)==(2*$ref_popsizes->[$pop]) ){
						warn "Error: genotype string not defined or wrong length for $chrom:", $pos+1, ", population $pop, setting to missing!\n";
						$gtstring='.' x (2*$ref_popsizes->[$pop]);
						}
					INDEX: for my $index ( 0..2*$ref_popsizes->[$pop]-1 ){
						my $allele=substr($gtstring, 0, 1, "");	# chew up genotype string.
						if (defined $allele && $allele ne "."){
							my $state;
							if ($allelemode eq 'refalt'){ $state=$allele }
							elsif ($allelemode eq 'ancder'){
								if ($ancestral==0){ $state=$allele }	# keeps as it is.
								elsif ($allele==$ancestral){ $state=0 }	# switch allele states between ancestral and reference.
								elsif ($allele==0){ $state=$ancestral }
								else { $state=$allele }
								}
							else { $state=$ref_genotypes->getBase($chrom, $pos, $allele) }
							push @{ $states[$pop] }, $state;
							$alleles{$allele}=1;
							++$validal;
							}
						else { push @{ $states[$pop] }, "?" }
						}
					if ($validal<2*$ref_minind->[$pop]){	# not enough alleles callable in a population
						@{ $states[$pop] }=("?") x (2*$ref_popsizes->[$pop]);
						}
					}

				if (scalar(keys %alleles)>1){	# check if position is really polymorphic.
					push @positions, ($pos+1);
					for my $pop (0..$#states){
						unless (@{$states[$pop]}==2*$ref_popsizes->[$pop]){ die "Wrong number of allelic states for population $pop in array at position $chrom:$pos! ", scalar(@{$states[$pop]}), " vs. ", 2*$ref_popsizes->[$pop], "\n" }
						my $index=0;
						$haplotypes[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add polymorphic position to haplotype string.
						}
					}

				if (@positions==$snplimit){
					print "Submitting to PHASE ...\n";
					my $inputfile="phasing_" . $chrom . "_" . ($windowstart+1) . "_" . ($pos+1) . ".inp";
					my $ref_inphased=$self->callPHASE($outfolder, $inputfile, "subpoplabels.txt", "results_phasing", \@haplotypes, $ref_groups, \@positions, $cluster, $starts, $iter, $nhaps, $qual, 0);
					$self->concatHaplotypes(\@phasedhaps, $ref_inphased, $withoverlap ? $overlap : 0);
					splice(@positions, 0, $overlap) if ($withoverlap);	# remove overlap from position array.
					unless (@positions==length($phasedhaps[0][0]) ){
						die "Length of phased haplotype strings (", length($phasedhaps[0][0]), ") does not match length of position array (", scalar(@positions), ")!\n";
						}
					$ref_genotypes->phaseGenotypes($chrom, \@positions, \@phasedhaps, $inpute);

					for my $ref_pop (@haplotypes){	# reset haplotype array but keep overlap snps.
						for my $hap (@$ref_pop){
							unless (length($hap)==$snplimit){ die "Haplotye array for $chrom:", $windowstart+1, "-", $pos+1, " has wrong length before resetting! ", length($hap), " vs. $snplimit\n" }
							substr($hap, 0, $snplimit-$overlap, "");
							unless (length($hap)==$overlap){ die "Haplotye array for $chrom:", $windowstart+1, "-", $pos+1, " has not been reset correctly! ", length($hap), " vs. $overlap\n" }
							}						
						}
					@positions=splice(@positions, -$overlap);	# reset position array but keep overlap snps.
					unless (@positions==$overlap){ die "Position array for locus $chrom:", $windowstart+1, "-", $pos+1, " has not been reset correctly!\n" }
					$windowstart=$positions[0];	# take position of first snp as windowstart.
					$withoverlap=1;
					}
				}

			if (@positions){	# print last input file with number of snps < snplimit.
				print "Submitting to PHASE ...\n";
				my $inputfile="phasing_" . $chrom . "_" . ($windowstart+1) . "_" . ($end+1) . ".inp";
				my $ref_inphased=$self->callPHASE($outfolder, $inputfile, "subpoplabels.txt", "results_phasing", \@haplotypes, $ref_groups, \@positions, $cluster, $starts, $iter, $nhaps, $qual, 0);
				$self->concatHaplotypes(\@phasedhaps, $ref_inphased, $withoverlap ? $overlap : 0);
				splice(@positions, 0, $overlap) if ($withoverlap);	# remove overlap from position array.
				$ref_genotypes->phaseGenotypes($chrom, \@positions, \@phasedhaps, $inpute);
				}
			}
		}
	return;
	}


sub concatHaplotypes{
	my ($self, $ref_chromhaps, $ref_parthaps, $lenoverlap)=@_;
	POP: for my $pop ( 0..$#{ $ref_parthaps } ){
		unless ($lenoverlap){
			$ref_chromhaps->[$pop][$_]=$ref_parthaps->[$pop][$_] foreach (0..$#{ $ref_parthaps->[$pop] });
			next POP;
			}
		my $nind=@{ $ref_parthaps->[$pop] }/2;
		IND: for my $ind (0..$nind-1){
			my $parthap1=$ref_parthaps->[$pop][2*$ind];
			my $parthap2=$ref_parthaps->[$pop][2*$ind+1];
			my $overlap1=substr($ref_chromhaps->[$pop][2*$ind], -$lenoverlap);
			my $overlap2=substr($ref_chromhaps->[$pop][2*$ind+1], -$lenoverlap);
			print STDERR substr($parthap1, 0, $lenoverlap), "\n", substr($parthap2, 0, $lenoverlap), "\n$overlap1\n$overlap2\n";
			my ($notinverted, $inverted)=(0,0);
			POS: for my $pos (1..$lenoverlap){
				my $newbase1=substr($parthap1, 0, 1, "");
				my $newbase2=substr($parthap2, 0, 1, "");
				my $presentbase1=substr($overlap1, 0, 1, "");
				my $presentbase2=substr($overlap2, 0, 1, "");
				if ($newbase1 eq $presentbase1 && $newbase2 eq $presentbase2){ ++$notinverted unless ($newbase1 eq $newbase2) }
				elsif ($newbase1 eq $presentbase2 && $newbase2 eq $presentbase1){ ++$inverted }
				else { warn "Error: genotypes don't match!\n" }
				}
			if (length($overlap1) || length($overlap2) ){ warn "Error: wrong length of haplotype strings!\n" }
			if (!$inverted && !$notinverted){ warn "Warning: no heterozygous positions in overlap!\n" }
			elsif ($inverted && $notinverted){
				warn "Warning: unambiguous matching of haplotypes! Not inverted: $notinverted, inverted: $inverted.\n";
				}
			if ($notinverted>=$inverted){
				$ref_chromhaps->[$pop][2*$ind]=$parthap1;
				$ref_chromhaps->[$pop][2*$ind+1]=$parthap2;
				}
			else {
				$ref_chromhaps->[$pop][2*$ind]=$parthap2;
				$ref_chromhaps->[$pop][2*$ind+1]=$parthap1;
				}
			}
		}
	return;
	}


sub callPHASE{	# Attention: fastPHASE has a maximum of 500'000 characters limit per line in the input file!
	my ($self, $workdir, $infile, $subpopfile, $out_prefix, $ref_haplotypes, $ref_groups, $ref_positions, $cluster, $starts, $iter, $nhaps, $qualthreshold, $printpos, $dryrun)=@_;

	my %groupids;
	for my $group ( 0..$#{ $ref_groups } ){
		$groupids{$_}=$group foreach (@{ $ref_groups->[$group] });
		}
	unless (keys %groupids==@$ref_haplotypes){ die "Not all populations have a group membership defined!\n" }

	my $rowlength=400000;
	my $totsize=0;
	my (@popsizes, @labels);
	for my $pop (0..$#{ $ref_haplotypes } ){
		if (@{ $ref_haplotypes->[$pop] } % 2){ die "Number of haplotypes in population $pop is not divisible by 2!\n" }
		my $size=@{ $ref_haplotypes->[$pop] } / 2;
		push @popsizes, $size;
		$totsize+=$size;
		push @labels, ($groupids{$pop}+1) x $size;
		}

	chdir($workdir);

	open my $subpop, ">", $subpopfile or die "Could not open phase input file $subpopfile!\n";
	print $subpop join(' ', @labels);
	close $subpop;

	open my $output, ">", $infile or die "Could not open phase input file $infile!\n";
	print $output "$totsize\n";
	print $output scalar(@$ref_positions), "\n";
	if ($printpos){
		my $full_posstring="P " . join(" ", @$ref_positions) . "\n";
		if (length($full_posstring<=$rowlength) ){ print $output $full_posstring }
		else {
			print $output "P";
			my ($nchar, $index)=(0, 0);
			while ($index<@$ref_positions){
				my $posstring=" " . $ref_positions->[$index];
				$nchar+=length($posstring);
				if ($nchar<$rowlength){ print $output $posstring; ++$index }
				else {
					print $output "\n";
					$nchar=0;
					}
				}
			}
		}

	for my $ref_pop (@$ref_haplotypes){
		for my $hap (@$ref_pop){
			my $hapcopy=$hap;
			print $output substr($hapcopy, 0, $rowlength, ""), "\n" while (length($hapcopy)>0);
			}
		}
	close $output;

	my @clustersettings=split(" ", $cluster);
	my @args=("fastPHASE", "-n", @clustersettings, "-q${qualthreshold}", "-T${starts}", "-C${iter}", "-H${nhaps}", "-u${subpopfile}", "-o${out_prefix}", "$infile");
	my $ref_phased;

	unless ($dryrun){
		print STDERR "Calling PHASE with the following command line:\n", join(" ", @args), "\n";
		system(@args);
		if ($? == -1){
			die "PHASE failed to execute: $!\n";
			}
		elsif ($? & 127){
			printf STDERR "child died with signal %d, %s coredump\n",
			($? & 127), ($? & 128) ? 'with' : 'without';
			die "PHASE failed to execute: $!\n";
			}
		else {
			printf STDERR "child exited with value %d\n", $? >> 8;
			}
		$ref_phased=$self->readPHASE($workdir, $out_prefix, \@popsizes, ($qualthreshold>0) );
		for my $pop (0..$#popsizes){	# check correct length of sequence and haplotype strings.
			for my $ind ( 0..$#{$ref_phased->[$pop]} ){
				unless ( @$ref_positions==length($ref_phased->[$pop][$ind]) ){
					die "Haplotype string for population $pop, individual $ind has wrong length! ", scalar(@$ref_positions), " vs. ", length($ref_phased->[$pop][$ind]), "\n";
					}
				}
			}
		}

	return($ref_phased);
	}


sub readPHASE{
	my ($self, $workdir, $out_prefix, $ref_popsizes, $exclude)=@_;

	my @phased;
	my ($switch, $pop, $ind);
	my $phasefile="$out_prefix" . "_hapguess_switch.out";
	open my $input, "<", "$workdir/$phasefile" or die "Could not open phase file $workdir/$phasefile!\n$!\n";
	LINE: while (<$input>){
		if(/^BEGIN GENOTYPES/){ $switch=1 }
		elsif (/^END GENOTYPES/){ last LINE }
		elsif(/^#/ && $switch){
			if (defined $pop && $ind){ --$ind }
			elsif (defined $pop){
				unless (@{$phased[$pop]}==2*$ref_popsizes->[$pop]){ die "Population $pop has wrong number of haplotypes! ", scalar(@{$phased[$pop]}), " vs. ", 2*$ref_popsizes->[$pop], "\n" }
				++$pop;
				$ind=$ref_popsizes->[$pop]-1;
				}
			else { $pop=0; $ind=$ref_popsizes->[$pop]-1 }
			}
		elsif($switch && defined $pop){
			chomp;
			my @line=split(/\s+/, $_);
			my $haplotype;
			for my $pol (@line){
				if ($pol=~/\[(\S)\]/){ $haplotype.=$exclude ? "?" : $1 }
				else { $haplotype.=$pol }
				}
			push @{$phased[$pop]}, $haplotype;
			}
		}
	close $input;

	return(\@phased);
	}


sub printFastPHASE{
	my ($self, $outfolder, $out_prefix, $ref_positions, $ref_haplotypes, $id, $windowstart, $windowend)=@_;

	my $totsize=0;
	for my $pop ( 0..$#{ $ref_haplotypes } ){
		if ( @{ $ref_haplotypes->[$pop] } % 2 ){ die "Number of haplotypes in population $pop is not divisible by 2!\n" }
		$totsize+=( @{ $ref_haplotypes->[$pop] } / 2 );
		}

	print "Printing fastPhase input file for $id:$windowstart-$windowend ...\n";
	my $inputfile=$out_prefix . "_" . $id . "_" . $windowstart . "_" . $windowend . ".inp";
	open my $output, ">", "$outfolder/$inputfile" or die "Could not open phase input file $inputfile!\n";
	print $output "$totsize\n";
	print $output scalar(@$ref_positions), "\n";
	print $output "P ", join(" ", @$ref_positions), "\n";
	for my $ref_pop ( @{ $ref_haplotypes } ){
		print $output "$_\n" foreach (@$ref_pop);
		}
	close $output;
	return;
	}





1;

