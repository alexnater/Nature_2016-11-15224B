package Nexus;

use strict;
use warnings;
use Classes::Misc;


sub new{
	my $class=shift;
	my $self={};
	$self->{_Test}="This should appear!";
	bless $self, $class;
	return $self;
	}


sub setSubIndices{
	my ($self, $ref_subind, $haplotize)=@_;
	my @indices;
	for my $pop (@$ref_subind){
		if ($haplotize){ push @indices, [ map { 2*$_+int(rand(2)) } @$pop ] }
		else { push @indices, [ map { (2*$_, 2*$_+1) } @$pop ] }
		}
	$self->{_SubIndices}=\@indices;
	return;
	}


sub genotypestoNexus{
	my ($self, $regions, $refseq, $ref_genotypes, $ref_missing, $npops, $popsizes, $groups, $mincov, $minind, $haplotize, $polarize, $ref_regionstomask)=@_;

	print STDERR "Generating nexus data...\n";

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;
	my @haplotypes;
	my $totsize;
	my $totsegs=0;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($chrom) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($chrom) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $undefanc, $notbiallelic)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless (defined $ref_genotypes->getValue($chrom, $pos) );
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.

				my $ancestral;
				my $nalleles=$ref_genotypes->NAlleles($chrom, $pos);
				if ($nalleles>2){ ++$notbiallelic; next POSITION }
				if ($polarize){
					$ancestral=$ref_genotypes->getAncestral($chrom, $pos);
					unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POSITION }
					unless ( $ref_genotypes->isAllele($chrom, $pos, $ancestral) ){ ++$undefanc; next POSITION }
					}
				else { $ancestral=0 }
				my @allalleles;
				my $totminnumber=0;
				my $offset=$pos-$start;
				GROUP: for my $metapop ( 0..$#{$groups} ){
					for my $subpop ( @{$groups->[$metapop]} ){
						my @validalleles;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=defined $ref_missing ? $ref_missing->getValue($chrom, $start, $allind, $offset) : $ref_genotypes->getCoverage($chrom, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($haplotize){
								my $index=2*$ind+int(rand(2));
								my $allele=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
								unless (defined $allele && $allele ne '.'){ next IND }
								push @validalleles, $allele;
								}
							else {
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								++$index;
								my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								push @validalleles, ($allele1, $allele2);
								}
							}
						my $minnumber=$haplotize ? $minind->[$subpop] : 2*$minind->[$subpop];
						$totminnumber+=$minnumber;
						unless (@validalleles>=$minnumber){ next POSITION }
						my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $minnumber);
						my @recoded=map { $_ eq $ancestral ? 0 : 1 } @$ref_suballeles;
						push @allalleles, @recoded;
						}
					}
				if (defined $totsize){
					die "Unequal sample sizes $totsize vs. $totminnumber!\n" unless $totsize==$totminnumber;
					}
				else { $totsize=$totminnumber }

				my $allelesum=Misc::sum(@allalleles);
				next POSITION if ($allelesum==0 || $allelesum==$totminnumber);	# check if subsampled alleles are still polymorphic.
				POP: for my $metapop ( 0..$#{$groups} ){
					my $minnumber=$haplotize ? $groupsizes[$metapop] : 2*$groupsizes[$metapop];
					IND: for my $ind (0..$minnumber-1){
						die "Too few alleles in array!\n" unless @allalleles;
						$haplotypes[$metapop][$ind].=shift(@allalleles);
						}
					}
				die "Wrong number of alleles in array!\n" if @allalleles;
				++$totsegs;
				} continue { $prevpos=$pos }
			print STDERR "$chrom:$start-$end; hardmasked: $hardmasked, undefined ancestral state: $undefanc, not biallelic: $notbiallelic\n";
			}
		}

	$self->{_Haplotypes}=\@haplotypes;
	$self->{_TotalInd}=$totsize;
	$self->{_TotalSegs}=$totsegs;
	return;
	}


sub genotypestoMSMC{
	my ($self, $outfile, $regions, $refseq, $ref_genotypes, $ref_missing, $ref_popsizes, $ref_groups, $mincov, $ref_minhaps, $ref_regionstomask)=@_;

	print STDERR "Generating MSMC data...\n";

	my @popsum;
	my $rt=0;
	for my $popsize (@$ref_popsizes){
		push @popsum, $rt;
		$rt+=$popsize;
		}

	my @subindices;
	if (defined $self->{_SubIndices} ){ @subindices=@{ $self->{_SubIndices} } }
	else { 
		for my $popsize (@$ref_popsizes){ push @subindices, [ map { (2*$_, 2*$_+1) } (0..$popsize-1) ] }
		}
	print STDERR join(',', @$_), "\n" foreach (@subindices);
	my $totsize=0;
	$totsize+=@$_ foreach (@subindices);

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($chrom) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($chrom) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($validsince, $totsegs, $hardmasked)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				my $offset=$pos-$start;
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;

				my $allelestring;
				my %alleles;
				GROUP: for my $metapop ( 0..$#{ $ref_groups } ){
					my $validalleles=0;
					for my $subpop ( @{ $ref_groups->[$metapop] } ){
						INDEX: for my $index (@{ $subindices[$subpop] }){	# inefficient if not haplotized.
							my $allind=int($index/2)+$popsum[$subpop];
#							print STDERR "$subpop, $index, $popsum[$subpop], $allind\n";
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next POSITION }	# don't allow missing data based on coverage.
							if ($snp){
								my $allele=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
								if (defined $allele && $allele ne 'N'){ ++$validalleles; ++$alleles{$allele} }
								else { $allele='?' }
								$allelestring.=$allele;
								}
							}
						}
					}
				++$validsince;
				if ($snp){
					unless (length($allelestring)==$totsize){ die "Invalid size of allelestring!\n" }
					if (keys %alleles>1){
						++$totsegs;
						print $output "$chrom\t$pos\t$validsince\t$allelestring\n";
						$validsince=0;
						}
					}					
				} continue { $prevpos=$pos }
			print STDERR "$chrom:$start-$end; total segregating sites: $totsegs, hardmasked: $hardmasked\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub genotypestoARGweaver{
	my ($self, $outfile_prefix, $regions, $refseq, $ref_genotypes, $ref_missing, $ref_poplabels, $ref_popsizes, $ref_groups, $mincov, $ref_mingrouphaps, $minpropvalid)=@_;

	print STDERR "Generating ARGweaver data...\n";

	my @popsum;
	my $rt=0;
	for my $popsize (@$ref_popsizes){
		push @popsum, $rt;
		$rt+=$popsize;
		}

	my @subindices;
	if (defined $self->{_SubIndices} ){ @subindices=@{ $self->{_SubIndices} } }
	else { 
		for my $popsize (@$ref_popsizes){ push @subindices, [ map { (2*$_, 2*$_+1) } (0..$popsize-1) ] }
		}
	print STDERR join(',', @$_), "\n" foreach (@subindices);
	my $totsize=0;
	$totsize+=@$_ foreach (@subindices);

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $maskstart=$start;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my ($totsegs, $validsites)=(0, 0);
			my $prevpos=$start;
			my $outfile=$outfile_prefix . "_" . $chrom . ":" . $start . "-" . $end . ".sites";
			my $maskfile=$outfile_prefix . "_" . $chrom . ":" . $start . "-" . $end . "_mask.bed";
			open my $output, ">", $outfile or die "Could not open output file!\n";
			open my $mask, ">", $maskfile or die "Could not open mask file!\n";
			print $output "NAMES";
			for my $subpop (0..$#subindices){
				print $output "\t" . $ref_poplabels->[$subpop] . "_hap" . $_ foreach (@{ $subindices[$subpop] })
				}
			print $output "\nREGION\t", $chrom, "\t", $start+1, "\t", $end+1, "\n";

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );
				next POSITION if ($refbase eq 'N');
				my $allelestring;
				my %alleles;
				my $offset=$pos-$start;

				GROUP: for my $metapop ( 0..$#{ $ref_groups } ){
					my $validgrouphaps=0;
					for my $subpop ( @{ $ref_groups->[$metapop] } ){
						INDEX: for my $index (@{ $subindices[$subpop] }){	# inefficient if not haplotized.
							my $allind=int($index/2)+$popsum[$subpop];
#							print STDERR "$subpop, $index, $popsum[$subpop], $allind\n";
							my $allele='N';
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							if (defined $coverage && $coverage>=$mincov){
								++$validgrouphaps;
								if (defined $ref_genotypes->getValue($chrom, $pos) ){
									$allele=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
									if (defined $allele && $allele ne 'N'){ ++$alleles{$allele} }
									}
								else { $allele=$refbase; ++$alleles{$allele} }
								}
							$allelestring.=$allele;								
							}
						}
					unless ($validgrouphaps>=$ref_mingrouphaps->[$metapop]){ next POSITION }
					}

				unless (length($allelestring)==$totsize){ die "Invalid size of allelestring!\n" }
				my $nalleles=scalar(keys %alleles);
				if ($nalleles){
					++$validsites;
					if ($pos-$maskstart){ print $mask "$chrom\t$maskstart\t$pos\n" }	# output in bed format!
					$maskstart=$pos+1;
					}
				if ($nalleles>1){
					++$totsegs;
					print $output $pos+1, "\t$allelestring\n";
					}					
				}
			print $mask "$chrom\t$maskstart\t", $end+1, "\n" if ($end-$maskstart>0);	# output in bed format!
			close $output;
			close $mask;
			unless ($validsites/$locuslength>=$minpropvalid){ rename $outfile, $outfile . "_invalid" }
			print STDERR "$chrom:$start-$end; total valid sites: $validsites, total segregating sites: $totsegs, hardmasked: ", $locuslength-$validsites, "\n";
			}
		}
	return;
	}


sub genotypestoIBS{
	my ($self, $outfolder, $regions, $refseq, $ref_genotypes, $ref_missing, $ref_popsizes, $ref_groups, $ref_grouplabels, $mincov, $ref_mingrouphaps)=@_;

	print STDERR "Generating IBS data...\n";

	my @popsum;
	my $rt=0;
	for my $popsize (@$ref_popsizes){
		push @popsum, $rt;
		$rt+=$popsize;
		}

	my @subindices;
	if (defined $self->{_SubIndices} ){ @subindices=@{ $self->{_SubIndices} } }
	else { 
		for my $popsize (@$ref_popsizes){ push @subindices, [ map { (2*$_, 2*$_+1) } (0..$popsize-1) ] }
		}
	print STDERR join(',', @$_), "\n" foreach (@subindices);

	my @groupsizes;
	my $totsize=0;
	for my $metapop ( 0..$#{ $ref_groups } ){
		for my $subpop ( @{ $ref_groups->[$metapop] } ){
			$groupsizes[$metapop]+=@{ $subindices[$subpop] };
			$totsize+=@{ $subindices[$subpop] };
			}
		}

	# open the population-specific output and mask files.
	my @outputs;
	my @ids=sort { Misc::expand($a) cmp Misc::expand($b) } $regions->allKeys();
	my $idstring=@ids>1 ? $ids[0] . "_" . $ids[-1] : $ids[0];
	for my $metapop ( 0..$#{ $ref_groups } ){
		my $outfile=$outfolder . "/" . $ref_grouplabels->[$metapop] . "_" . $idstring . ".popdata";
		open $outputs[$metapop], ">>", $outfile or die "Could not open output file $outfile!\n";
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		my ($chromnr, $suffix)=$chrom=~/\D*(\d+)(a|b)?/i;
		$chromnr=23 if (defined $suffix && $suffix eq 'b');
		my $maskfile=$outfolder . "/chrom" . $chromnr . ".txt";
		open my $mask, ">>", $maskfile or die "Could not open mask file $maskfile!\n";
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $maskstart=$start;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my ($totsegs, $validsites)=(0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );
				next POSITION if ($refbase eq 'N');
				my @allelestrings;
				my %alleles;
				my $offset=$pos-$start;

				GROUP: for my $metapop ( 0..$#{ $ref_groups } ){
					my $validgrouphaps=0;
					for my $subpop ( @{ $ref_groups->[$metapop] } ){
						INDEX: for my $index (@{ $subindices[$subpop] }){	# inefficient if not haplotized.
							my $allind=int($index/2)+$popsum[$subpop];
#							print STDERR "$subpop, $index, $popsum[$subpop], $allind\n";
							my $allele='N';
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							if (defined $coverage && $coverage>=$mincov){
								++$validgrouphaps;
								if (defined $ref_genotypes->getValue($chrom, $pos) ){
									$allele=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
									if (defined $allele && $allele ne 'N'){ ++$alleles{$allele} }
									}
								else { $allele=$refbase; ++$alleles{$allele} }
								}
							$allelestrings[$metapop].=$allele;								
							}
						}
					unless (length($allelestrings[$metapop])==$groupsizes[$metapop]){ die "Invalid size of allelestring!\n" }
					unless ($validgrouphaps>=$ref_mingrouphaps->[$metapop]){ next POSITION }
					}

				my $nalleles=scalar(keys %alleles);
				if ($nalleles){
					++$validsites;
					if ($pos-$maskstart){ print $mask $maskstart+1, "\t$pos\n" }	# output not in bed format.
					$maskstart=$pos+1;
					}
				if ($nalleles>1){
					++$totsegs;
					print { $outputs[$_] } "$chromnr\t", $pos+1, "\t$refbase\t$allelestrings[$_]\n" foreach (0..$#allelestrings);
					}					
				}
			print $mask $maskstart+1, "\t", $end+1, "\n" if ($end-$maskstart>0);	# output not in bed format.
			print STDERR "$chrom:$start-$end; total valid sites: $validsites, total segregating sites: $totsegs, hardmasked: ", $locuslength-$validsites, "\n";
			}
		close $mask;
		}
	close $_ foreach (@outputs);
	return;
	}


sub genotypestoGPhocs{
	my ($self, $outfile, $regions, $ref_fasta, $ref_genotypes, $ref_missing, $ref_minind, $mincov, $minprop_valid, $subsample)=@_;

	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum);
	my $totsize=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $totsize;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$totsize+=$temp;
		print STDERR "$temp, $totsize\n";
		}

	my $locusnumber=0;
#	my $nloci=0;
#	$nloci+=$regions->getSize($_) foreach ( $regions->allKeys() );

	open my $output, ">", $outfile or die "Could not open output file!\n";
#	print $output "$nloci\n\n";

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ($regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my @seq;
			my @validsites_per_pop=(0) x @popsizes;
			my $printedsites=0;

			my (@indlist, $subtotsize);
			if ( $subsample && defined $self->getreftoIndividuals('locus', $chrom, $locus) ){
				my $iterator=$self->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				for my $pop (0..$#popsizes){
					my %indid=map { int($_/2) => 1 } @{ $iterator->[$pop] };	# determine individual id from index and remove duplicates if not haplotized.
					push @indlist, [ sort {$a<=>$b} keys %indid ];
					$subtotsize+=@{ $indlist[$pop] };
					}
				}
			else {
				for my $popsize (@popsizes){ push @indlist, [ 0..$popsize-1 ] }
				$subtotsize=$totsize;
				}

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($refbase eq 'N'){	# if reference has a N at the current position, skip position.
#					for my $pop (0..$#popsizes){ $seq[$pop][$_].='N' foreach (0..$#{ $indlist[$pop] }) }
					next POSITION;
					}
				my $offset=$pos-$start;
				my @states_by_pop=[] x @popsizes;
				my $valid=0;

				POP: for my $pop (0..$#popsizes){
					my $validind=0;
					my @states;
					IND: for my $ind (0..$#{ $indlist[$pop] }){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						unless (defined $coverage && $coverage>=$mincov){ push(@states, 'N'); next IND }	# missing data, replace reference base at offset with N.
						if (defined $ref_genotypes->getValue($chrom, $pos) ){	# if position is variable and covered, obtain individual genotypes.
							my $index=2*$ind;
							my $allele1=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
							unless (defined $allele1 && $allele1 ne 'N'){ push(@states, 'N'); next IND }
							my $allele2=$ref_genotypes->getBase($chrom, $pos, $pop, ++$index);
							unless (defined $allele2 && $allele2 ne 'N'){ push(@states, 'N'); next IND }
							if ($allele1 eq $allele2){ push(@states, $allele1) }
							else { push( @states, Misc::iupac($allele1, $allele2) ) }
							}
						else {	# if position not variable but covered, retain reference base.
							push(@states, $refbase);
							}
						++$validind;
						}
					unless( @states==@{ $indlist[$pop] }){ die "Wrong length of states array at $chrom:$pos for population $pop! (", scalar(@states), " vs. ", scalar(@{ $indlist[$pop] }), ")\n" }
					$states_by_pop[$pop]=\@states;
					if ($validind>=$ref_minind->[$pop]){ ++$validsites_per_pop[$pop]; $valid=1 }	# enough individuals callable in at least one population.
					}

				if ($valid){
					for my $pop (0..$#popsizes){
						$seq[$pop][$_].=$states_by_pop[$pop][$_] foreach (0..$#{ $states_by_pop[$pop] });
						}
					++$printedsites;
					}
				}

			for my $pop (0..$#popsizes){
				unless ($validsites_per_pop[$pop]/$locuslength>=$minprop_valid){
					print STDERR "Warning: Locus $chrom:$start-$end, population $pop has not enough valid sites (only $validsites_per_pop[$pop] out of $locuslength)\n";
					next LOCUS;
					}
				}
			print STDERR "Locus $chrom:$start-$end, all population have enough valid sites.\n";
			print $output $chrom, "_locus", ++$locusnumber, " $subtotsize $printedsites\n";
			for my $pop (0..$#popsizes){ 
				for my $ind ( 0..$#{ $seq[$pop] } ){
					unless ( length($seq[$pop][$ind])==$printedsites){ die "Wrong length of sequence string for $chrom:$start-$end in population $pop for individual $ind! (", length($seq[$pop][$ind]), " vs. $printedsites)\n" }
					print $output "pop${pop}_${ind}\t$seq[$pop][$ind]\n";
					}
				}
			print $output "\n";
			}
		}
	close $output;

	open my $header, ">", "${outfile}.header" or die "Could not open header file!\n";
	print $header "$locusnumber\n";
	close $header;

	return;
	}


sub genotypesforHtest{
	my ($self, $outfile, $regions, $refseq, $ref_genotypes, $ref_missing, $npops, $ref_popsizes, $ref_groups, $mincov, $ref_minind, $chrommap, $phased, $ref_regionstomask)=@_;

	print STDERR "Generating genotype data ...\n";

	my @popsum;
	my $rt=0;
	foreach (@$ref_popsizes){ push @popsum, $rt; $rt+=$_ }
	my @groupsizes=map { Misc::sum(@{ $ref_popsizes }[@$_]) } (@$ref_groups);
	my @minnumber=$phased ? map { 2*$_ } @$ref_minind : @$ref_minind;
	my $totsize=$phased ? 2*Misc::sum(@groupsizes) : Misc::sum(@groupsizes);
	my @miss=$phased ? ('N', 'N') : ('N');

	open my $output, ">>", $outfile or die "Could not open output file $outfile!\n";

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $nsnps, $filtered)=(0, 0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless (defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my @alleles=$ref_genotypes->getAlleles($id, $pos);
				if (@alleles<2){ next POSITION }	# only consider sites with exactly two alleles.
				elsif (@alleles>2){ ++$notbiallelic; next POSITION }
				my @genotypes;
				my $offset=$pos-$start;
				GROUP: for my $metapop ( 0..$#{$ref_groups} ){
					for my $subpop ( @{$ref_groups->[$metapop]} ){
						my $validgenotypes=0;
						IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=defined $ref_missing ? $ref_missing->getValue($id, $start, $allind, $offset) : $ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ push @genotypes, @miss; next IND }
							my $index=2*$ind;
							my $allele1=$ref_genotypes->getBase($id, $pos, $subpop, $index);
							unless (defined $allele1 && $allele1 ne '.'){ push @genotypes, @miss; next IND }
							my $allele2=$ref_genotypes->getBase($id, $pos, $subpop, ++$index);
							unless (defined $allele2 && $allele2 ne '.'){ push @genotypes, @miss; next IND }
							my @gt;
							if ($phased){ @gt=($allele1, $allele2) }
							elsif ($allele1 eq $allele2){ @gt=($allele1) }
							else { @gt=('.') }
							push @genotypes, @gt;
							++$validgenotypes;
							}
						unless ( $validgenotypes>=$minnumber[$subpop] ){ ++$filtered; next POSITION }
						}
					}
				unless (@genotypes==$totsize){ warn "Wrong number of genotypes in array at $id:$pos, ", scalar(@genotypes), " vs. $totsize, skipping position!\n"; next POSITION }
				my ($chrom, $chrompos);
				if (defined $chrommap){
					($chrom, $chrompos)=$chrommap->findPosonChrom($id, $pos);
					}
				else {
					($chrom, $chrompos)=($id, $pos);
					}
				print $output "${chrompos},", join(',', @genotypes), "\n";
				++$nsnps;
				} continue { $prevpos=$pos }
			print STDERR "$id:$start-$end; number of SNPs: $nsnps, hardmasked: $hardmasked, not biallelic: $notbiallelic, filtered: $filtered\n";
			}
		}

	close $output;
	return;
	}


sub genotypestoPCA{
	my ($self, $outfile, $snpposfile, $regions, $refseq, $ref_genotypes, $ref_missing, $npops, $popsizes, $groups, $mincov, $minind, $chrommap, $ref_regionstomask)=@_;

	print STDERR "Generating nexus data ...\n";

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my @groupsizes=map { Misc::sum(@{$popsizes}[ @$_ ]) } @$groups;
	my $totsize=Misc::sum(@groupsizes);

	open my $output, ">>", $outfile or die "Could not open output file $outfile!\n";
	open my $snppos, ">>", $snpposfile or die "Could not open output file $snpposfile!\n";

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $nsnps)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless (defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my @alleles=$ref_genotypes->getAlleles($id, $pos);
				if (@alleles<2){ next POSITION }	# only consider sites with exactly two alleles.
				elsif (@alleles>2){ ++$notbiallelic; next POSITION }
				my @genotypes;
				my $offset=$pos-$start;
				GROUP: for my $metapop ( 0..$#{$groups} ){
					for my $subpop ( @{$groups->[$metapop]} ){
						my $validgenotypes=0;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=defined $ref_missing ? $ref_missing->getValue($id, $start, $allind, $offset) : $ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ push @genotypes, 9; next IND }
							my $index=2*$ind;
							my $allele1=$ref_genotypes->getValue($id, $pos, $subpop, $index);
							unless (defined $allele1 && $allele1 ne '.'){ push @genotypes, 9; next IND }
							my $allele2=$ref_genotypes->getValue($id, $pos, $subpop, ++$index);
							unless (defined $allele2 && $allele2 ne '.'){ push @genotypes, 9; next IND }
							my $genotype;
							if ($allele1 eq $allele2){
								if ($allele1 eq $alleles[0]){ $genotype=0 }
								elsif ($allele1 eq $alleles[1]){ $genotype=2 }
								else { warn "Invalid allele $allele1 at $id:$pos, skipping position!\n"; next POSITION }
								}
							else { $genotype=1 }
							push @genotypes, $genotype;
							++$validgenotypes;
							}
						unless ( $validgenotypes>=$minind->[$subpop] ){ next POSITION }
						}
					}
				unless (@genotypes==$totsize){ warn "Wrong number of genotypes in array at $id:$pos, ", scalar(@genotypes), " vs. $totsize, skipping position!\n"; next POSITION }
				my ($chrom, $chrompos);
				if (defined $chrommap){
					($chrom, $chrompos)=$chrommap->findPosonChrom($id, $pos);
					}
				else {
					($chrom, $chrompos)=($id, $pos);
					}
				print $output "$chrom\t$chrompos\t", join(' ', @genotypes), "\n";
				print $snppos "$id\t$pos\n";
				++$nsnps;
				} continue { $prevpos=$pos }
			print STDERR "$id:$start-$end; number of SNPs: $nsnps, hardmasked: $hardmasked, not biallelic: $notbiallelic\n";
			}
		}

	close $output;
	close $snppos;
	return;
	}


sub genotypestoREHH{
	my ($self, $mapfile, $hapfile_prefix, $regions, $refseq, $ancref, $ref_genotypes, $npops, $popsizes, $minind, $mincov)=@_;

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my $totsize=Misc::sum(@$popsizes);

	my @haps;
	open my $map, ">>", $mapfile or die "Could not open map file $mapfile!\n";

	ID: for my $id ( $regions->allKeys() ){
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($id) ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($id, $start, $end, 1);
			my ($hardmasked, $undefanc, $notbiallelic, $nsnps)=(0, 0, 0, 0);

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1, "") );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.
				next POSITION unless (defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.

				my $ancallele=$ref_genotypes->getAlleleforBase($id, $pos, $ancbase);
				my ($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($id, $pos, $ancallele);
				unless (defined $ancestral){ ++$undefanc; next POSITION }
				unless (@derived){ next POSITION }	# only consider sites with exactly two alleles.
				if (@derived>1){ ++$notbiallelic; next POSITION }	# maximum one derived allele permitted.
				my @alleles_bypop;

				POP: for my $pop (0..$npops-1){
					my @alleles;
					my $validgenotypes=0;
					IND: for my $ind (0..$popsizes->[$pop]-1){
						my $coverage=$ref_genotypes->getCoverage($id, $pos, $pop, $ind);
						unless (defined $coverage && $coverage>=$mincov){ push @alleles, (0, 0); next IND }
						my $index=2*$ind;
						my $allele1=$ref_genotypes->getValue($id, $pos, $pop, $index);
						unless (defined $allele1 && $allele1 ne '.'){ push @alleles, (0, 0); next IND }
						my $allele2=$ref_genotypes->getValue($id, $pos, $pop, ++$index);
						unless (defined $allele2 && $allele2 ne '.'){ push @alleles, (0, 0); next IND }
						push @alleles, ( ($allele1==$ancestral) ? 1 : 2, ($allele2==$ancestral) ? 1 : 2);
						++$validgenotypes;
						}
					unless (@alleles==2*$popsizes->[$pop]){ die "Wrong number of bases in array at $id:$pos for population $pop, ", scalar(@alleles), " vs. ", 2*$popsizes->[$pop], ", skipping position!\n" }
					unless ( $validgenotypes>=$minind->[$pop] ){ next POSITION }
					push @alleles_bypop, \@alleles;
					}

				for my $pop (0..$npops-1){
					my @alleles=@{ $alleles_bypop[$pop] };
					$haps[$pop][$_].=" " . $alleles[$_] foreach (0..$#alleles);
					}
				my $snpname="SNP" . ++$nsnps;
				my $chromnr=$1 if $id=~/\D*(\d+)(a|b)?/i;
				if (defined $2 && $2 eq 'b'){ $chromnr=23 }
				unless ($chromnr){ die "Unable to assign chromosome number for $id!\n" }
				print $map "$snpname $chromnr ", $pos+1, " 1 2\n";
				}

			print STDERR "$id:$start-$end; number of SNPs: $nsnps, hardmasked: $hardmasked, not biallelic: $notbiallelic\n";
			}
		}
	close $map;

	for my $i (0..$npops-1){
		open my $output, ">", $hapfile_prefix . ".pop" . $i or die "Could not open output file!\n";
		print $output "hap${_}${haps[$i][$_]}\n" foreach (0..$#{ $haps[$i] });
		close $output;
		}

	return;
	}


sub genotypestoPLINK{
	my ($self, $mapfile, $regions, $refseq, $ref_genotypes, $npops, $popsizes, $mincov, $minind, $ref_regionstomask)=@_;

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my $totsize=Misc::sum(@$popsizes);

	open my $map, ">>", $mapfile or die "Could not open map file $mapfile!\n";

	ID: for my $id ( $regions->allKeys() ){
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($id) ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring;
			if (defined $refseq){ $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1) }
			else { my $tempstring="A" x $locuslength; $ref_seqstring=\$tempstring }
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $nsnps)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless (defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my @alleles=$ref_genotypes->getAlleles($id, $pos);
				if (@alleles<2){ next POSITION }	# only consider sites with exactly two alleles.
				elsif (@alleles>2){ ++$notbiallelic; next POSITION }
				my @bases;
				POP: for my $pop (0..$npops-1){
					my $validgenotypes=0;
					IND: for my $ind (0..$popsizes->[$pop]-1){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_genotypes->getCoverage($id, $pos, $pop, $ind);
						unless (defined $coverage && $coverage>=$mincov){ push @bases, (0, 0); next IND }
						my $index=2*$ind;
						my $allele1=$ref_genotypes->getBase($id, $pos, $pop, $index);
						unless (defined $allele1 && $allele1 ne 'N'){ push @bases, (0, 0); next IND }
						my $allele2=$ref_genotypes->getBase($id, $pos, $pop, ++$index);
						unless (defined $allele2 && $allele2 ne 'N'){ push @bases, (0, 0); next IND }
						push @bases, ($allele1, $allele2);
						++$validgenotypes;
						}
					unless ( $validgenotypes>=$minind->[$pop] ){ next POSITION }
					}
				unless (@bases==2*$totsize){ warn "Wrong number of bases in array at $id:$pos, ", scalar(@bases), " vs. $totsize, skipping position!\n"; next POSITION }
				$self->{_PLINK}[$_].=" " . shift(@bases) . " " . shift(@bases) foreach (0..$totsize-1);
				if (@bases){ die "Wrong length of base array for $id:$pos!\n" }
				print $map "$id\tSNP\t0.0\t$pos\n";
				++$nsnps;
				} continue { $prevpos=$pos }
			$self->{_TotalSNPs}+=$nsnps;
			print STDERR "$id:$start-$end; number of SNPs: $nsnps, hardmasked: $hardmasked, not biallelic: $notbiallelic\n";
			}
		}
	close $map;
	return;
	}

sub genotypestoPLINK_byChrom{
	my ($self, $outfile, $mapfile, $regions, $refseq, $ref_genotypes, $npops, $popsizes, $mincov, $minind, $chrommap, $ref_regionstomask)=@_;

	print STDERR "Generating PLINK file ...\n";

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my $totsize=Misc::sum(@$popsizes);

	my %snps;
	my $prevchrom;

	ID: for my $id ( $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allValues($id) ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $nsnps)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless (defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites. This has to be placed after the refseq chewing!
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my @alleles=$ref_genotypes->getAlleles($id, $pos);
				if (@alleles<2){ next POSITION }	# only consider sites with exactly two alleles.
				elsif (@alleles>2){ ++$notbiallelic; next POSITION }
				my @bases;
				POP: for my $pop (0..$npops-1){
					my $validgenotypes=0;
					IND: for my $ind (0..$popsizes->[$pop]-1){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_genotypes->getCoverage($id, $pos, $pop, $ind);
						unless (defined $coverage && $coverage>=$mincov){ push @bases, (0, 0); next IND }
						my $index=2*$ind;
						my $allele1=$ref_genotypes->getBase($id, $pos, $pop, $index);
						unless (defined $allele1 && $allele1 ne 'N'){ push @bases, (0, 0); next IND }
						my $allele2=$ref_genotypes->getBase($id, $pos, $pop, ++$index);
						unless (defined $allele2 && $allele2 ne 'N'){ push @bases, (0, 0); next IND }
						push @bases, ($allele1, $allele2);
						++$validgenotypes;
						}
					unless ( $validgenotypes>=$minind->[$pop] ){ next POSITION }
					}
				unless (@bases==2*$totsize){ warn "Wrong number of bases in array at $id:$pos, ", scalar(@bases), " vs. $totsize, skipping position!\n"; next POSITION }
				my ($chrom, $chrompos);
				if (defined $chrommap){
					($chrom, $chrompos)=$chrommap->findPosonChrom($id, $pos);
					}
				else {
					($chrom, $chrompos)=($id, $pos);
					}
				$prevchrom//=$chrom;
				if ($prevchrom ne $chrom){ die "Function assumes that all regions belong to the same chromosome!\n" }
				$snps{$chrompos}=\@bases;
				++$nsnps;
				} continue { $prevpos=$pos }
			print STDERR "$id:$start-$end; number of SNPs: $nsnps, hardmasked: $hardmasked, not biallelic: $notbiallelic\n";
			}
		}

	my @sortedsnps;
	open my $map, ">>", $mapfile or die "Could not open map file $mapfile!\n";
	for my $pos (sort {$a<=>$b} keys %snps){
		print $map "$prevchrom\tSNP\t0.0\t$pos\n";
		$sortedsnps[$_].=" " . shift(@{ $snps{$pos} }) . " " . shift(@{ $snps{$pos} }) foreach (0..$totsize-1);
		if (@{ $snps{$pos} }){ die "Wrong length of base array for $prevchrom:$pos!\n" }
		}
	close $map;

	open my $output, ">>", $outfile or die "Could not open output file $outfile!\n";
	print $output "ID1 ID2 0 0 0 0", $_, "\n" foreach (@sortedsnps);
	close $output;
	return;
	}


sub genotypestoNexusOld{
	my ($self, $ref_groups, $ref_loci, $ref_genotypes, $ref_missing, $haplotize, $mode)=@_;
	my $npops=$ref_missing->getNPop();
	my @haplotypes;
	my $totsize;
	my $totsegs=0;

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;

			my $iterator;
			if ($mode eq 'locus' && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus)){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				print STDERR "Locus $chrom: $locus in locus mode\n";
				}
			elsif ($mode eq 'position'){ print STDERR "Locus $chrom: $locus in position mode\n" }
			else { die "Information about missing data mode not available!\n" }

			SITE: for my $seg ( sort {$a<=>$b} $ref_genotypes->allKeys($chrom) ){
				if ($seg<$start){ next SITE }
				elsif ($seg>=$start && $seg<=$end){
					next SITE unless ($ref_genotypes->NAlleles($chrom, $seg)==2);
					next SITE if ($ref_genotypes->{_Alleles}{$chrom}{$seg}>3);

					if ($mode eq 'position'){$iterator=$ref_loci->getreftoIndividuals('position', $chrom, $seg)}
					if ($npops!=scalar(@{$iterator})){die "Wrong information about number of populations!\n"}
					my $offset=$seg-$start;
					my $outstate=$ref_genotypes->getAncestral($chrom, $seg);
					$outstate=4 unless defined $outstate;
					my $rt=0; my $subrt=0;

					GROUP: for my $metapop ( 0..$#{$ref_groups} ){
						my $metaind=0;
						for my $subpop ( @{$ref_groups->[$metapop]} ){
							my $gtstring=$ref_genotypes->getValue($chrom, $seg, $subpop);
							my @newiterator;
							if ($haplotize){
								for (my $ind=0;$ind<$#{$iterator->[$subpop]};$ind=$ind+2){
									push @newiterator, $iterator->[$subpop][$ind]+int(rand(2));
									}
								}
							else {@newiterator=( @{$iterator->[$subpop]} )}
							for my $ind (0..$#newiterator){
								my $index=$newiterator[$ind];
								my $allind=int($index/2)+$rt;
								my $gt=substr($gtstring, $index, 1);
								my $miss=$ref_missing->getValue($chrom, $start, $allind, $offset);
								if (!defined $gt){ $haplotypes[$metapop][$metaind].='?' }
								elsif (!defined $miss){ $haplotypes[$metapop][$metaind].='?' }
								elsif ($miss==0 || $gt eq '.'){ $haplotypes[$metapop][$metaind].='?' }
								else { $haplotypes[$metapop][$metaind].=$gt }	# ancestral state needed?
								} continue { ++$metaind; ++$subrt }
							$rt+=$ref_missing->getNInd($subpop);
							}
						}
					if (defined $totsize){
						die "Unequal sample sizes $totsize vs. $subrt!\n" unless $totsize==$subrt;
						}
					else { $totsize=$subrt }
					++$totsegs;
					}
				else { last SITE }
				}
			}
		}
	$self->{_Haplotypes}=\@haplotypes;
	$self->{_TotalInd}=$totsize;
	$self->{_TotalSegs}=$totsegs;
	return;
	}


sub summarizeNexus{
	my ($self, $folder, $prefix, $glob, @idlist)=@_;

	my @filelist;
	if ($glob){ @filelist=<$folder/${prefix}_*.nex> }
	else {
		for my $id (@idlist){
			my $filename="$folder" . "/" . "$prefix" . '_' . "$id" . '.nex';
			push @filelist, $filename;
			}
		}
#	print STDERR join(', ', @filelist), "\n";

	my $totalind;

	for my $filename (@filelist){
		print STDERR "Summarizing Nexus file $filename.\n";
		open my $fh, "<", $filename or die "Could not open $filename!\n";
		my ($pop, $ind, $length);
		my $rt=0;
		LINE: while (<$fh>){
			if (/^pop(\d+)_ind(\d+)/){
				$pop=$1; $ind=$2;
				chomp $_;
				my $haplotype=( split(/\s+/, $_) )[1];
				if ( defined $length && $length!=length($haplotype) ){ die "Wrong length of haplotype in file $filename: Pop: $pop, Ind: $ind!\n" }
				$length=length($haplotype) unless defined $length;
				$self->{_Haplotypes}[$pop][$ind].=$haplotype;
				++$rt;
				}
			}
		close $fh;
		if ( defined $totalind && $totalind!=$rt ){ die "Wrong number of individuals in file $filename!\n" }
		$totalind=$rt unless defined $totalind;
		$self->{_TotalSegs}+=$length;		
		}
	$self->{_TotalInd}=$totalind;
	return;
	}


sub summarizeNexusOld{
	my ($self, $folder, $prefix, $glob, @idlist)=@_;

	my @filelist;
	if ($glob){ @filelist=<$folder/${prefix}_*.nex> }
	else {
		for my $id (@idlist){
			my $filename="$folder" . "/" . "$prefix" . '_' . "$id" . '.nex';
			push @filelist, $filename;
			}
		}
#	print STDERR join(', ', @filelist), "\n";

	my $totalind;

	for my $filename (@filelist){
		print STDERR "Summarizing Nexus file $filename.\n";
		open my $fh, "<", $filename or die "Could not open $filename!\n";
		my ($pop, $ind, $length);
		my $rt=0;
		my $switch=0;
		LINE: while (<$fh>){
			if ($switch){
				chomp $_;
				my $haplotype=$_;
				if ( defined $length && $length!=length($haplotype) ){ die "Wrong length of haplotype in file $filename: Pop: $pop, Ind: $ind!\n" }
				$length=length($haplotype) unless defined $length;
				$self->{_Haplotypes}[$pop][$ind].=$haplotype;
				$switch=0;
				}

			elsif (/^pop(\d+)_ind(\d+)/){
				$pop=$1; $ind=$2;
				++$rt;
				$switch=1;
				}
			}
		close $fh;
		if ( defined $totalind && $totalind!=$rt ){ die "Wrong number of individuals in file $filename!\n" }
		$totalind=$rt unless defined $totalind;
		$self->{_TotalSegs}+=$length;		
		}
	$self->{_TotalInd}=$totalind;
	return;
	}


sub printNexus{
	my ($self, $outfile, $ref_poplabels)=@_;

	open my $output, ">", $outfile or die "Could not open nexus file $outfile!\n";

	print $output "#NEXUS\n\nBegin data;\n\tDimensions ntax=$self->{_TotalInd} nchar=$self->{_TotalSegs};\n\tFormat datatype=STANDARD symbols=\"01\" gap=- missing=?;\n\tMatrix\n";
	for my $pop ( 0..$#{$self->{_Haplotypes}} ){
		my $poplabel=defined $ref_poplabels ? $ref_poplabels->[$pop] : "pop" . $pop;	# use population abbreviations if available, otherwise population index.
		my $ind=0;
		for my $haplo ( @{ $self->{_Haplotypes}[$pop] } ){ print $output $poplabel, "_ind", $ind++, "\t$haplo\n" }
		}
	print $output "\t;\nEnd;\n";
	close $output;

	return;
	}

sub printPLINK{
	my ($self, $outfile)=@_;
	open my $output, ">", $outfile or die "Could not open output file $outfile!\n";
	for my $ind (0..$#{ $self->{_PLINK} }){
		if($self->{_PLINK}[$ind]=~/[^ATGC0 ]/){ warn "Illegal character detected for individual $ind!\n" }
		unless (length($self->{_PLINK}[$ind])/4==$self->{_TotalSNPs}){ warn "Wrong number of SNPs for individual $ind! $self->{_TotalSNPs} vs. ", length($self->{_PLINK}[$ind])/4, ".\n" } 
		print $output "ind", $ind, " ind", $ind, " 0 0 0 0", $self->{_PLINK}[$ind], "\n";
		}
	close $output;
	return;
	}



1;


