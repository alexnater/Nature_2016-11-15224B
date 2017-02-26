package Trees;

use strict;
use warnings;
use Classes::Misc;
use Bio::TreeIO;
use IO::String;


sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	return $self;
	}


sub readTrees{
	my ($self, $treefile, $format, $maxtrees)=@_;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	if ($treefile eq "none"){ return }
	unless ($format eq 'newick' || $format eq 'bed'){
		warn "Format $format not recognized!\n";
		return;
		}
	my @trees;

	open my $input, "<", $treefile or die "Could not open tree file $treefile!\n$!\n";
	my $ntrees=0;
	while (<$input>){
		if ($format eq 'newick' && m/^\(.*;/){
			push @trees, { 'tree'=>$1 };
			++$ntrees;
			}
		elsif ($format eq 'bed' && m/^(\w+)\s(\d+)\s(\d+)\s(\(.*;)/){
			push @trees, { 'id'=>$1, 'start'=>$2, 'end'=>($3-1), 'tree'=>$4 };
			++$ntrees;
			}
		last if ($ntrees>=$maxtrees);
		}
	print STDERR "Read a total of $ntrees trees.\n";
	close $input;

	$self->{_Trees}=\@trees;
	return;
	}


sub genotypestoIndhaps{	# Under construction! Are non-biallelic loci an issue here?
	my ($self, $outfolder, $out_prefix, $ref_loci, $refseq, $ref_genotypes, $ref_missing, $ref_groups, $ref_groupnames, $ref_outgroups, $maxnotbiallelic, $mincov, $ref_minind, $excluderefN)=@_;

	$self->{_Species}=$ref_groupnames;
	$self->{_Outgroups}=$ref_outgroups;
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @groupsizes, @popsum);
	my $totsize=0;
	for my $ref ( @$ref_groups ){
		my $groupsize=0;
		for my $subpop ( @$ref ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($subpop);
			push @popsizes, $temp;
			$groupsize+=$temp;
			$totsize+=$temp;
			print "$temp, $totsize\n";
			}
		push @groupsizes, $groupsize;
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my @seq;

			my $iterator;
			if ( defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus) ){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				}
			else {
				warn "No subsampling information found, using all haplotypes.\n";
				my @indlist;
				for my $popsize (@popsizes){ push @indlist, [ 0..$popsize-1 ] }
				$iterator=\@indlist;
				}

			my $notbiallelic=0;

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
					for my $group (0..$#groupsizes){
						$seq[$group][$_].="?" foreach (0..$groupsizes[$group]-1);
						}
					next POSITION;
					}
				if ($snp){	# maybe more efficient to run through all positions first and check if too many not biallelic sites?
					my $nalleles=$ref_genotypes->NAlleles($chrom, $pos);
					if ($nalleles>2){
						++$notbiallelic;
						if ($notbiallelic>$maxnotbiallelic){
							$ref_loci->setExcluded($chrom, $locus, 1);
							print "Locus $chrom - $locus is not biallelic at position $pos. ${notbiallelic}th occurrence, locus excluded!\n";
							next LOCUS;
							}
						else {
							my $stringref=$ref_genotypes->getValue($chrom, $pos);
							for my $pop (0..$#popsizes){
								$stringref->[$pop]='.' x length($stringref->[$pop]);	# set genotype data to missing.
								}
							for my $group (0..$#groupsizes){
								$seq[$group][$_].="?" foreach (0..$groupsizes[$group]-1);
								}
							print "Locus $chrom - $locus is not biallelic at position $pos. ${notbiallelic}th occurrence, so still included.\n";
							next POSITION;
							}
						}
					}

				my %alleles;

				GROUP: for my $group (0..$#{$ref_groups}){
					my $ind=0;
					POP: for my $subpop ( @{$ref_groups->[$group]} ){
						my @bases;
						my $validal=0;
						IND: for my $index ( @{$iterator->[$subpop]} ){
							my $allind=int($index/2)+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){	# missing data, replace reference base at offset with N.
								push @bases, "?";
								next IND;
								}
							if ($snp){	# if position is covered and variable, get genotypes.
								my $allele=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
								if (defined $allele && $allele ne "N"){
									push @bases, $allele;
									$alleles{$allele}=1;
									++$validal;
									}
								else { push @bases, "?"; }
								}
							else {	# if position not variable but covered, add reference base.
								push @bases, $refbase;
								++$validal;
								}
							}
						unless (@bases==$popsizes[$subpop]){ die "Wrong number of bases for population $subpop in array at position $chrom:$pos! ", scalar(@bases), " vs. ", $popsizes[$subpop], "\n" }
						if ($validal<$ref_minind->[$subpop]){	# not enough alleles callable in a population, set base to N. Better leave these columns out to avoid lengthy warning messges in RAxML.
							$seq[$group][$ind++].="?" foreach (@bases);
							}
						else {
							$seq[$group][$ind++].=$_ foreach (@bases);
							}
						}
					}

				if (scalar(keys %alleles)>2){	# check if position is not biallelic.
					for my $group (0..$#groupsizes){
						warn "$chrom:$pos is not biallelic, setting all bases to missing!\n";
						for my $ind (0..$groupsizes[$group]-1){ substr($seq[$group][$ind], -1, 1, "?") }	# replace last base of haplotype string with missing.
						}
					}
				}

			for my $group (0..$#groupsizes){	# check correct length of haplotype strings.
				for my $ind (0..$groupsizes[$group]-1){
					unless ( length($seq[$group][$ind])==$locuslength ){
						die "Locus $chrom:$start-$end has unequal length of haplotype strings and locus length! ", length($seq[$group][$ind]), " vs. $locuslength\n";
						}
					}
				}

			my $ref_treestring=$self->callRAxML( $outfolder, "TEST", \@seq, int(rand(100000)) );
			$self->{_Haplotypes}{$chrom}[$locus]=\@seq;
			}
		}
	return;
	}


sub treeConcordancePhased{
	my ($self, $outfolder, $ref_loci, $refseq, $ref_genotypes, $ref_missing, $ref_minind, $ref_poplabels, $ref_groups, $ref_groupnames, $ref_outgroups_phase, $ref_outgroups_tree, $niter, $mincov, $excluderefN, $minpropvalid, $minpropcomplete, $subsample, $verbose)=@_;

	my $fastPhaseLimit=500000;	# maximum number of polymorphic positions for fastPhase.
	my @ingroups=@$ref_groups;
	my (@outpops_phase, @outpops_tree);
	if (defined $ref_outgroups_phase){
		splice( @ingroups, -(scalar @$ref_outgroups_phase) );	# remove outgroups from list of groups.
		@outpops_phase=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_phase ];
		}
	if (defined $ref_outgroups_tree){	
		@outpops_tree=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_tree ];
		}
	$self->{_Species}=$ref_groupnames;
	$self->{_Outgroups}=\@outpops_phase;

	my $npops=$ref_missing->getNPop();
	my (@popsizes, @groupsizes, @popsum);
	my $totsize=0;
	for my $ref ( @$ref_groups ){
		my $groupsize=0;
		for my $subpop ( @$ref ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($subpop);
			push @popsizes, $temp;
			$groupsize+=$temp;
			$totsize+=$temp;
			print "$temp, $totsize\n";
			}
		push @groupsizes, $groupsize;
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my $teststring=$$ref_seqstring;
			my $hardmasked=($teststring=~tr/Nn//);
			print "Processing $chrom:$start-$end. Found a total of $hardmasked hardmasked sites.\n";
			$ref_loci->{_Loci}{$chrom}[$locus]{'covered'}=$locuslength-$hardmasked;	# only count the sites as covered that are not hardmasked.
			if ($locuslength-$hardmasked < $locuslength*$minpropvalid){
				print "$chrom:$start-$end has to many hardmasked sites! Skipping locus.\n";
				next LOCUS;
				}

			my (@seqs, @haplotypes, @positions);

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				if ($snp && $ref_genotypes->NAlleles($chrom, $pos)>2){	# fastPHASE does not accept non-biallelic SNPs.
					warn "Position $chrom:", $pos+1, " is not bi-allelic, skipping position!\n" if $verbose;
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my (@states, %alleles);

				POP: for my $pop (0..$#popsizes){
					my $validal=0;
					IND: for my $ind (0..$popsizes[$pop]-1){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						unless (defined $coverage && $coverage>=$mincov){	# missing data, replace reference base at offset with N.
							push @{$states[$pop]}, ("?", "?");
							next IND;
							}
						if ($snp){	# if position is covered and variable but covered, get genotypes.
							for my $index (2*$ind, 2*$ind+1){
								my $allele=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
								if (defined $allele && $allele ne "N"){
									push @{$states[$pop]}, $allele;
									$alleles{$allele}=1;
									++$validal;
									}
								else { push @{$states[$pop]}, "?" }
								}
							}
						else {	# if position not variable but covered, add reference base to string.
							push @{$states[$pop]}, ($refbase, $refbase);
							$validal+=2;
							}
						}
					if ($validal<2*$ref_minind->[$pop]){	# not enough alleles callable in a population.
						@{ $states[$pop] }=("?") x (2*$popsizes[$pop]);
						}
					}

				for my $pop (0..$#states){
					unless (@{$states[$pop]}==2*$popsizes[$pop]){
						die "Wrong number of allelic states for population $pop in array at position $chrom:$pos! ", scalar(@{$states[$pop]}), " vs. ", 2*$popsizes[$pop], "\n";
						}
					my $index=0;
					$seqs[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add bases to sequence string.
					}
				if (keys %alleles>1){	# check if position is really polymorphic.
					push @positions, $offset;
					for my $pop (0..$#states){
						my $index=0;
						$haplotypes[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add polymorphic position to haplotype string.
						}
					}
				}

			unless (@positions){ warn "Locus $chrom:", $start+1, "-", $end+1, "has no polymorphic positions, skip locus!\n"; next LOCUS }
			if (@positions > $fastPhaseLimit){ warn "Locus $chrom:$start-$end has to many polymorphic positions (", scalar(@positions), " vs. $fastPhaseLimit)!\n"; next LOCUS }

			for my $pop (0..$#popsizes){	# check correct length of sequence and haplotype strings.
				for my $ind (0..2*$popsizes[$pop]-1){
					unless ( length($seqs[$pop][$ind])==$locuslength ){
						die "Locus $chrom:$start-$end has unequal length of sequence strings and locus length! ", length($seqs[$pop][$ind]), " vs. $locuslength\n";
						}
					unless ( length($haplotypes[$pop][$ind])==@positions ){
						die "Locus $chrom:$start-$end has unequal length of haplotype and position arrays! ", length($haplotypes[$pop][$ind]), " vs. ", scalar(@positions), "\n";
						}
					}
				}

			print "Submitting ingroups to PHASE ...\n";
			my $inputfile="phasing_" . $chrom . "_" . ($start+1) . "_" . ($end+1) . ".inp";
			my $ref_inphased=$self->callPHASE($outfolder, $inputfile, "subpoplabels.txt", "results_phasing", \@haplotypes, \@ingroups, \@positions, 10, 80, 0);

			my @phased=@$ref_inphased;
			if (@outpops_phase){
				print "Submitting outgroups to random phasing ...\n";
				my @outhaps=@haplotypes[ @outpops_phase ];
				my $ref_outphased=$self->getRandomHap(\@outhaps);
				push @phased, @$ref_outphased;
				}

			print "Preparing phased sequence strings ...\n";
			for my $offset (@positions){	# replace polymorphic positions in seq strings with correctly phased genotypes.
				for my $pop (0..$#popsizes){
					for my $ind (0..2*$popsizes[$pop]-1){
						substr( $seqs[$pop][$ind], $offset, 1, substr($phased[$pop][$ind], 0, 1, "") );	# chew up phased string and replace genotype in sequence string.
						}
					}
				}

			print "Subsampling phased sequence strings ...\n";
			my $iterator;
			if ( $subsample && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus) ){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				}
			else {
				warn "No subsampling selected, using all haplotypes.\n";
				my @indlist;
				for my $pop (0..$#popsizes){ push @indlist, [ 0..$#{$seqs[$pop]} ] }
				$iterator=\@indlist;
				}

			my @subseqs_by_pop=map { [ @{ $seqs[$_] }[ @{$iterator->[$_]} ] ] } (0..$#seqs);
			my @subseqs_by_group;
			for my $group (0..$#groupsizes){
				for my $subpop ( @{$ref_groups->[$group]} ){
					push @{$subseqs_by_group[$group]}, \( @{ $seqs[$subpop] }[ @{$iterator->[$subpop]} ] );	# push only references to the sequence string to avoid excessive copying?
					}
				}

			my @outseqs;
			for my $pop (@outpops_tree){
				push @outseqs, $ref_poplabels->[$pop] . "_ind" . $_ for ( 0..$#{ $subseqs_by_pop[$pop] } );
				}

			my @results4species=(0) x 21;

			print "Calling RAxML on complete gene tree ...\n";
			my $infile="haplotypes_". $chrom . ":" . $start . "-" . $end . "_complete.fasta";
			$self->printFasta($outfolder . "/" . $infile, \@subseqs_by_pop, $ref_poplabels, 100);
			my $ref_treestring=$self->callRAxML( $outfolder, $infile, "TEST", \@outseqs, int(rand(100000)) );
#			unlink $outfolder . "/" . $infile;	# delete fasta file.
			my ($ref_dist, $ref_gsi)=$self->calculateTreeStats($ref_treestring, $ref_poplabels, $ref_groups, $ref_groupnames, $ref_outgroups_phase);
			$ref_loci->setValue($chrom, $locus, 'tree', $$ref_treestring);
			$ref_loci->setValue($chrom, $locus, 'gsi', $ref_gsi);
			$ref_loci->setValue($chrom, $locus, 'dist', $ref_dist);

			my $iter=1;
			my $maxiter=4*$niter;
			REP: while ($iter <= $niter && $maxiter--){
				print "Working on iteration $iter of $niter ...\n";
				my @haps=map { ${ $_->[ int( rand( scalar(@$_) ) ) ] } } @subseqs_by_group;	# haps array now contains one random haplotype per species.
				my @redhaps;
				my $emptycols=0;
				my $completecols=0;

				OFFSET: for my $offset (0..$locuslength-1){
					my @bases;
					my $data=0;
					my $complete=1;
					SPECIES: for my $hap (@haps){
						my $base=substr($hap, 0, 1, "");
						if ($base eq "?"){ $complete=0 }
						else { $data=1 }
						push @bases, $base;
						}
					if ($data){
						$redhaps[$_][0].=$bases[$_] foreach (0..$#bases);
						}
					else { ++$emptycols }
					++$completecols if ($complete);
					}
				print "Locus length: $locuslength, complete columns: $completecols, empty columns: $emptycols\n" if $verbose;
				for my $ref_hap (@redhaps){	# check correct length of haplotype strings.
					unless ( length($ref_hap->[0])==$locuslength-$emptycols ){
						die "Locus $chrom:$start-$end has unequal length of haplotype strings and locus length! ", length($ref_hap->[0]), " vs. ", $locuslength-$emptycols, "\n";
						}
					unless ($locuslength-$emptycols >= $locuslength*$minpropvalid){
						print "$chrom:$start-$end in iteration $iter has only ", $locuslength-$emptycols, " sites with data, minimum of ", $locuslength*$minpropvalid, " required! Skipping iteration. $maxiter tries left.\n";
						next REP;
						}
					unless ($completecols >= $locuslength*$minpropcomplete){
						print "$chrom:$start-$end in iteration $iter has only $completecols sites with complete data, minimum of ", $locuslength*$minpropcomplete, " required! Skipping iteration. $maxiter tries left.\n";
						next REP;
						}
					}

				my @outseqs;
				for my $group ( @{ $ref_outgroups_tree } ){
					push @outseqs, $ref_groupnames->[$group] . "_ind" . $_ for ( 0..$#{ $redhaps[$group] } );
					}

				print "Calling RAxML on randomly subsampled gene tree ...\n";
				my $infile="haplotypes_". $chrom . ":" . $start . "-" . $end . "_iter" . $iter . ".fasta";
				$self->printFasta($outfolder . "/" . $infile, \@redhaps, $ref_groupnames, 100);
				my $ref_treestring=$self->callRAxML( $outfolder, $infile, "TEST", \@outseqs, int(rand(100000)) );
				unlink $outfolder . "/" . $infile;	# delete fasta file.
				$self->processTrees6Ind($ref_treestring, \@results4species, $ref_groupnames);
				++$iter;
				}
			if ($results4species[16]>0){
				$results4species[17]/=$results4species[16];
				$results4species[18]/=$results4species[16];
				$results4species[19]/=$results4species[16];
				$results4species[20]=$results4species[18]/$results4species[19];
				}
			else {
				$results4species[17]="NA";
				$results4species[18]="NA";
				$results4species[19]="NA";
				$results4species[20]="NA";
				}
			$ref_loci->setValue($chrom, $locus, 'concordance', \@results4species);
			}
		}
	return;
	}

sub treeConcordancePhasedBootstrap{
	my ($self, $outfolder, $ref_loci, $refseq, $ref_genotypes, $ref_missing, $ref_minind, $ref_poplabels, $ref_groups, $ref_groupnames, $ref_outgroups_phase, $ref_outgroups_tree, $niter, $mincov, $excluderefN, $minpropvalid, $minpropcomplete, $subsample, $bootstrapreps, $verbose)=@_;

	my $fastPhaseLimit=500000;	# maximum number of polymorphic positions for fastPhase.
	my @ingroups=@$ref_groups;
	my (@outpops_phase, @outpops_tree);
	if (defined $ref_outgroups_phase){
		splice( @ingroups, -(scalar @$ref_outgroups_phase) );	# remove outgroups from list of groups.
		@outpops_phase=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_phase ];
		}
	if (defined $ref_outgroups_tree){	
		@outpops_tree=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_tree ];
		}
	$self->{_Species}=$ref_groupnames;
	$self->{_Outgroups}=\@outpops_phase;

	my $npops=$ref_missing->getNPop();
	my (@popsizes, @groupsizes, @popsum);
	my $totsize=0;
	for my $ref ( @$ref_groups ){
		my $groupsize=0;
		for my $subpop ( @$ref ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($subpop);
			push @popsizes, $temp;
			$groupsize+=$temp;
			$totsize+=$temp;
			print "$temp, $totsize\n";
			}
		push @groupsizes, $groupsize;
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my $teststring=$$ref_seqstring;
			my $hardmasked=($teststring=~tr/Nn//);
			print STDERR "Processing $chrom:$start-$end. Found a total of $hardmasked hardmasked sites.\n";
			$ref_loci->{_Loci}{$chrom}[$locus]{'covered'}=$locuslength-$hardmasked;	# only count the sites as covered that are not hardmasked.
			if ($locuslength-$hardmasked < $locuslength*$minpropvalid){
				warn "$chrom:$start-$end has to many hardmasked sites! Skipping locus.\n";
				next LOCUS;
				}

			my (@seqs, @haplotypes, @positions);

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				if ($snp && $ref_genotypes->NAlleles($chrom, $pos)>2){	# fastPHASE does not accept non-biallelic SNPs.
					warn "Position $chrom:", $pos+1, " is not bi-allelic, skipping position!\n" if $verbose;
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my (@states, %alleles);
				if ($snp){
					unless ($ref_genotypes->getBase($chrom, $pos, 0) eq $refbase){ die "Reference bases do not match at $chrom:$pos!\n" }
					}

				POP: for my $pop (0..$#popsizes){
					my $validal=0;
					IND: for my $ind (0..$popsizes[$pop]-1){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						unless (defined $coverage && $coverage>=$mincov){	# missing data, replace reference base at offset with N.
							push @{$states[$pop]}, ("?", "?");
							next IND;
							}
						if ($snp){	# if position is covered and variable, get genotypes.
							for my $index (2*$ind, 2*$ind+1){
								my $allele=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
								if (defined $allele && $allele ne "N"){
									push @{$states[$pop]}, $allele;
									$alleles{$allele}=1;
									++$validal;
									}
								else { push @{$states[$pop]}, "?" }
								}
							}
						else {	# if position not variable but covered, add reference base to string.
							push @{$states[$pop]}, ($refbase, $refbase);
							$validal+=2;
							}
						}
					if ($validal<2*$ref_minind->[$pop]){	# not enough alleles callable in a population.
						@{ $states[$pop] }=("?") x (2*$popsizes[$pop]);
						}
					}

				for my $pop (0..$#states){
					unless (@{$states[$pop]}==2*$popsizes[$pop]){
						die "Wrong number of allelic states for population $pop in array at position $chrom:$pos! ", scalar(@{$states[$pop]}), " vs. ", 2*$popsizes[$pop], "\n";
						}
					my $index=0;
					$seqs[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add bases to sequence string.
					}
				if (keys %alleles>1){	# check if position is really polymorphic.
					push @positions, $offset;
					for my $pop (0..$#states){
						my $index=0;
						$haplotypes[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add polymorphic position to haplotype string.
						}
					}
				}

			unless (@positions){ warn "Locus $chrom:", $start+1, "-", $end+1, "has no polymorphic positions, skip locus!\n"; next LOCUS }
			if (@positions > $fastPhaseLimit){ warn "Locus $chrom:$start-$end has to many polymorphic positions (", scalar(@positions), " vs. $fastPhaseLimit)!\n"; next LOCUS }

			for my $pop (0..$#popsizes){	# check correct length of sequence and haplotype strings.
				for my $ind (0..2*$popsizes[$pop]-1){
					unless ( length($seqs[$pop][$ind])==$locuslength ){
						die "Locus $chrom:$start-$end has unequal length of sequence strings and locus length! ", length($seqs[$pop][$ind]), " vs. $locuslength\n";
						}
					unless ( length($haplotypes[$pop][$ind])==@positions ){
						die "Locus $chrom:$start-$end has unequal length of haplotype and position arrays! ", length($haplotypes[$pop][$ind]), " vs. ", scalar(@positions), "\n";
						}
					}
				}

			print STDERR "Submitting ingroups to PHASE ...\n";
			my $inputfile="phasing_" . $chrom . "_" . ($start+1) . "_" . ($end+1) . ".inp";
			my $ref_inphased=$self->callPHASE($outfolder, $inputfile, "subpoplabels.txt", "results_phasing", \@haplotypes, \@ingroups, \@positions, 10, 80, 0);

			my @phased=@$ref_inphased;
			if (@outpops_phase){
				print STDERR "Submitting outgroups to random phasing ...\n";
				my @outhaps=@haplotypes[ @outpops_phase ];
				my $ref_outphased=$self->getRandomHap(\@outhaps);
				push @phased, @$ref_outphased;
				}

			print STDERR "Preparing phased sequence strings ...\n";
			for my $offset (@positions){	# replace polymorphic positions in seq strings with correctly phased genotypes.
				for my $pop (0..$#popsizes){
					for my $ind (0..2*$popsizes[$pop]-1){
						substr( $seqs[$pop][$ind], $offset, 1, substr($phased[$pop][$ind], 0, 1, "") );	# chew up phased string and replace genotype in sequence string.
						}
					}
				}

			print STDERR "Subsampling phased sequence strings ...\n";
			my $iterator;
			if ( $subsample && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus) ){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				}
			else {
				warn "No subsampling selected, using all haplotypes.\n";
				my @indlist;
				for my $pop (0..$#popsizes){ push @indlist, [ 0..$#{$seqs[$pop]} ] }
				$iterator=\@indlist;
				}

			my @subseqs_by_pop=map { [ @{ $seqs[$_] }[ @{$iterator->[$_]} ] ] } (0..$#seqs);
			my @seqnames_by_group;
			for my $group (0..$#groupsizes){
				for my $subpop ( @{$ref_groups->[$group]} ){
					push @{ $seqnames_by_group[$group] }, map { $ref_poplabels->[$subpop] . "_ind" . $_ } ( 0..$#{ $subseqs_by_pop[$subpop] } );
					}
				}

			my @outseqs;
			for my $pop (@outpops_tree){
				push @outseqs, $ref_poplabels->[$pop] . "_ind" . $_ for ( 0..$#{ $subseqs_by_pop[$pop] } );
				}

			print STDERR "Calling RAxML on complete gene tree ...\n";
			my $infile="haplotypes_". $chrom . ":" . $start . "-" . $end . "_complete.fasta";
			$self->printFasta($outfolder . "/" . $infile, \@subseqs_by_pop, $ref_poplabels, 100);
			my ($ref_treestring, $treefilename, $bfilename)=$self->callRAxML( $outfolder, $infile, "TEST", \@outseqs, int(rand(100000)), $bootstrapreps );
#			unlink $outfolder . "/" . $infile;	# delete fasta file.
			my ($ref_dist, $ref_gsi)=$self->calculateTreeStats($ref_treestring, $ref_poplabels, $ref_groups, $ref_groupnames, $ref_outgroups_phase);
			$ref_loci->setValue($chrom, $locus, 'tree', $$ref_treestring);
			$ref_loci->setValue($chrom, $locus, 'gsi', $ref_gsi);
			$ref_loci->setValue($chrom, $locus, 'dist', $ref_dist);

			my @results4species=(0) x 52;
			my @bsthresholds=(50,70,85);
			@bsthresholds=map { $_/100*$bootstrapreps } @bsthresholds;

			REP: for my $iter (1..$niter){
				my @keepnodeids=map { $_->[ int( rand( scalar(@$_) ) ) ] } @seqnames_by_group;
				print STDERR "Iteration $iter: ", join(",", @keepnodeids), "\n" if $verbose;
				system("nw_prune -v ${outfolder}/${bfilename} @keepnodeids 1>${outfolder}/bootstrap_subsampled.txt 2>/dev/null");
				my $btreestring=`nw_prune -v ${outfolder}/${treefilename} @keepnodeids | nw_support - bootstrap_subsampled.txt 2>/dev/null`;
				if (@keepnodeids==5){ $self->processTrees5IndBootstrap(\$btreestring, \@results4species, \@keepnodeids, \@bsthresholds) }
				else { $self->processTrees6IndBootstrap(\$btreestring, \@results4species, \@keepnodeids, \@bsthresholds) }
				}
			my $valit=$results4species[16]-$results4species[15];
			if ($valit>0){
				$results4species[32]/=$valit;
				$results4species[33]/=$valit;
				$results4species[$_]/=$valit foreach (17..31);
				}
			else {
				$results4species[32]="NA";
				$results4species[33]="NA";
				}
			$ref_loci->setValue($chrom, $locus, 'concordance', \@results4species);
			}
		}
	return;
	}


sub treeConcordance{
	my ($self, $outfolder, $ref_loci, $refseq, $ref_genotypes, $ref_missing, $ref_groups, $ref_groupnames, $ref_comparisons, $ref_outgroups, $niter, $mincov, $excluderefN, $minpropvalid, $minpropcomplete, $subsample)=@_;

	$self->{_Species}=$ref_groupnames;
	$self->{_Outgroups}=$ref_outgroups;
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @groupsizes, @popsum);
	my $totsize=0;
	for my $ref ( @$ref_groups ){
		my $groupsize=0;
		for my $subpop ( @$ref ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($subpop);
			push @popsizes, $temp;
			$groupsize+=$temp;
			$totsize+=$temp;
			print "$temp, $totsize\n";
			}
		push @groupsizes, $groupsize;
		}

#	$ref_loci->{_Trees}{'comparisons'}=$ref_comparisons;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my $teststring=$$ref_seqstring;
			my $hardmasked=($teststring=~tr/Nn//);
			$ref_loci->{_Loci}{$chrom}[$locus]{'covered'}=$locuslength-$hardmasked;	# only count the sites as covered that are not hardmasked.
			if ($locuslength-$hardmasked < $locuslength*$minpropvalid){
				print "$chrom:$start-$end has $hardmasked hardmasked sites! Skipping locus.\n";
				next LOCUS;
				}

			my $iterator;
			if ( $subsample && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus) ){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				}
			else {
				warn "No subsampling information found, using all haplotypes.\n";
				my @indlist;
				for my $popsize (@popsizes){ push @indlist, [ 0..2*$popsize-1 ] }
				$iterator=\@indlist;
				}
			my @indices;	# array containing the subpop id and individual number of all individuals for each group.
			for my $group (0..$#groupsizes){
				for my $subpop ( @{$ref_groups->[$group]} ){
					my %indid=map { int($_/2) => 1 } @{$iterator->[$subpop]};	# determine individual id from index and remove duplicates if not haplotized.
					push @{$indices[$group]}, [ $subpop, $_ ] foreach sort {$a<=>$b} keys %indid;
					}
				}

			my @results=(0) x (scalar(@$ref_comparisons)+1);
			my @results4species=(0) x 21;
			my $iter=1;
			my $maxiter=4*$niter;
			REP: while ($iter <= $niter && $maxiter--){
				print "Working on iteration $iter of $niter ...\n";
				my $seqstring=$$ref_seqstring;	# make hard copy of reference string.
				my @seq;
				my @inds=map { $_->[ int( rand( scalar(@$_) ) ) ] } @indices;	# inds array now contains one random subpop/indid pair for each species.
				my $emptycols=0;
				my $completecols=0;

				POSITION: for my $pos ($start..$end){
					my $offset=$pos-$start;
					my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
					my $refbase=uc( substr($seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
					if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
						++$emptycols;
						next POSITION;
						}
					my $data=0;
					my $complete=1;
					my @bases=("?") x scalar(@groupsizes);
	
					GROUP: for my $group (0..$#groupsizes){
						my $subpop=$inds[$group][0];
						my $ind=$inds[$group][1];
						my $allind=$ind + $popsum[$subpop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						unless (defined $coverage && $coverage>=$mincov){ $complete=0; next GROUP }
						if ($snp){	# if position is covered and variable, get genotypes.
							my $index=2*$ind + int( rand(2) );	# randomly select one of the two chromosomes.
							my $allele=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
							if (defined $allele && $allele ne "N"){
								$bases[$group]=$allele;
								$data=1;
								}
							}
						else {	# if position not variable but covered, add reference base.
							$bases[$group]=$refbase;
							$data=1;
							}
						}
					if ($data){
						$seq[$_][0].=$bases[$_] foreach (0..$#bases);
						}
					else { ++$emptycols }
					++$completecols if ($complete);
					}

				for my $group (0..$#groupsizes){	# check correct length of haplotype strings.
					unless ( length($seq[$group][0])==$locuslength-$emptycols ){
						die "Locus $chrom:$start-$end has unequal length of haplotype strings and locus length! ", length($seq[$group][0]), " vs. ", $locuslength-$emptycols, "\n";
						}
					unless ($locuslength-$emptycols >= $locuslength*$minpropvalid){
						print "$chrom:$start-$end in iteration $iter has only ", $locuslength-$emptycols, " sites with data, minimum of ", $locuslength*$minpropvalid, " required! Skipping iteration. $maxiter tries left.\n";
						next REP;
						}
					unless ($completecols >= $locuslength*$minpropcomplete){
						print "$chrom:$start-$end in iteration $iter has only $completecols sites with complete data, minimum of ", $locuslength*$minpropcomplete, " required! Skipping iteration. $maxiter tries left.\n";
						next REP;
						}
					}

				my @outseqs;
				for my $group ( @{ $ref_outgroups } ){
					push @outseqs, $ref_groupnames->[$group] . "_ind" . $_ for ( 0..$#{ $seq[$group] } );
					}

				print "Calling RAxML ...\n";
				my $infile="haplotypes_". $chrom . ":" . $start . "-" . $end . "_iter" . $iter . ".fasta";
				$self->printFasta($outfolder . "/" . $infile, \@seq, $ref_groupnames, 100);
				my $ref_treestring=$self->callRAxML( $outfolder, $infile, "TEST", \@outseqs, int(rand(100000)) );
				unlink $outfolder . "/" . $infile;
				$self->processTrees6Ind($ref_treestring, \@results4species, $ref_groupnames);
				++$iter;
				}
			$results4species[17]/=$results4species[16];
			$results4species[18]/=$results4species[16];
			$results4species[19]/=$results4species[16];
			$results4species[20]=$results4species[18]/$results4species[19];
			$ref_loci->{_Loci}{$chrom}[$locus]{'concordance'}=\@results4species;
			}
		}
	return;
	}


sub genotypestoFastPhase{
	my ($self, $outfolder, $out_prefix, $ref_loci, $refseq, $ref_genotypes, $ref_popsizes, $ref_groups, $ref_poplabels, $ref_minind, $excluderefN, $snplimit, $overlap, $allelemode, $verbose)=@_;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my (@haplotypes, @positions);
			my ($windowstart, $windowend)=($start, $end);

			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless( defined $ref_genotypes->getValue($chrom, $pos) );
				next POSITION if ($excluderefN && $refbase eq 'N');	# if reference has a N at the current position, skip position.
				if ($ref_genotypes->NAlleles($chrom, $pos)>2){	# fastPHASE does not accept non-biallelic SNPs.
					warn "Position $chrom:", $pos+1, " is not bi-allelic, skipping position!\n" if $verbose;
					next POSITION;
					}
				my $ancestral;
				if ($allelemode eq 'ancder'){
					$ancestral=$ref_genotypes->getAncestral($chrom, $pos);
					unless (defined $ancestral && $ancestral<4){
						warn "Ancestral state not defined for $chrom:", $pos+1, ", skipping position!\n" if $verbose;
						next POSITION;
						} 
					}
				my (@states, %alleles);

				GROUP: for my $group ( 0..$#{ $ref_groups } ){
					POP: for my $pop ( @{ $ref_groups->[$group] } ){
						my $validal=0;
						my $gtstring=$ref_genotypes->getValue($chrom, $pos, $pop);
						unless ( $gtstring && length($gtstring)==(2*$ref_popsizes->[$pop]) ){
							warn "Genotype string not defined or wrong length for $chrom:", $pos+1, ", population $pop, setting to missing!\n" if $verbose;
							$gtstring='.' x (2*$ref_popsizes->[$pop]);
							}
						INDEX: for my $index ( 0..2*$ref_popsizes->[$pop]-1 ){
							my $allele=substr($gtstring, 0, 1, "");	# chew up genotype string.
							if (defined $allele && $allele ne "."){
								my $state;
								if ($allelemode eq 'refalt'){ $state=$allele }
								elsif ($allelemode eq 'ancder'){
									if ($ancestral==0){ $state=$allele }
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
					unless (@positions){ warn "Locus $chrom:", $windowstart+1, "-", $pos+1, "has no polymorphic positions, skip locus!\n"; next LOCUS }
					$self->printFastPHASE($outfolder, $out_prefix, \@positions, \@haplotypes, $chrom, $windowstart+1, $pos+1);

					for my $ref_pop (@haplotypes){	# reset haplotype array but keep overlap snps.
						for my $hap (@$ref_pop){
							substr($hap, 0, $snplimit-$overlap, "");
							unless (length($hap)==$overlap){ die "Haplotye array for $chrom:", $windowstart+1, "-", $pos+1, "has not been reset correctly!\n" }
							}						
						}
					@positions=splice(@positions, $snplimit-$overlap);	# reset position array but keep overlap snps.
					unless (@positions==$overlap){ die "Position array for locus $chrom:", $windowstart+1, "-", $pos+1, "has not been reset correctly!\n" }
					$windowstart=@positions ? $positions[0] : $pos+1;	# if there is a window overlap, take position of first snp as windowstart.
					}
				}

			if (@positions){	# print last input file with number of snps < snplimit.
				$self->printFastPHASE($outfolder, $out_prefix, \@positions, \@haplotypes, $chrom, $windowstart+1, $end+1);
				}

			if (@$ref_groups>1){	# if more than one group of populations in the same file, a subpopulation file needs to be printed.
				my @labels;
				for my $group ( 0..$#{ $ref_groups } ){
					push @labels, ($group+1) x $ref_popsizes->[$_] foreach ( @{ $ref_groups->[$group] } );
					}
				print "Printing subpopulation file ...\n";
				my $subpopfile="subpopulations_" . $chrom . "_" . ($start+1) . "_" . ($end+1) . ".txt";
				open my $subpop, ">", $subpopfile or die "Could not open phase input file $subpopfile!\n";
				print $subpop join(' ', @labels);
				close $subpop;
				}
			}
		}
	return;
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


sub genotypestoPhase{	# seperate method for haplotized sequences?
	my ($self, $outfolder, $tempdir, $out_prefix, $ref_loci, $refseq, $chrommap, $ref_genotypes, $ref_missing, $ref_minind, $ref_poplabels, $ref_groups, $ref_outgroups_phase, $ref_outgroups_tree, $mincov, $minpropvalid, $subsample, $excluderefN, $storehaps, $printfasta, $generatetrees, $verbose)=@_;

	my $fastPhaseLimit=500000;	# maximum number of polymorphic positions for fastPhase.
	my @ingroups=@$ref_groups;
	my (@outpops_phase, @outpops_tree);
	if (defined $ref_outgroups_phase){
		splice( @ingroups, -(scalar @$ref_outgroups_phase) );	# remove outgroups from list of groups.
		@outpops_phase=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_phase ];
		}
	if (defined $ref_outgroups_tree){	
		@outpops_tree=map { @$_ } @{ $ref_groups }[ @$ref_outgroups_tree ];
		}

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

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys() ){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if ( $ref_loci->getExcluded($chrom, $locus) );
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my $teststring=$$ref_seqstring;
			my $hardmasked=($teststring=~tr/Nn//);
			print STDERR "Processing $chrom:$start-$end. Found a total of $hardmasked hardmasked sites.\n";
			$ref_loci->{_Loci}{$chrom}[$locus]{'covered'}=$locuslength-$hardmasked;	# only count the sites as covered that are not hardmasked.
			if ($locuslength-$hardmasked < $locuslength*$minpropvalid){
				print STDERR "$chrom:$start-$end has to many hardmasked sites! Skipping locus.\n";
				next LOCUS;
				}

			my (@seqs, @haplotypes, @positions);

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				if ($snp && $ref_genotypes->NAlleles($chrom, $pos)>2){	# fastPHASE does not accept non-biallelic SNPs.
					warn "Position $chrom:", $pos+1, " is not bi-allelic, skipping position!\n" if $verbose;
					for my $pop (0..$#popsizes){ $seqs[$pop][$_].="?" foreach (0..2*$popsizes[$pop]-1) }
					next POSITION;
					}
				my (@states, %alleles);
				if ($snp){
					unless ($ref_genotypes->getBase($chrom, $pos, 0) eq $refbase){ die "Reference bases do not match at $chrom:$pos!\n" }
					}

				POP: for my $pop (0..$#popsizes){
					my $validal=0;
					IND: for my $ind (0..$popsizes[$pop]-1){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						unless (defined $coverage && $coverage>=$mincov){	# missing data, replace reference base at offset with N.
							push @{$states[$pop]}, ("?", "?");
							next IND;
							}
						if ($snp){	# if position is covered and variable, get genotypes.
							for my $index (2*$ind, 2*$ind+1){
								my $allele=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
								if (defined $allele && $allele ne "N"){
									push @{$states[$pop]}, $allele;
									$alleles{$allele}=1;
									++$validal;
									}
								else { push @{$states[$pop]}, "?" }
								}
							}
						else {	# if position not variable but covered, add reference base to string.
							push @{$states[$pop]}, ($refbase, $refbase);
							$validal+=2;
							}
						}
					if ($validal<2*$ref_minind->[$pop]){	# not enough alleles callable in a population.
						@{ $states[$pop] }=("?") x (2*$popsizes[$pop]);
						}
					}

				for my $pop (0..$#states){
					unless (@{$states[$pop]}==2*$popsizes[$pop]){
						die "Wrong number of allelic states for population $pop in array at position $chrom:$pos! ", scalar(@{$states[$pop]}), " vs. ", 2*$popsizes[$pop], "\n";
						}
					my $index=0;
					$seqs[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add bases to sequence string.
					}
				if (keys %alleles>1){	# check if position is really polymorphic.
					push @positions, $offset;
					for my $pop (0..$#states){
						my $index=0;
						$haplotypes[$pop][$index++].=$_ foreach @{ $states[$pop] };	# add polymorphic position to haplotype string.
						}
					}
				}

			unless (@positions){ warn "Locus $chrom:", $start+1, "-", $end+1, "has no polymorphic positions, skip locus!\n"; next LOCUS }
			if (@positions > $fastPhaseLimit){ warn "Locus $chrom:$start-$end has to many polymorphic positions (", scalar(@positions), " vs. $fastPhaseLimit)!\n"; next LOCUS }

			for my $pop (0..$#popsizes){	# check correct length of sequence and haplotype strings.
				for my $ind (0..2*$popsizes[$pop]-1){
					unless ( length($seqs[$pop][$ind])==$locuslength ){
						die "Locus $chrom:$start-$end has unequal length of sequence strings and locus length! ", length($seqs[$pop][$ind]), " vs. $locuslength\n";
						}
					unless ( length($haplotypes[$pop][$ind])==@positions ){
						die "Locus $chrom:$start-$end has unequal length of haplotype and position arrays! ", length($haplotypes[$pop][$ind]), " vs. ", scalar(@positions), "\n";
						}
					}
				}

			print STDERR "Submitting ingroups to PHASE ...\n";
			my $inputfile="phasing_" . $chrom . "_" . ($start+1) . "_" . ($end+1) . ".inp";
			my $ref_inphased=$self->callPHASE($tempdir, $inputfile, "subpoplabels.txt", "results_phasing", \@haplotypes, \@ingroups, \@positions, 10, 80, 0);

			my @phased=@$ref_inphased;
			if (@outpops_phase){
				print STDERR "Submitting outgroups to random phasing ...\n";
				my @outhaps=@haplotypes[ @outpops_phase ];
				my $ref_outphased=$self->getRandomHap(\@outhaps);
				push @phased, @$ref_outphased;
				}

			print STDERR "Preparing phased sequence strings ...\n";
			for my $offset (@positions){	# replace polymorphic positions in seq strings with correctly phased genotypes.
				for my $pop (0..$#popsizes){
					for my $ind (0..2*$popsizes[$pop]-1){
						substr( $seqs[$pop][$ind], $offset, 1, substr($phased[$pop][$ind], 0, 1, "") );	# chew up phased string and replace genotype in sequence string.
						}
					}
				}

			print STDERR "Subsampling phased sequence strings ...\n";
			my $iterator;
			if ( $subsample && defined $ref_loci->getreftoIndividuals('locus', $chrom, $locus) ){
				$iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
				unless ($npops==scalar(@{$iterator})){ die "Wrong information about number of populations!\n" }
				}
			else {
				warn "No subsampling selected, using all haplotypes.\n";
				my @indlist;
				for my $pop (0..$#popsizes){ push @indlist, [ 0..$#{$seqs[$pop]} ] }
				$iterator=\@indlist;
				}
			my @subseqs=map { [ @{ $seqs[$_] }[ @{$iterator->[$_]} ] ] } (0..$#seqs);
			next LOCUS unless ($printfasta || $generatetrees || $storehaps);

			my $infile=$out_prefix . "_". $chrom . ":" . $start . "-" . $end . ".fasta";
			if (defined $chrommap){
				my ($trchrom1, $trstart)=$chrommap->findPosonChrom($chrom, $start);
				my ($trchrom2, $trend)=$chrommap->findPosonChrom($chrom, $end);
				if ($trchrom1 eq $trchrom2){
					my @pos=sort {$a<=>$b} ($trstart, $trend);
					$infile="haplotypes_". $chrom . ":" . $start . "-" . $end . "_" . $trchrom1 . ":" . $pos[0] . "-" . $pos[1] . ".fasta";
					}
				}

			if ($storehaps){
				$self->{_Loci}{$chrom}[$locus]{'haplotypes'}=\@subseqs;
				}

			if ($printfasta){
				print STDERR "Printing Fasta file ...\n";
				$self->printFasta($outfolder . "/" . $infile, \@subseqs, $ref_poplabels, 100);
				}

			if ($generatetrees){
				my @outseqs;
				for my $outpop (@outpops_tree){
					push @outseqs, $ref_poplabels->[$outpop] . "_ind" . $_ for ( 0..$#{ $subseqs[$outpop] } );
					}

				print STDERR "Calling RAxML ...\n";
				$self->printFasta($tempdir . "/" . $infile, \@subseqs, $ref_poplabels, 100);
				my $ref_treestring=$self->callRAxML( $tempdir, $infile, "TEST", \@outseqs, int(rand(100000)) );
				unlink $tempdir . "/" . $infile;
				$ref_loci->setValue($chrom, $locus, 'tree', $$ref_treestring);
				}
			}
		}
	return;
	}


sub printsubsampledTrees{
	my ($self, $outfile, $ref_popsizes, $ref_poplabels, $ref_groups, $ref_groupnames, $niter, $verbose)=@_;

	my @groupsizes=map { Misc::sum(@{ $ref_popsizes }[@$_]) } @$ref_groups;
	my @seqnames_by_group;
	for my $group (0..$#groupsizes){
		for my $subpop ( @{$ref_groups->[$group]} ){
			push @{ $seqnames_by_group[$group] }, map { $ref_poplabels->[$subpop] . "_ind" . $_ } (0..$ref_popsizes->[$subpop]-1);
			}
		}

	my $curtree=0;
	TREE: for my $ref_treehash ( @{ $self->{_Trees} } ){
		print STDERR "Working on tree ", ++$curtree, " ...\n";
		my $treestring=$ref_treehash->{'tree'};
		REP: for my $iter (1..$niter){
			my @keepnodeids=map { $_->[ int( rand( scalar(@$_) ) ) ] } @seqnames_by_group;
			print STDERR "Iteration $iter: ", join(",", @keepnodeids), "\n" if $verbose;
			`echo "$treestring" | nw_prune -v - @keepnodeids 1>>$outfile 2>/dev/null`;
#			my $treestring=`echo $treestring | nw_prune -v - @keepnodeids 2>/dev/null`;
			}
		}

	return;
	}


sub speciesTreeBootstrap{
	my ($self, $outfolder, $workdir, $treefile, $ref_ntaxa, $ref_poplabels, $ref_groups, $ref_groupnames, $ref_outgroups, $niter, $nsamples, $replacement, $verbose)=@_;

	my $subtreefile="subsampled.trees";
	my $ctrlfile="control.txt";
	my @taxa;
	for my $species ( 0..$#{ $ref_groups } ){
		for my $pop ( @{ $ref_groups->[$species] } ){
			my @poptaxa=map { $ref_poplabels->[$pop] . "_ind" . $_ } (0..$ref_ntaxa->[$pop]-1);
#			push @taxa, \@poptaxa;
			push @{ $taxa[$species] }, @poptaxa;
			}
		}
	my @outgroups=map { $ref_groupnames->[$_] } @$ref_outgroups;

	open my $output, ">", "$workdir/$ctrlfile" or die "Could not open control file $workdir/$ctrlfile!\n$!\n";
	print $output "$subtreefile\n0\n-1\n$nsamples ", scalar(@taxa), "\n";
	print $output "$ref_groupnames->[$_] ", scalar(@{ $taxa[$_] }), " ", join(" ", @{ $taxa[$_] }), "\n" foreach (0..$#taxa);
#	print $output "$ref_poplabels->[$_] ", scalar(@{ $taxa[$_] }), " ", join(" ", @{ $taxa[$_] }), "\n" foreach (0..$#taxa);
	print $output "0\n";
	close $output;

	open my $input, "<", $treefile or die "Could not open tree file $treefile!\n$!\n";
	my @trees=<$input>;
	close $input;

	my @results4species=(0) x 21;
	my @treestrings;
	for my $iter (0..$niter-1){
		my $ref_subrows=$replacement ? [ map { int(rand(@trees)) } (1..$nsamples) ] : Misc::samplewithoutRep([0..$#trees], $nsamples);
		open my $output, ">", "$workdir/$subtreefile" or die "Could not open subsampled tree file $workdir/$subtreefile!\n$!\n";
		print $output join('', @trees[@$ref_subrows]), "\n";
		close $output;
		my $ref_treestring=$self->callMPest($workdir, $subtreefile, $ctrlfile, $ref_groupnames, \@outgroups);
#		my $ref_treestring=$self->callMPest($workdir, $subtreefile, $ctrlfile, $ref_poplabels, \@outgroups);
		push @treestrings, $$ref_treestring;
		$self->processTreesBootstrap($ref_treestring, \@results4species, \@taxa, $ref_poplabels, $ref_groups, $ref_groupnames);
		}

	return(\@results4species, \@treestrings);
	}


sub simulatedTrees{
	my ($self, $outfile, $ref_ntaxa, $ref_poplabels, $ref_groups, $ref_groupnames, $niter, $printpertree, $printtreestrings)=@_;

	my %taxa;
	my @taxabyspecies;
	my $taxacount=0;
	for my $species ( 0..$#{ $ref_groups } ){
		my @speciestaxa;
		if (defined $ref_poplabels){
			for my $pop ( @{ $ref_groups->[$species] } ){
				my @poptaxa=map { $ref_poplabels->[$pop] . "_ind" . $_ } (0..$ref_ntaxa->[$pop]-1);
				push @speciestaxa, @poptaxa;
				}
			}
		else {
			@speciestaxa=map { $ref_groupnames->[$species] . "_ind" . $_ } (0..$ref_ntaxa->[$species]-1);
			}
		$taxa{++$taxacount}=$_ foreach @speciestaxa;
		push @taxabyspecies, \@speciestaxa;
		}

#	warn "Wrong number of taxa in hash: ", $taxacount, " vs. ", scalar(@taxabyspecies), "!\n" unless ($taxacount==@taxabyspecies);

	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n$!\n";
	open my $outstream, ">", "${outfile}.trees" or die "Could not open outfile ${outfile}.trees!\n$!\n";
	my $outtree=new Bio::TreeIO(-fh => $outstream, -format => "newick");

	my @allresults=(0) x 21;
	my $curtree=0;
	TREE: for my $ref_treehash ( @{ $self->{_Trees} } ){
		my $treestring=$ref_treehash->{'tree'};
		$treestring=~s/([\(\,])$_:/${1}$taxa{$_}:/ foreach (1..$taxacount);
#		print STDERR "$treestring\n"; 
		my $treeio=IO::String->new($treestring);
		my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");
		my $tree=$input->next_tree;
		my @results4species=(0) x 21;

		ITER: for my $iter (0..$niter-1){
			my $treecopy=$tree->clone();
			warn "Tree is not binary!\n" unless ( $treecopy->is_binary() );
			my @keepnodeids=map { $_->[ int( rand( scalar(@$_) ) ) ] } @taxabyspecies;
			print STDERR "Tree $curtree, iteration $iter: ", join(",", @keepnodeids), "\n";
			my %internal_nodes;
			my %external_nodes;
			for my $extnodeid (@keepnodeids){
				my $extnode=$treecopy->find_node(-id => $extnodeid);
				my @lineagenodes=$treecopy->get_lineage_nodes($extnode);
				$external_nodes{ $extnode->internal_id }=$extnode;
				$internal_nodes{ $_->internal_id }=$_ foreach (@lineagenodes);
#				my @keepids=map { $_->internal_id } (@lineagenodes, $extnode);
#				print STDERR "$extnodeid: ", join(",", @keepids), "\n";
				}
			my @keep_internal_ids=sort {$a<=>$b} (keys %internal_nodes, keys %external_nodes);
#			print STDERR "Tree $curtree, iteration $iter: ", join(",", @keep_internal_ids), "\n";
			$treecopy->splice(-keep_internal_id => \@keep_internal_ids, -preserve_lengths => 1);
			Misc::contract_linear_paths($treecopy, 1);	# original method does not support preservation of branch lengths.
			my @all_nodes=$treecopy->get_nodes;
			warn "Wrong number of nodes!\n" unless (@all_nodes==2*scalar(@keepnodeids)-1);
			$outtree->write_tree($treecopy) if $printtreestrings;
#			print STDERR "Tree $curtree, iteration $iter: ", join(",", map { $_->internal_id } $treecopy->get_nodes), "\n";

			my @ingroups=@keepnodeids;
			my $outgroup1=pop @ingroups;
			my $outgroup2=pop @ingroups;

			my @innodes=grep { defined $_->id && $_->id !~ /$outgroup1/ && $_->id !~ /$outgroup2/ } $treecopy->get_leaf_nodes;
			my $root=$treecopy->get_root_node;
			my $outnode1=$treecopy->find_node(-id => $outgroup1);
			my $outnode2=$treecopy->find_node(-id => $outgroup2);
#			print "Outnode1: ", $outnode1->id, "\nOutnode2: ", $outnode2->id, "\n";

			unless ( $outnode1->ancestor->internal_id == $root->internal_id ){
				warn "Tree $curtree is not correctly rooted:\n$treestring\n";
				++$results4species[15];
				next ITER;
				}

			unless ( $treecopy->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
				warn "Tree $curtree is not monophyletic for ingroups in:\n$treestring\n";
				++$results4species[15];
				next ITER;
				}

			my @speciesnodes;
			for my $species (0..$#ingroups){
				my @nodes=$treecopy->find_node(-id => $ingroups[$species]);
				if (@nodes>1){ warn "Multiple nodes found for $ref_groupnames->[$species], skipping tree!\n"; ++$results4species[15]; next ITER } 
				elsif (@nodes){ push @speciesnodes, @nodes }
				else { warn "No nodes found for $ref_groupnames->[$species], skipping tree!\n"; ++$results4species[15]; next ITER } 
				}

			my (@pair1, @pair2);
			my $pairancestor;
			FIRST: for my $group1 (0..2){
				SECOND: for my $group2 ($group1+1..3){
					if ( $speciesnodes[$group1]->ancestor->internal_id == $speciesnodes[$group2]->ancestor->internal_id ){
						$pairancestor=$speciesnodes[$group1]->ancestor;
						@pair1=($group1, $group2);
						@pair2=grep { $_ != $group1 && $_ != $group2 } (0..3);
						print STDERR "Identified pairings for tree $iter: ", join(',', @pair1), " and ", join(',', @pair2), "\n";
						last FIRST;
						}
					}
				}
			unless (@pair1==2 && @pair2==2){
				warn "No species pairing found for tree $iter:\n$treestring\n";
				++$results4species[15];
				next ITER;
				}

			if ( $speciesnodes[ $pair2[0] ]->ancestor->internal_id == $speciesnodes[ $pair2[1] ]->ancestor->internal_id ){
				print STDERR "Tree $iter is balanced with species pairings ", join(',', @pair1), " and ", join(',', @pair2), "\n";
				my $index=$self->findTreeConf($pair1[0], $pair1[1], 4);
				++$results4species[$index];
				next ITER;
				}
			else{
				print STDERR "Tree $iter is not balanced, with internal species pairing ", join(',', @pair1), ".\n";
				for my $species ($pair2[0], $pair2[1]){
#					print STDERR $pairancestor->ancestor->internal_id, ",", $speciesnodes[$species]->ancestor->internal_id, "\n";
#					print STDERR join(",", map { $_->internal_id } $speciesnodes[$species]->ancestor->each_Descendent), "\n";
#					print STDERR join(",", map { $_->internal_id } $pairancestor->each_Descendent), "\n";
#					print STDERR join(",", map { $_->internal_id } $pairancestor->ancestor->each_Descendent), "\n";
					if ( $pairancestor->ancestor->internal_id == $speciesnodes[$species]->ancestor->internal_id ){
						my $rootpos=( $species==$pair2[0] ) ? $pair2[1] : $pair2[0];
						my $index=$self->findTreeConf($pair1[0], $pair1[1], $rootpos);
						++$results4species[$index];
						next ITER;
						}
					}
				}
			warn "No tree topology found for tree $iter in:\n$treestring\n";
			++$results4species[15];
			} continue { ++$results4species[16] }

		print $output $ref_treehash->{'id'} || 'unmapped', "\t", $ref_treehash->{'start'} || $curtree, "\t", $ref_treehash->{'end'}+1 || $curtree,
						"\t0\t", join("\t", @results4species), "\n" if ($printpertree);
		$allresults[$_]+=$results4species[$_] foreach (0..20);
		} continue { ++$curtree }
	print join("\t", @allresults), "\n";
	return;
	}


sub simulatedTrees6Species{
	my ($self, $outfolder, $ref_groups, $ref_groupnames)=@_;

	my @results4species=(0) x 21;
	my $curtree=0;
	TREE: for my $ref_treehash ( @{ $self->{_Trees} } ){
		print STDERR "Working on tree ", ++$curtree, " ...\n";
		my $treestring=$ref_treehash->{'tree'};
		$treestring=~s/([\(\,])$_:/${1}$ref_groupnames->[$_-1]_ind0:/ foreach (1..6);
		print STDERR "$treestring\n"; 
		$self->processTrees6Ind(\$treestring, \@results4species, $ref_groupnames);
		}

	print join("\t", @results4species), "\n";
	return;
	}


sub callPHASE{
	my ($self, $workdir, $infile, $subpopfile, $out_prefix, $ref_haplotypes, $ref_groups, $ref_positions, $iter, $qualthreshold, $dryrun)=@_;

	my $totsize=0;
	my (@popsizes, @labels);
	for my $group ( 0..$#{ $ref_groups } ){
		for my $pop ( @{ $ref_groups->[$group] } ){
			if (@{ $ref_haplotypes->[$pop] } % 2){ die "Number of haplotypes in population $pop is not divisible by 2!\n" }
			my $size=@{ $ref_haplotypes->[$pop] } / 2;
			push @popsizes, $size;
			$totsize+=$size;
			push @labels, ($group+1) x $size;
			}
		}

	chdir($workdir);

	open my $subpop, ">", $subpopfile or die "Could not open phase input file $subpopfile!\n";
	print $subpop join(' ', @labels);
	close $subpop;

	open my $output, ">", $infile or die "Could not open phase input file $infile!\n";
	print $output "$totsize\n";
	print $output scalar(@$ref_positions), "\n";
	print $output "P ", join(" ", @$ref_positions), "\n";

	for my $group ( 0..$#{ $ref_groups } ){
		for my $pop ( @{ $ref_groups->[$group] } ){
			print $output "$_\n" foreach ( @{$ref_haplotypes->[$pop]} );
			}
		}
	close $output;

	my @args=("fastPHASE", "-n", "-K10", "-q${qualthreshold}", "-T${iter}", "-u${subpopfile}", "-o${out_prefix}", "$infile");
	my $ref_phased;

	unless ($dryrun){
		system(@args);
		if ($? == -1){
			print "failed to execute: $!\n";
			}
		elsif ($? & 127){
			printf "child died with signal %d, %s coredump\n",
			($? & 127), ($? & 128) ? 'with' : 'without';
			}
		else {
			printf "child exited with value %d\n", $? >> 8;
			}
		$ref_phased=$self->readPHASE($workdir, $out_prefix, \@popsizes, 0);
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


sub getRandomHap{
	my ($self, $ref_haplotypes)=@_;

	my @phased;
	for my $pop ( 0..$#{$ref_haplotypes} ){
		if (@{ $ref_haplotypes->[$pop] } % 2){ die "Number of haplotypes in population $pop is not divisible by 2!\n" }
		my $size=@{ $ref_haplotypes->[$pop] } / 2;
		print "Pop: $pop, size: $size\n";
		for (my $ind=0; $ind<$size; ++$ind){
			my ($haplotype1, $haplotype2);
			my @seqs=( $ref_haplotypes->[$pop][2*$ind], $ref_haplotypes->[$pop][2*$ind+1] );	# make hard copies of both unphased haplotypes.
			unless ( length($seqs[0])==length($seqs[1]) ){ die "Unequal length of unphased haplotype strings for population $pop, individual $ind!\n" }
			while (length($seqs[0])>0){
				my @bases=( substr($seqs[0], 0, 1, ""), substr($seqs[1], 0, 1, "") );	# extract first bases in strings and chew them up for faster access.
				if ($bases[0] eq $bases[1]){	# if both bases identical, take the first one for both haplotype strings.
					$haplotype1.=$bases[0];
					$haplotype2.=$bases[0];
					}
				else {
					my $rindex=int( rand(2) );
					$haplotype1.=$bases[ $rindex ];	# randomly select one of the two bases for the first haplotype string.
					$haplotype2.=$bases[ 1-$rindex ];	# take the other base for the second haplotype string.
					}
				}
			unless ( length($haplotype1)==length($ref_haplotypes->[$pop][2*$ind]) ){
				die "Unequal length of phased and unphased haplotype strings for population $pop, individual $ind (", length($haplotype1), " vs. ", length($ref_haplotypes->[$pop][2*$ind]), ")!\n";
				}
			push @{$phased[$pop]}, $haplotype1;

			unless ( length($haplotype2)==length($ref_haplotypes->[$pop][2*$ind]) ){
				die "Unequal length of phased and unphased haplotype strings for population $pop, individual $ind (", length($haplotype2), " vs. ", length($ref_haplotypes->[$pop][2*$ind]), ")!\n";
				}
			push @{$phased[$pop]}, $haplotype2;
			}
		}
	return(\@phased);
	}


sub callRAxML{
	my ($self, $workdir, $infile, $out_postfix, $ref_outgroups, $seed, $nrepbootstrap)=@_;

	chdir($workdir);
	my $bseed=$seed + int(rand(100000));
	my @args=("raxmlHPC-AVX", "-mGTRGAMMA", "-s${infile}", "-n${out_postfix}", "-p${seed}");
	if ($nrepbootstrap){ push @args, ("-fa", "-#${nrepbootstrap}", "-x${bseed}") }
	print STDERR join(' ', @args), "\n";
	system(@args);
	if ($? == -1){
		print "failed to execute: $!\n";
		}
	elsif ($? & 127){
		printf "child died with signal %d, %s coredump\n",
		($? & 127), ($? & 128) ? 'with' : 'without';
		}
	else {
		printf "child exited with value %d\n", $? >> 8;
		}

	my $treefile="RAxML_bestTree" . "\." . $out_postfix;
	my $bootstrapfile="RAxML_bootstrap" . "\." . $out_postfix;
	my $rootedtreefile="RAxML_ML_rooted.txt";
	my $rootedbsfile="RAxML_bootstrap_rooted.txt";
	my $treestring;
	if ($nrepbootstrap){
		system("nw_reroot -l $bootstrapfile @$ref_outgroups 1>$rootedbsfile 2>/dev/null");	# reroot bootstrap trees with outgroup taxa.
		system("nw_reroot -l $treefile @$ref_outgroups 1>$rootedtreefile 2>/dev/null");	# reroot ML tree with outgroup taxa.
		$treestring=`nw_support $rootedtreefile $rootedbsfile 2>/dev/null`;	# annotate ML tree.
		}
	else {
		$treestring=`nw_reroot -l $treefile @$ref_outgroups 2>/dev/null`;	# reroot tree with outgroup taxa.
		}
	my @outfiles=<*.${out_postfix}>;
	print STDERR "Removing ", join(', ', @outfiles), " ...\n";
	unlink @outfiles;	# remove all output files so that RAxML does not complain if file already exists.
	print STDERR "$treestring\n";

	return(\$treestring, $rootedtreefile, $rootedbsfile);
	}


sub callMPest{
	my ($self, $workdir, $treefile, $ctrlfile, $ref_taxalist, $ref_outgroups)=@_;

	my %taxanumbers=map { $ref_taxalist->[$_]=>($_+1) } (0..$#{ $ref_taxalist });
	my @outgroups=map { $taxanumbers{$_} } (@$ref_outgroups);
	die "Outgroups not present in taxalist!\n" unless (@outgroups==@$ref_outgroups);

	chdir($workdir);
	my @args=("mpest", $ctrlfile);
	print STDERR join(' ', @args), "\n";
	system(@args);
	if ($? == -1){
		print "failed to execute: $!\n";
		}
	elsif ($? & 127){
		printf "child died with signal %d, %s coredump\n",
		($? & 127), ($? & 128) ? 'with' : 'without';
		}
	else {
		printf "child exited with value %d\n", $? >> 8;
		}

	my $outfile=$treefile . ".tre";
	my $treestring=qx(grep 'tree mpest' $outfile | grep -o '\(.*\);' | nw_reroot -l - @outgroups 2>/dev/null);	# reroot tree with outgroup taxa.
	$treestring=~s/([\(\,])$_:/${1}$ref_taxalist->[$_-1]:/ foreach ( 1..scalar(@$ref_taxalist) );	# to avoid replacing 12 with 2.
	print STDERR "$treestring\n";

	return(\$treestring);
	}


sub processTrees{
	my ($self, $ref_treestring, $ref_results, $ref_comparisons)=@_;

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");
	my $tree;

	while($tree=$input->next_tree){
		my @testnodes;
		for my $comp (@$ref_comparisons){
			push @testnodes, [ grep { defined $_->id && ($_->id =~ /$comp->[0]_/ || $_->id =~ /$comp->[1]_/) } $tree->get_nodes ];
			}

		my @outnodes=grep { defined $_->id && ($_->id =~ /HYP_/ || $_->id =~ /parv_/) } $tree->get_nodes;
		my $outgroup=$tree->get_lca(-nodes => \@outnodes );

		for my $comp (0..$#testnodes){
			if ($tree->is_paraphyletic(-nodes => $testnodes[$comp], -outgroup => $outgroup) == 0 ){ ++$ref_results->[$comp+1] }
			++$ref_results->[0];
#			if ($tree->is_paraphyletic(-nodes => $ref_nodes, -outgroup => $outgroup) > 0 ){
#				print "$comp is paraphyletic for outgroup ", $outgroup->id || "unlabelled node", "\n";
#				}
#			elsif ($tree->is_paraphyletic(-nodes => $ref_nodes, -outgroup => $outgroup) < 0 ){ print "$comp is not monophyletic for outgroup ", $outgroup->id, "\n" }
#			else {
#				++$ref_results->{$comp};
#				print "$comp is not paraphyletic for outgroup ", $outgroup->id || "unlabelled node", "\n";
#				}
			}
		}
	return;
	}


sub processTrees4Species{
	my ($self, $ref_treestring, $ref_results, $ref_grouplabels)=@_;

	# $ref_results->[ 3 unrooted trees x 5 root positions (A,B,C,D,middle), non monophyletic ingroups, total ntrees ]

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");

	my $outgroup1=$ref_grouplabels->[-1];
	my $outgroup2=$ref_grouplabels->[-2];
	my @ingroups=@{ $ref_grouplabels }[0..3];

	my $tree;
	my $ntrees=1;

	TREE: while($tree=$input->next_tree){
		my @innodes=grep { defined $_->id && $_->id !~ /${outgroup1}_/ && $_->id !~ /${outgroup2}_/ } $tree->get_leaf_nodes;
		my $root=$tree->get_root_node;
		my @outnodes2=grep { defined $_->id && $_->id =~ /${outgroup2}_/ } $tree->get_leaf_nodes;
		my $outnode2;
		if (@outnodes2>1){ $outnode2=$tree->get_lca(-nodes => \@outnodes2) }
		else { $outnode2=$outnodes2[0] }
		print "Outnode2: ", $outnode2->id, "\n";

		unless ( $tree->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
			warn "Tree $ntrees is not monophyletic for ingroups in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my @speciesnodes;
		for my $species (0..3){
			push @speciesnodes, [ grep { defined $_->id && $_->id =~ /$ingroups[$species]_/ } $tree->get_leaf_nodes ];
			}
		my (@pair1, @pair2);
		my $pairnodes;
		FIRST: for my $group1 (0..2){
			SECOND: for my $group2 ($group1+1..3){
				$pairnodes=[ @{ $speciesnodes[$group1] }, @{ $speciesnodes[$group2] } ];
				print "Pairnodes: ", join(",", map { $_->id } @$pairnodes), "\n";
				if ($tree->is_paraphyletic(-nodes => $pairnodes, -outgroup => $outnode2) == 0 ){
					@pair1=($group1, $group2);
					@pair2=grep { $_ != $group1 && $_ != $group2 } (0..3);
					print "Identified pairings for tree $ntrees: ", join(',', @pair1), " and ", join(',', @pair2), "\n";
					last FIRST;
					}
				}
			}
		unless (@pair1==2 && @pair2==2){
			warn "No species pairing found for tree $ntrees in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my $testnodes=[ @{ $speciesnodes[ $pair2[0] ] }, @{ $speciesnodes[ $pair2[1] ] } ];
		print "Testnodes: ", join(",", map { $_->id } @$testnodes), "\n";
		if ($tree->is_paraphyletic(-nodes => $testnodes, -outgroup => $outnode2) == 0 ){
			warn "Tree $ntrees is balanced with species pairings ", join(',', @pair1), " and ", join(',', @pair2), " in:\n$$ref_treestring\n";
			my $index=$self->findTreeConf($pair1[0], $pair1[1], 4);
			++$ref_results->[$index];
			next TREE;
			}
		else{
			warn "Tree $ntrees is not balanced, with internal species pairing ", join(',', @pair1), ".\n";
			for my $species ($pair2[0], $pair2[1]){
				my $trippelnodes=[ @$pairnodes, @{ $speciesnodes[$species] } ];
				if ($tree->is_paraphyletic(-nodes => $trippelnodes, -outgroup => $outnode2) == 0 ){
					my $rootpos=( $species==$pair2[0] ) ? $pair2[1] : $pair2[0];
					my $index=$self->findTreeConf($pair1[0], $pair1[1], $rootpos);
					++$ref_results->[$index];
					next TREE;
					}
				}
			}
		warn "No tree topology found for tree $tree in:\n$$ref_treestring\n";
		++$ref_results->[15];
		} continue { ++$ntrees; ++$ref_results->[16] }
	return;
	}


sub processTrees6Ind{
	my ($self, $ref_treestring, $ref_results, $ref_grouplabels)=@_;

	# $ref_results->[ 3 unrooted trees x 5 root positions (A,B,C,D,middle), non monophyletic ingroups, total ntrees ]

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");

	my $outgroup1=$ref_grouplabels->[-1];
	my $outgroup2=$ref_grouplabels->[-2];
	my @ingroups=@{ $ref_grouplabels }[0..3];

	my $tree;
	my $ntrees=1;

	TREE: while ($tree=$input->next_tree){
		my @innodes=grep { defined $_->id && $_->id !~ /${outgroup1}_/ && $_->id !~ /${outgroup2}_/ } $tree->get_leaf_nodes;
		my $lcanode=$tree->get_lca(-nodes => \@innodes);
		my $lcaheight=$lcanode->height;
		my $root=$tree->get_root_node;
		my $outnode1=$tree->find_node(-id => "${outgroup1}_ind0");
		my $outnode2=$tree->find_node(-id => "${outgroup2}_ind0");
#		print "Outnode1: ", $outnode1->id, "\nOutnode2: ", $outnode2->id, "\n";

		unless ( $outnode1->ancestor->internal_id == $root->internal_id ){
			warn "Tree $ntrees is not correctly rooted:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		unless ( $tree->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
			warn "Tree $ntrees is not monophyletic for ingroups in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my @speciesnodes;
		for my $species (0..3){
			push @speciesnodes, $tree->find_node(-id => "${ingroups[$species]}_ind0");
			}

		my $lcadiv=0;	# get average distance from ingroups lca to ingroup tips.
		$lcadiv+=$tree->distance(-nodes => [ $lcanode, $_ ] ) foreach (@speciesnodes);
		$lcadiv/=scalar(@speciesnodes);
		my $outdiv=0;	# get average distance from outgroup to ingroup tips.
		$outdiv+=$tree->distance(-nodes => [ $outnode1, $_ ] ) foreach (@speciesnodes);
		$outdiv/=scalar(@speciesnodes);
		$ref_results->[17]+=$lcaheight;
		$ref_results->[18]+=$lcadiv;
		$ref_results->[19]+=$outdiv;

		my (@pair1, @pair2);
		my $pairancestor;
		FIRST: for my $group1 (0..2){
			SECOND: for my $group2 ($group1+1..3){
				if ( $speciesnodes[$group1]->ancestor->internal_id == $speciesnodes[$group2]->ancestor->internal_id ){
					$pairancestor=$speciesnodes[$group1]->ancestor;
					@pair1=($group1, $group2);
					@pair2=grep { $_ != $group1 && $_ != $group2 } (0..3);
					print STDERR "Identified pairings for tree $ntrees: ", join(',', @pair1), " and ", join(',', @pair2), "\n";
					last FIRST;
					}
				}
			}
		unless (@pair1==2 && @pair2==2){
			warn "No species pairing found for tree $ntrees in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		if ( $speciesnodes[ $pair2[0] ]->ancestor->internal_id == $speciesnodes[ $pair2[1] ]->ancestor->internal_id ){
			print STDERR "Tree $ntrees is balanced with species pairings ", join(',', @pair1), " and ", join(',', @pair2), "\n";
			my $index=$self->findTreeConf($pair1[0], $pair1[1], 4);
			++$ref_results->[$index];
			next TREE;
			}
		else{
			print STDERR "Tree $ntrees is not balanced, with internal species pairing ", join(',', @pair1), ".\n";
			for my $species ($pair2[0], $pair2[1]){
				if ( $pairancestor->ancestor->internal_id == $speciesnodes[$species]->ancestor->internal_id ){
					my $rootpos=( $species==$pair2[0] ) ? $pair2[1] : $pair2[0];
					my $index=$self->findTreeConf($pair1[0], $pair1[1], $rootpos);
					++$ref_results->[$index];
					next TREE;
					}
				}
			}
		warn "No tree topology found for tree $ntrees in:\n$$ref_treestring\n";
		++$ref_results->[15];
		} continue { ++$ntrees; ++$ref_results->[16] }
	return;
	}


sub processTrees5IndBootstrap{
	my ($self, $ref_treestring, $ref_results, $ref_grouplabels, $ref_bsthresholds)=@_;

	# $ref_results->[ 1 unrooted trees x 3 root positions (A,B,C), non monophyletic ingroups, total ntrees ]
	# BC,A; AC,B; AB,C

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick", -internal_node_id => 'bootstrap');

	my $outgroup1=$ref_grouplabels->[-1];
	my $outgroup2=$ref_grouplabels->[-2];
	my @ingroups=@{ $ref_grouplabels }[0..2];

	my $tree;
	my $ntrees=1;

	TREE: while ($tree=$input->next_tree){
		my @speciesnodes;
		for my $species (0..2){
			push @speciesnodes, $tree->find_node(-id => "$ingroups[$species]");
			}
		unless (@speciesnodes==3){ warn "Tree $ntrees has more than one leaf per species:\n$$ref_treestring\n"; next TREE }
		my $root=$tree->get_root_node;
		my $outnode1=$tree->find_node(-id => "$outgroup1");
		my $outnode2=$tree->find_node(-id => "$outgroup2");
		my $roottestnode;
		if ($outnode1->ancestor->internal_id == $outnode2->ancestor->internal_id){ $roottestnode=$outnode1->ancestor }
		else { $roottestnode=$outnode1 }

		unless ( $roottestnode->ancestor->internal_id == $root->internal_id ){
			warn "Tree $ntrees is not correctly rooted:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		unless ( $tree->is_monophyletic(-nodes => \@speciesnodes, -outgroup => $outnode2) ){
			warn "Tree $ntrees is not monophyletic for ingroups in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my (@pair1, @basal);
		my $bootstrap=0;
		my $pairancestor;
		FIRST: for my $group1 (0..1){
			SECOND: for my $group2 ($group1+1..2){
				if ( $speciesnodes[$group1]->ancestor->internal_id == $speciesnodes[$group2]->ancestor->internal_id ){
					$pairancestor=$speciesnodes[$group1]->ancestor;
					@pair1=($group1, $group2);
					@basal=grep { $_ != $group1 && $_ != $group2 } (0..2);
					$bootstrap=$pairancestor->bootstrap;
					last FIRST;
					}
				}
			}
		if (@pair1==2 && @basal==1 && $pairancestor->ancestor->internal_id == $speciesnodes[ $basal[0] ]->ancestor->internal_id){
			++$ref_results->[ $basal[0] ];
			$ref_results->[$basal[0]+17]+=$bootstrap;
			$ref_results->[32]+=$bootstrap;
			$ref_results->[33]+=$bootstrap;
			++$ref_results->[$basal[0]+34] if ($bootstrap>=$ref_bsthresholds->[1]);
			print STDERR "Tree $ntrees has an internal species pairing ", join(',', @pair1), " (bootstrap support: $bootstrap) and basal group $basal[0].\n";
			}
		else {
			warn "No species pairing found for tree $ntrees in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}
		} continue { ++$ntrees; ++$ref_results->[16] }
	return;
	}

sub processTrees6IndBootstrap{
	my ($self, $ref_treestring, $ref_results, $ref_grouplabels, $ref_bsthresholds)=@_;

	# $ref_results->[ 3 unrooted trees x 5 root positions (A,B,C,D,middle), non monophyletic ingroups, total ntrees ]

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick", -internal_node_id => 'bootstrap');

	my $outgroup1=$ref_grouplabels->[-1];
	my $outgroup2=$ref_grouplabels->[-2];
	my @ingroups=@{ $ref_grouplabels }[0..3];

	my $tree;
	my $ntrees=1;

	TREE: while ($tree=$input->next_tree){
		my @speciesnodes;
		for my $species (0..3){
			push @speciesnodes, $tree->find_node(-id => "$ingroups[$species]");
			}
		unless (@speciesnodes==4){ warn "Tree $ntrees has more than one leaf per species:\n$$ref_treestring\n"; next TREE }
		my $root=$tree->get_root_node;
		my $outnode1=$tree->find_node(-id => "$outgroup1");
		my $outnode2=$tree->find_node(-id => "$outgroup2");
		my $roottestnode;
		if ($outnode1->ancestor->internal_id == $outnode2->ancestor->internal_id){ $roottestnode=$outnode1->ancestor }
		else { $roottestnode=$outnode1 }

		unless ( $roottestnode->ancestor->internal_id == $root->internal_id ){
			warn "Tree $ntrees is not correctly rooted:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		unless ( $tree->is_monophyletic(-nodes => \@speciesnodes, -outgroup => $outnode2) ){
			warn "Tree $ntrees is not monophyletic for ingroups in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my (@pair1, @pair2);
		my ($bootstrap1, $bootstrap2)=(0, 0);
		my $pairancestor;
		FIRST: for my $group1 (0..2){
			SECOND: for my $group2 ($group1+1..3){
				if ( $speciesnodes[$group1]->ancestor->internal_id == $speciesnodes[$group2]->ancestor->internal_id ){
					$pairancestor=$speciesnodes[$group1]->ancestor;
					@pair1=($group1, $group2);
					@pair2=grep { $_ != $group1 && $_ != $group2 } (0..3);
					$bootstrap1=$pairancestor->bootstrap;
					print STDERR "Identified pairings for tree $ntrees: ", join(',', @pair1), " and ", join(',', @pair2), "\n";
					last FIRST;
					}
				}
			}
		unless (@pair1==2 && @pair2==2){
			warn "No species pairing found for tree $ntrees in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		if ( $speciesnodes[ $pair2[0] ]->ancestor->internal_id == $speciesnodes[ $pair2[1] ]->ancestor->internal_id ){
			$bootstrap2=$speciesnodes[ $pair2[0] ]->ancestor->bootstrap;
			print STDERR "Tree $ntrees is balanced with species pairings ", join(',', @pair1), " (bootstrap support: $bootstrap1) and ", join(',', @pair2), " (bootstrap support: $bootstrap2)\n";
			my $index=$self->findTreeConf($pair1[0], $pair1[1], 4);
			++$ref_results->[$index];
			$ref_results->[$index+17]+=($bootstrap1+$bootstrap2)/2;
			$ref_results->[32]+=($bootstrap1+$bootstrap2)/2;
			$ref_results->[33]+=($bootstrap1<$bootstrap2) ? $bootstrap1 : $bootstrap2;
			++$ref_results->[$index+34] if ($bootstrap1>=$ref_bsthresholds->[1] && $bootstrap2>=$ref_bsthresholds->[1]);
			for my $thres (0..$#{ $ref_bsthresholds }){ ++$ref_results->[$thres+49] if ($bootstrap1>=$ref_bsthresholds->[$thres] && $bootstrap2>=$ref_bsthresholds->[$thres]) }
			next TREE;
			}
		else{
			for my $species ($pair2[0], $pair2[1]){
				if ( $pairancestor->ancestor->internal_id == $speciesnodes[$species]->ancestor->internal_id ){
					$bootstrap2=$speciesnodes[$species]->ancestor->bootstrap;
					print STDERR "Tree $ntrees is not balanced, with internal species pairing ", join(',', @pair1), " (bootstrap support: $bootstrap1). Internal tripplet ", join(',', (@pair1, $species)), " (bootstrap support: $bootstrap2)\n";
					my $rootpos=( $species==$pair2[0] ) ? $pair2[1] : $pair2[0];
					my $index=$self->findTreeConf($pair1[0], $pair1[1], $rootpos);
					++$ref_results->[$index];
					$ref_results->[$index+17]+=($bootstrap1+$bootstrap2);
					$ref_results->[32]+=($bootstrap1+$bootstrap2);
					$ref_results->[33]+=($bootstrap1<$bootstrap2) ? $bootstrap1 : $bootstrap2;
					++$ref_results->[$index+34] if ($bootstrap1>=$ref_bsthresholds->[1] && $bootstrap2>=$ref_bsthresholds->[1]);
					for my $thres (0..$#{ $ref_bsthresholds }){ ++$ref_results->[$thres+49] if ($bootstrap1>=$ref_bsthresholds->[$thres] && $bootstrap2>=$ref_bsthresholds->[$thres]) }
					next TREE;
					}
				}
			}
		warn "No tree topology found for tree $ntrees in:\n$$ref_treestring\n";
		++$ref_results->[15];
		} continue { ++$ntrees; ++$ref_results->[16] }
	return;
	}

sub processTreesAll{	# under construction!!!
	my ($self, $outfolder, $ref_treestring, $ref_results, $ref_grouplabels, $ref_outgroups)=@_;

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");
	my $tree;

	my @outgroups=@{ $ref_grouplabels }[ @$ref_outgroups ];
	my @ingroups=@$ref_grouplabels;
	splice( @ingroups, -(scalar @$ref_outgroups) );	# remove outgroups from list of groups.

	my @comparisons;
	for my $group1 ( 0..$#ingroups-1 ){
		for my $group2 ( $group1+1..$#ingroups ){ push @comparisons, ( $ingroups[$group1], $ingroups[$group2] ) }
		}
	my $ncomps=@comparisons;

	while($tree=$input->next_tree){
		my @testnodes;
		for my $group1 ( 0..$#ingroups-1 ){
			for my $group2 ( $group1+1..$#ingroups ){
				push @testnodes, [ grep { defined $_->id && ($_->id =~ /$ingroups[$group1]_/ || $_->id =~ /$ingroups[$group2]->[1]_/) } $tree->get_nodes ];
				}
			}

		my @outnodes;
		push @outnodes, ( grep { defined $_->id && ($_->id =~ /($_)_/) } $tree->get_nodes ) foreach @outgroups;
		my $outgroup=$tree->get_lca(-nodes => \@outnodes );

		COMP: for my $comp (0..$#testnodes){
			if ($tree->is_paraphyletic(-nodes => $testnodes[$comp], -outgroup => $outgroup) == 0 ){
				if ($tree->is_paraphyletic(-nodes => $testnodes[-($comp+1)], -outgroup => $outgroup) == 0 ){
					++$ref_results->[ ($comp*3)+1 ];	# tree is balanced.
					}
				last COMP;
				}
			else { print "Species $comparisons[$comp][0] and species $comparisons[$comp][1] do not form a monophyletic group.\n" }
			} continue { ++$ref_results->[0] }
		}
	return;
	}


sub processTreesBootstrap{
	my ($self, $ref_treestring, $ref_results, $ref_taxa, $ref_poplabels, $ref_species, $ref_specieslabels)=@_;

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");

	my @ingroups=@$ref_specieslabels;
	my $outgroup1=pop @ingroups;
	my $outgroup2=pop @ingroups;

	my $tree;
	my $ntrees=1;

	TREE: while ($tree=$input->next_tree){
		my @innodes=grep { defined $_->id && $_->id !~ /$outgroup1/ && $_->id !~ /$outgroup2/ } $tree->get_leaf_nodes;
		my $root=$tree->get_root_node;
		my $outnode1=$tree->find_node(-id => $outgroup1);
		my $outnode2=$tree->find_node(-id => $outgroup2);
#		print "Outnode1: ", $outnode1->id, "\nOutnode2: ", $outnode2->id, "\n";

		unless ( $outnode1->ancestor->internal_id == $root->internal_id ){
			warn "Tree $ntrees is not correctly rooted:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		unless ( $tree->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
			warn "Tree $ntrees is not monophyletic for ingroups in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		my @speciesnodes;
		for my $species (0..$#ingroups){
			my @nodes=$tree->find_node(-id => $ingroups[$species]);
			unless (@nodes){
				for my $pop (@{ $ref_species->[$species] }){
					my $poplabel=$ref_poplabels->[$pop];
					push @nodes, $tree->find_node(-id => $poplabel);
					}
				}

			if (@nodes>1){
				if ($tree->is_paraphyletic(-nodes => \@nodes, -outgroup => $outnode2) == 0 ){
					push @speciesnodes, $tree->get_lca(-nodes => \@nodes);
					++$ref_results->[$species+17];
					}
				else { warn "Species $ref_specieslabels->[$species] is not monophyletic, skipping tree!\n"; ++$ref_results->[15]; next TREE }					
				}
			elsif (@nodes){ push @speciesnodes, @nodes }
			else { warn "No nodes found for $ref_specieslabels->[$species], skipping tree!\n"; ++$ref_results->[15]; next TREE } 
			}

		my (@pair1, @pair2);
		my $pairancestor;
		FIRST: for my $group1 (0..2){
			SECOND: for my $group2 ($group1+1..3){
				if ( $speciesnodes[$group1]->ancestor->internal_id == $speciesnodes[$group2]->ancestor->internal_id ){
					$pairancestor=$speciesnodes[$group1]->ancestor;
					@pair1=($group1, $group2);
					@pair2=grep { $_ != $group1 && $_ != $group2 } (0..3);
					print STDERR "Identified pairings for tree $ntrees: ", join(',', @pair1), " and ", join(',', @pair2), "\n";
					last FIRST;
					}
				}
			}
		unless (@pair1==2 && @pair2==2){
			warn "No species pairing found for tree $ntrees in:\n$$ref_treestring\n";
			++$ref_results->[15];
			next TREE;
			}

		if ( $speciesnodes[ $pair2[0] ]->ancestor->internal_id == $speciesnodes[ $pair2[1] ]->ancestor->internal_id ){
			warn "Tree $ntrees is balanced with species pairings ", join(',', @pair1), " and ", join(',', @pair2), "\n";
			my $index=$self->findTreeConf($pair1[0], $pair1[1], 4);
			++$ref_results->[$index];
			next TREE;
			}
		else{
			warn "Tree $ntrees is not balanced, with internal species pairing ", join(',', @pair1), ".\n";
			for my $species ($pair2[0], $pair2[1]){
				if ( $pairancestor->ancestor->internal_id == $speciesnodes[$species]->ancestor->internal_id ){
					my $rootpos=( $species==$pair2[0] ) ? $pair2[1] : $pair2[0];
					my $index=$self->findTreeConf($pair1[0], $pair1[1], $rootpos);
					++$ref_results->[$index];
					next TREE;
					}
				}
			}
		warn "No tree topology found for tree $ntrees in:\n$$ref_treestring\n";
		++$ref_results->[15];
		} continue { ++$ntrees; ++$ref_results->[16] }
	return;
	}


sub calculateTreeStats{
	my ($self, $ref_treestring, $ref_poplabels, $ref_groups, $ref_grouplabels, $ref_outgroups)=@_;

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");

	my @outgroups=@{ $ref_grouplabels }[ @$ref_outgroups ];
	my @ingroups=@$ref_groups;
	splice( @ingroups, -(scalar @$ref_outgroups) );	# remove outgroups from list of groups.

	my $tree=$input->next_tree;
	my @innodes=grep { defined $_->id && $_->id !~ /${outgroups[-1]}_/ && $_->id !~ /${outgroups[-2]}_/ } $tree->get_leaf_nodes;
	my $root=$tree->get_root_node;	# this is the unproper root.
	my $lcanode=$tree->get_lca(-nodes => \@innodes);
	my @outnodes1=grep { defined $_->id && $_->id =~ /${outgroups[-1]}_/ } $tree->get_leaf_nodes;
	my @outnodes2=grep { defined $_->id && $_->id =~ /${outgroups[-2]}_/ } $tree->get_leaf_nodes;
	my $outlca=$tree->get_lca(-nodes => [ @outnodes1, @outnodes2 ]);

	my ($outnode1, $outnode2, $realroot);
	if ($outlca eq $root){
		print STDERR "Tree has paraphyletic outgroups.\n";
		$outnode1=@outnodes1>1 ? $tree->get_lca(-nodes => \@outnodes1): $outnodes1[0];
		$outnode2=@outnodes2>1 ? $tree->get_lca(-nodes => \@outnodes2): $outnodes2[0];
		$realroot=$outnode2->ancestor;
		}
	elsif ($outlca->ancestor eq $root){
		print STDERR "Tree has monophyletic outgroups.\n";
		$outnode1=$outlca;
		$outnode2=$outlca;
		$realroot=$lcanode;
		}
	else {
		warn "Tree is not correctly rooted:\n$$ref_treestring\n";
		warn $outlca->internal_id, "/", $outlca->ancestor->internal_id, " vs. ", $root->internal_id, "\n";
		return;
		}

	unless ( $realroot->ancestor eq $root ){
		warn "Wrong definition of real root!\n";
		warn $realroot->ancestor->internal_id, " vs. ", $root->internal_id, "\n";
		return;
		}

	unless ( $tree->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
		warn "Tree is not monophyletic for ingroups in:\n$$ref_treestring\n";
		return;
		}

	my $lcadiv=0;	# get average distance from ingroups lca to ingroup tips.
	$lcadiv+=$tree->distance(-nodes => [ $lcanode, $_ ] ) foreach (@innodes);
	$lcadiv/=scalar(@innodes);
	my $rootdiv=0;	# get average distance from proper root to ingroup tips.
	$rootdiv+=$tree->distance(-nodes => [ $realroot, $_ ] ) foreach (@innodes);
	$rootdiv/=scalar(@innodes);
	my ($outdiv, $n)=(0, 0);	# get average distance from outgroup to ingroup tips.
	for my $outnode (@outnodes1){
		$outdiv+=$tree->distance(-nodes => [ $outnode, $_ ] ) foreach (@innodes);
		++$n;
		}
	$outdiv/=$n;
	my @dist=($lcanode->height, $realroot->height, $root->height, $lcadiv, $rootdiv, $outdiv);
#	print join(',', @dist), "\n";

	my @speciesnodes;
	push @speciesnodes, \@innodes;	# first element is a reference to list of all ingroup leaf nodes.
	for my $ref_species (@ingroups){
		my @nodes;
		for my $poplabel ( @{ $ref_poplabels }[ @$ref_species ] ){
#			print "Poplabel: $poplabel\n";
			push @nodes, ( grep { defined $_->id && $_->id =~ /${poplabel}_/ } $tree->get_leaf_nodes );
			}
		push @speciesnodes, \@nodes;
		}

	my @gsi;
	my $maxdenom;
	SET: for my $ref_nodes (@speciesnodes){
#		print join(',', map { $_->id } @$ref_nodes), "\n";
		my $lca=$tree->get_lca(-nodes => $ref_nodes);
		my %pathnodes;
		my $denom=0;
		for my $node (@$ref_nodes){
			my $anc=$node;
			while ($anc ne $lca){
				$anc=$anc->ancestor;
				$pathnodes{$anc}=$anc;
				}
			}
		for my $node (values %pathnodes){
			$denom+=scalar($node->each_Descendent)-1;
			}
		unless (defined $maxdenom){	# the first iteration of the loop is over all ingroup nodes and sets the maximum denominator.
			$maxdenom=$denom;
			print "Maximum denominator: $maxdenom for ", scalar(@$ref_nodes), " taxa.\n";
			next SET;
			}
		my $gs=(@$ref_nodes-1) / $denom;
		my $mings=(@$ref_nodes-1) / $maxdenom;
		my $gsi=($gs-$mings) / (1-$mings);
		print "gs: ", sprintf("%.4f", $gs), ", mings: ", sprintf("%.4f", $mings), ", gsi: ", sprintf("%.4f", $gsi), "\n";
		push @gsi, $gsi;
		}

	return(\@dist, \@gsi);
	}


sub calculateTreeDxy{
	my ($self, $ref_treestring, $ref_poplabels, $ref_groups, $ref_grouplabels, $ref_outgroups)=@_;

	my $treeio=IO::String->new($$ref_treestring);
	my $input=new Bio::TreeIO(-fh => $treeio, -format => "newick");

	my @outgroups=@{ $ref_grouplabels }[ @$ref_outgroups ];
	my @ingroups=@$ref_groups;
	splice( @ingroups, -(scalar @$ref_outgroups) );	# remove outgroups from list of groups.

	my $tree=$input->next_tree;
	my @innodes=grep { defined $_->id && $_->id !~ /${outgroups[-1]}_/ && $_->id !~ /${outgroups[-2]}_/ } $tree->get_leaf_nodes;
	my $root=$tree->get_root_node;	# this is the unproper root.
	my $lcanode=$tree->get_lca(-nodes => \@innodes);
	my @outnodes1=grep { defined $_->id && $_->id =~ /${outgroups[-1]}_/ } $tree->get_leaf_nodes;
	my @outnodes2=grep { defined $_->id && $_->id =~ /${outgroups[-2]}_/ } $tree->get_leaf_nodes;
	my $outnode1=@outnodes1>1 ? $tree->get_lca(-nodes => \@outnodes1): $outnodes1[0];
	my $outnode2=@outnodes2>1 ? $tree->get_lca(-nodes => \@outnodes2): $outnodes2[0];
	my $realroot=$outnode2->ancestor;

	unless ( $outnode1->ancestor eq $root ){
		warn "Tree is not correctly rooted:\n$$ref_treestring\n";
		warn $outnode1->ancestor->internal_id, " vs. ", $root->internal_id, "\n";
		return;
		}

	unless ( $realroot->ancestor eq $root ){
		warn "Wrong definition of real root!\n";
		warn $realroot->ancestor->internal_id, " vs. ", $root->internal_id, "\n";
		return;
		}

	unless ( $tree->is_monophyletic(-nodes => \@innodes, -outgroup => $outnode2) ){
		warn "Tree is not monophyletic for ingroups in:\n$$ref_treestring\n";
		return;
		}

	my $rootdiv=0;	# get average distance from proper root to ingroup tips.
	$rootdiv+=$tree->distance(-nodes => [ $realroot, $_ ] ) foreach (@innodes);
	$rootdiv/=scalar(@innodes);
	my $outdiv=0;	# get average distance from outgroup to ingroup tips.
	$outdiv+=$tree->distance(-nodes => [ $outnode1, $_ ] ) foreach (@innodes);
	$outdiv/=scalar(@innodes);

	my @speciesnodes;
	for my $ref_species (@ingroups){
		my @nodes;
		for my $poplabel ( @{ $ref_poplabels }[ @$ref_species ] ){
			push @nodes, ( grep { defined $_->id && $_->id =~ /${poplabel}_/ } $tree->get_leaf_nodes );
			}
		push @speciesnodes, \@nodes;
		}

	my (@meandxy, @mindxy, @maxdxy);
	my (@meanreldxy, @minreldxy, @maxreldxy);
	GROUP1: for my $group1 (0..$#speciesnodes-1){
		GROUP2: for my $group2 ($group1..$#speciesnodes){
			my ($sum, $max, $ncomp)=(0, 0, 0);
			my $min;
			for my $node1 (@{ $speciesnodes[$group1] }){
				for my $node2 (@{ $speciesnodes[$group2] }){
					my $dist=$tree->distance(-nodes => [ $node1, $node2 ] );
					$min=$dist if (defined $min && $dist<$min);
					$max=$dist if ($dist>$max);
					$sum+=$dist;
					++$ncomp;
					}
				}
			if ($ncomp){
				push @meandxy, $sum/$ncomp;
				push @mindxy, $min;
				push @maxdxy, $max;
				push @meanreldxy, $sum/($ncomp*$outdiv);
				push @minreldxy, $min/$outdiv;
				push @maxreldxy, $max/$outdiv;
				}
			}
		}

	return(\@meandxy, \@mindxy, \@maxdxy, \@meanreldxy,\@minreldxy, \@maxreldxy);
	}


sub parseTreestring{	# under construction!
	my ($self, $ref_treestring, $ref_grouplabels, $ref_outgrouplabels)=@_;

	my $locus=-1;
	my $treestring=$$ref_treestring;
	if ( substr($treestring, 0, 1, "") ne '(' ){ die "Wrongly formated newick tree string!\n" }
	my $level=1;
	my $sampleid=999;
	my $bl=0;
	CURSOR: while ( $level ){
		my $cursor=substr($treestring, 0, 1, "");
		if ( $cursor eq '(' ){ ++$level }
		elsif ( $cursor eq ')' ){ --$level }
		elsif ( $cursor eq ',' ){ next CURSOR }
		elsif ( $cursor eq ':' ){ next CURSOR }
		else { }
		}
	return;
	}


sub findTreeConf{
	my ($self, $value1, $value2, $rootpos)=@_;
	unless ($value2>$value1){
		if ($value1==$value2){ die "Invalid species indices specified!\n" }
		else { ($value1, $value2)=($value2, $value1) }
		}
	my $tree;
	if ($value1==0){
		if ($value2==1){ $tree=0 }	#	AB,CD
		elsif ($value2==2){ $tree=1 }	# AC,BD
		elsif ($value2==3){ $tree=2 }	# AD,BC
		}
	elsif ($value1==1){
		if ($value2==2){ $tree=2 }	# AD,BC
		elsif ($value2==3){ $tree=1 }	# AC,BD
		}
	elsif ($value1==2){
		if ($value2==3){ $tree=0 }	# AB,CD
		}
	return( ($tree * 5) + $rootpos );
	}


sub printFasta{
	my ($self, $outfile, $ref_seq, $ref_names, $rowlength)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file $outfile!\n" }

	for my $pop ( 0..$#{$ref_seq} ){
		for my $ind (0..$#{$ref_seq->[$pop]} ){
			print $output ">", $ref_names->[$pop], "_ind", $ind, "\n";
			my $seqstring=$ref_seq->[$pop][$ind];
			my $seqlength=length($seqstring);
			while ($seqlength>0){
				print $output substr($seqstring, 0, $rowlength, ""), "\n";	# chew up sequence string.
				$seqlength-=$rowlength;
				}
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }

	return;
	}



1;

