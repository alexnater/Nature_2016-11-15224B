package SFS;

use strict;
use warnings;
use Classes::Misc;


sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	return $self;
	}


sub generateSFSperWindow{
	my ($self, $regions, $refseq, $genotypes, $ref_missing, $groups, $mincov, $minind, $folded, $haplotize, $filtercpg, $print_persite_report, $notbycat, $ref_regionstomask)=@_;

	print "Generating SFS ...\n";
	my $npops=$ref_missing->getNPop();
	my $popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		$rt+=$temp;
		}
	my @groupsizes;
	for my $ref_group (@$groups){
		push @groupsizes, Misc::sum(@{$minind}[ @$ref_group ]);
		}

	my @arraysizes;
	for my $group ( 0..$#{$groups} ){
		if ($haplotize && $folded){ $arraysizes[$group]=int($groupsizes[$group]/2) }
		elsif ($haplotize || $folded){ $arraysizes[$group]=$groupsizes[$group] }
		else { $arraysizes[$group]=2*$groupsizes[$group] }
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		my ($prevend, $firststart);
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my (@sfs, @validpos);
			for my $group ( 0..$#{$groups} ){
				$sfs[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] for 0..4;
				$validpos[$group]=0;
				}
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($chrom) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($chrom) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $undefanc, $notbiallelic, $undefcat, $unknowncpg, $insufcov)=(0, 0, 0, 0, 0, 0);
			print "$chrom, $locus, $firststart, $ref_range->[0]-$ref_range->[-1]\n";
			my $prevpos=$start; my ($refbase, $prevbase);
			POSITION: for my $pos ( @$ref_range ){
				my $offset=$pos-$firststart;
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				$refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
				my ($ancestral, $derived, $cat);
				if ($snp){
					my $nalleles=$genotypes->NAlleles($chrom, $pos);
					if ($nalleles>2){ ++$notbiallelic; next POSITION }
					if (!$folded){
						$ancestral=$genotypes->getAncestral($chrom, $pos);
						$derived=$ancestral==0 ? 1 : 0;
						unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POSITION }
						if (!$notbycat){
							my $ancbase=$genotypes->getBase($chrom, $pos, $ancestral);
							my $derbase=$genotypes->getBase($chrom, $pos, $derived);
							$cat=Misc::getWeakStrong($ancbase, $derbase);
							if ($ancbase eq 'C' && $derbase eq 'T'){
								unless ($end-$pos){ ++$unknowncpg; next POSITION }	# last base of locus, next base not known.
								if ( defined $genotypes->getValue($chrom, $pos+1) ){ ++$unknowncpg; next POSITION }	# next base is polymorphic.
								elsif ( uc( substr($$ref_seqstring, 1, 1) ) eq 'G' ){ $cat=4 }
								}
							elsif ($ancbase eq 'G' && $derbase eq 'A'){
								unless ($pos-$start){ ++$unknowncpg; next POSITION }	# first base of locus, previous base not known.
								if ( defined $genotypes->getValue($chrom, $pos-1) ){ ++$unknowncpg; next POSITION }	# previous base was polymorphic.
								elsif ($prevbase eq 'C' ){ $cat=4 }
								}
							unless (defined $cat){ ++$undefcat; next POSITION }
							}
						else { $cat=0 }
						}
					else { $ancestral=0; $cat=0 }
					}
				GROUP: for my $metapop ( 0..$#{$groups} ){
					my $nder=0;
					for my $subpop ( @{$groups->[$metapop]} ){
						my @validalleles;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @validalleles, $allele;
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @validalleles, ($allele1, $allele2);
									}
								}
							elsif ($haplotize){ push @validalleles, 0 }
							else { push @validalleles, (0, 0) }
							}
						my $minnumber=$haplotize ? $minind->[$subpop] : 2*$minind->[$subpop];
						unless (@validalleles>=$minnumber){ ++$insufcov; next POSITION }
						if ($snp){
							my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $minnumber);
							for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
							}
						}
					if ($folded && $nder>$arraysizes[$metapop] ){ $nder=2*$arraysizes[$metapop]-$nder }
					if ($nder){ ++$sfs[$metapop][$cat][$nder] }
					++$validpos[$metapop];
					}
				} continue { $prevpos=$pos; $prevbase=$refbase }
			$regions->{_Loci}{$chrom}[$locus]{'SFS'}=\@sfs;
			$regions->{_Loci}{$chrom}[$locus]{'covered'}=\@validpos;
			$prevend=$end;
			print "$chrom:$start-$end; hardmasked: $hardmasked, undefined ancestral state: $undefanc, not biallelic: $notbiallelic, undefined category: $undefcat, unkown cpg: $unknowncpg, insufficent coverage: $insufcov\n";
			}
		}
	return;
	}


sub printCatReport{
	my ($self, $outfile, $windowfile, $regions, $refseq, $ancref, $chrommap, $ref_genotypes, $popsizes, $groups, $mincov, $maxcov, $ref_minind, $haplotize, $subsample, $maskcpg, $countcpg, $ref_regionstomask, $ref_annotation)=@_;

	my @type_of_change=('W2W', 'W2S', 'S2W', 'S2S', 'CpG2TpG');
	if (defined $maxcov && $maxcov>255){
		warn "Maximum coverage is set to a value higher than 255. Current framework cannot store coverage values higher than 255. Maximum coverage is set to 255!\n";
		$maxcov=255;
		}
	else { $maxcov=255 } 

	open my $output, ">>", $outfile or die "Could not open output file $outfile!$!\n";
	open my $windowstats, ">>", $windowfile or die "Could not open output file $windowfile!$!\n";

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my ($startchrom, $endchrom, $trstart, $trend);
			if (defined $chrommap){ 
				($startchrom, $trstart)=$chrommap->findPosonChrom($id, $start);
				($endchrom, $trend)=$chrommap->findPosonChrom($id, $end);
				if (defined $startchrom && defined $endchrom && $startchrom ne $endchrom){ die "Start and end position of window $id:$start-$end are not on the same chromosome!\n" }
				}
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $undefanc, $undefcat, $unknowncpg, $filtered)=(0, 0, 0, 0, 0, 0);
			my $nsnps_tot=0;
			my $prevpos=$start;
			my ($refbase, $ancbase, $prev_refbase, $prev_ancbase);
			my %refbasecounts=( 'G' => 0, 'C' => 0, 'A' => 0, 'T' => 0, 'N' => 0 );
			my %ancbasecounts=( 'G' => 0, 'C' => 0, 'A' => 0, 'T' => 0, 'N' => 0 );
			my ($refcpgcount, $anccpgcount, $totrefcpg, $totanccpg)=(0, 0, 0, 0);
			my ($wasrefcpg, $wasanccpg)=(0, 0);

			POSITION: for my $pos ( @$ref_range ){
				my $gap=$pos-$prevpos;
				substr($$ref_seqstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				substr($$ref_ancstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				$refbase=uc( substr($$ref_seqstring, 0, 1) );
				$ancbase=uc( substr($$ref_ancstring, 0, 1) );

				if ($countcpg && $pos){
					if ($pos==$start){
						my $ref_prev_refbase=$refseq->randomAccess($id, $pos-1, $pos-1, 1);
						$prev_refbase=$$ref_prev_refbase;
						my $ref_prev_ancbase=$ancref->randomAccess($id, $pos-1, $pos-1, 1);
						$prev_ancbase=$$ref_prev_ancbase;
						}
					++$totrefcpg unless ($prev_refbase eq 'N' || $refbase eq 'N');
					++$totanccpg unless ($prev_ancbase eq 'N' || $ancbase eq 'N');
					if ($refbase eq 'G' && $prev_refbase eq 'C'){ ++$refcpgcount; $wasrefcpg=1 }
					else { $wasrefcpg=0 }
					if ($ancbase eq 'G' && $prev_ancbase eq 'C'){ ++$anccpgcount; $wasanccpg=1 }
					else { $wasanccpg=0 }
					}

				if ($gap>1){ $prev_refbase="N"; $prev_ancbase="N"; $hardmasked+=$gap-1; $refbasecounts{'N'}+=$gap-1; $ancbasecounts{'N'}+=$gap-1 }
				if ($refbase eq 'N' ){ ++$hardmasked; $ancbase='N'; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.
				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				++$nsnps_tot;

				my $ancallele=$ref_genotypes->getAlleleforBase($id, $pos, $ancbase);
				my ($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($id, $pos, $ancallele);
				unless (defined $ancestral){ ++$undefanc; $ancbase='N'; next POSITION }
				unless (@derived){ $derived[0]=!$ancestral }
				if (@derived>1){ ++$notbiallelic; next POSITION }	# maximum one derived allele permitted.

#				$ancestral=$ref_genotypes->getAncestral($id, $pos);
#				if (!defined $ancestral && $ancbase ne 'N'){	# if ancestral state cannot be determined based on outgroup vcf genotypes.
#					my @bases=$ref_genotypes->getBase($id, $pos);
#					foreach (0..$#bases){ if ($bases[$_] eq $ancbase){ $ancestral=$_; last } }
#					}
#				my $ancgenobase=$ref_genotypes->getBase($id, $pos, $ancestral);
#				if ($ancgenobase ne $ancbase){
#					if ($ancgenobase ne 'N' && $ancbase ne 'N'){
#						warn "$id:$pos - Ancestral bases don't match: $ancgenobase vs. $ancbase! Skipping position.\n";
#						++$undefanc; $ancbase='N';
#						next POSITION;
#						}
#					else { $ancbase=$ref_genotypes->getBase($id, $pos, $ancestral) }
#					}

				my $derbase=$ref_genotypes->getBase($id, $pos, $derived[0]);
				my $cat=Misc::getWeakStrong($ancbase, $derbase);
				if ($maskcpg){
					if ($ancbase eq 'C' && $derbase eq 'T'){
						unless ($end-$pos){ ++$unknowncpg; next POSITION }	# last base of locus, next base not known.
						if ( defined $ref_genotypes->getValue($id, $pos+1) ){ ++$unknowncpg; next POSITION }	# next base is polymorphic.
						elsif ( uc( substr($$ref_seqstring, 1, 1) ) eq 'G' ){ $cat=4 }
						}
					elsif ($ancbase eq 'G' && $derbase eq 'A'){
						unless ($pos-$start && defined $prev_ancbase){ ++$unknowncpg; next POSITION }	# first base of locus, previous base not known.
						if ( defined $ref_genotypes->getValue($id, $pos-1) ){ ++$unknowncpg; next POSITION }	# previous base was polymorphic.
						elsif ($prev_ancbase eq 'C' ){ $cat=4 }
						}
					}
				unless (defined $cat){ ++$undefcat; $refbase='N'; $ancbase='N'; next POSITION }

				my @nchrom=(0) x @$groups;
				my @nder=(0) x @$groups;
				my $ndertot=0;
				my $onevalidgroup=0;

				GROUP: for my $metapop ( 0..$#{$groups} ){
					for my $subpop ( @{$groups->[$metapop]} ){
						my @validalleles;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov && $coverage<=$maxcov){ next IND }
							if ($haplotize){
								my $index=2*$ind+int(rand(2));
								my $allele=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele && $allele ne '.'){ next IND }
								push @validalleles, $allele;
								}
							else {
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								++$index;
								my $allele2=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								push @validalleles, ($allele1, $allele2);
								}
							}

						my $minnumber=$haplotize ? $ref_minind->[$subpop] : 2*$ref_minind->[$subpop];
						unless (@validalleles>=$minnumber){ $nchrom[$metapop]=0; $nder[$metapop]=0; next GROUP }
						my $ref_suballeles;
						if ($subsample && @validalleles!=$minnumber){ $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $minnumber) }
						else { $ref_suballeles=\@validalleles }
						for my $allele (@$ref_suballeles){ ++$nder[$metapop] if ($allele ne $ancestral) }
						$nchrom[$metapop]+=@$ref_suballeles;
						}
					$onevalidgroup=1;
					}
				unless ($onevalidgroup){ ++$filtered; next POSITION }	# at least one group needs to be covered.
				$ndertot=Misc::sum(@nder);
				if ($ndertot){	# only print if derived alleles present in any population.
					my ($trchrom, $trpos);
					if (defined $chrommap){
						($trchrom, $trpos)=$chrommap->findPosonChrom($id, $pos);
						}
					++$trpos if (defined $trpos);	# change to 1-based indexing. 
					my $nchromtot=Misc::sum(@nchrom);
					print $output "$id\t", $pos+1, "\t", $trchrom // 'NA', "\t", $trpos // 'NA', "\t$type_of_change[$cat]";
					print $output "\t", $ref_annotation->{_Annotations}{$id}{$pos} // 'NA' if (defined $ref_annotation);
					print $output "\t", sprintf("%.4f", $ndertot/$nchromtot), "\t$ndertot\t", $nchromtot/2;
					if (@$groups>1){ print $output "\t", $nchrom[$_] ? sprintf("%.4f", $nder[$_]/$nchrom[$_]) . "\t" . $nder[$_] : "NA\tNA", "\t", $nchrom[$_]/2 foreach ( 0..$#{$groups} ) }
					print $output "\n";	# scaffold    position    chromosome    position    type_of_change    derived_allele_frequency    derived_allele_count    number_of_individuals_covered
					}
				} continue {
					++$refbasecounts{$refbase}; ++$ancbasecounts{$ancbase}; $prevpos=$pos; $prev_refbase=$refbase; $prev_ancbase=$ancbase;
					}
			my ($chromstart, $chromend)=($trend>=$trstart) ? ($trstart, $trend+1) : ($trstart+1, $trend);
			print $windowstats "$id\t$start\t", $end+1, "\t", $startchrom // 'NA', "\t", $chromstart // 'NA', "\t", $chromend // 'NA';
			print $windowstats "\t$ancbasecounts{'G'}\t$ancbasecounts{'C'}\t$ancbasecounts{'A'}\t$ancbasecounts{'T'}\t$ancbasecounts{'N'}";
			print $windowstats "\t$refbasecounts{'G'}\t$refbasecounts{'C'}\t$refbasecounts{'A'}\t$refbasecounts{'T'}\t$refbasecounts{'N'}";
			my ($totGC, $totAT)=($ancbasecounts{'G'}+$ancbasecounts{'C'}, $ancbasecounts{'A'}+$ancbasecounts{'T'});
			print $windowstats "\t", ($totGC || $totAT) ? sprintf("%.4f", $totGC/($totGC+$totAT)*100) : "NA";
			($totGC, $totAT)=($refbasecounts{'G'}+$refbasecounts{'C'}, $refbasecounts{'A'}+$refbasecounts{'T'});
			print $windowstats "\t", ($totGC || $totAT) ? sprintf("%.4f", $totGC/($totGC+$totAT)*100) : "NA";
			print $windowstats "\t$anccpgcount\t$totanccpg\t", ($totanccpg) ? sprintf("%.4f", $anccpgcount/$totanccpg*100) : "NA" if ($countcpg);
			print $windowstats "\t$refcpgcount\t$totrefcpg\t", ($totrefcpg) ? sprintf("%.4f", $refcpgcount/$totrefcpg*100) : "NA" if ($countcpg);
			print $windowstats "\n";
			print "$id:$start-$end; total number of SNPs: $nsnps_tot, hardmasked: $hardmasked, not biallelic: $notbiallelic, undefined ancestral state: $undefanc, undefined category: $undefcat, unkown cpg: $unknowncpg, filtered: $filtered\n";
			}
		}
	close $output;
	close $windowstats;
	return;
	}


sub printAnnotations{
	my ($self, $outfile, $regions, $ref_annotation, $refseq, $ancref, $chrommap, $ref_genotypes, $popsizes, $groups, $mincov, $ref_haplotize, $ref_subsamples, $ref_regionstomask, $onlyannotated)=@_;

	open my $output, ">>", $outfile or die "Could not open output file $outfile!$!\n";
	my $nsets=(defined $ref_subsamples) ? @$ref_subsamples : 1;
	my @type_of_change=('W2W', 'W2S', 'S2W', 'S2S');

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my ($startchrom, $endchrom, $trstart, $trend);
			if (defined $chrommap){ 
				($startchrom, $trstart)=$chrommap->findPosonChrom($id, $start);
				($endchrom, $trend)=$chrommap->findPosonChrom($id, $end);
				if (defined $startchrom && defined $endchrom && $startchrom ne $endchrom){ die "Start and end position of window $id:$start-$end are not on the same chromosome!\n" }
				}
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $undefanc, $undefcat, $filtered)=(0, 0, 0, 0, 0);
			my $nsnps_tot=0;
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				my $gap=$pos-$prevpos;
				substr($$ref_seqstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				substr($$ref_ancstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1) );
				if ($gap>1){ $hardmasked+=$gap-1 }
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.

				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				next POSITION if ($onlyannotated && !exists $ref_annotation->{_Annotations}{$id}{$pos});	# only consider exonic sites.
				++$nsnps_tot;

				my $ancallele=$ref_genotypes->getAlleleforBase($id, $pos, $ancbase);
				my ($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($id, $pos, $ancallele);
				unless (defined $ancestral){ warn "$id:$pos - no ancestral allele for ancestral base $ancbase!\n"; ++$undefanc; next POSITION }
				$derived[0]=!$ancestral unless (@derived);
				if (@derived>1){ ++$notbiallelic; next POSITION }	# maximum one derived allele permitted.
				my $derbase=$ref_genotypes->getBase($id, $pos, $derived[0]);
				my $cat=Misc::getWeakStrong($ancbase, $derbase);
				unless (defined $cat){ ++$undefcat; next POSITION }

				my (@nchrom, @nder);
				my (@ancestral, @fixedder, @polym);
				for my $set (0..$nsets-1){
					@{$nchrom[$set]}=(0) x @$groups;
					@{$nder[$set]}=(0) x @$groups;
					@{$ancestral[$set]}=(0) x @$groups;
					@{$fixedder[$set]}=(0) x @$groups;
					@{$polym[$set]}=(0) x @$groups;
					}
				my @validgroups=(1) x @$groups;
				my $onevalidgroup=0;

				GROUP: for my $metapop ( 0..$#{$groups} ){
					for my $subpop ( @{$groups->[$metapop]} ){
						my $genotypestring=$ref_genotypes->getValue($id, $pos, $subpop);
						if (length($genotypestring)!=2*$popsizes->[$subpop] || $genotypestring=~/\./){ next GROUP }
						for my $ind (0..$popsizes->[$subpop]-1){
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next GROUP }
							}
						}

					SET: for my $set (0..$nsets-1){
						for my $subpop ( @{$groups->[$metapop]} ){
							my @alleles;
							my @subsample=(defined $ref_subsamples) ? @{$ref_subsamples->[$set][$subpop]} : (1..$popsizes->[$subpop]);
							IND: for my $ind (@subsample){
								if ($ref_haplotize->[$set]){
									push @alleles, $ref_genotypes->getValue($id, $pos, $subpop, 2*($ind-1)+int(rand(2)) );
									}
								else {
									push @alleles, ($ref_genotypes->getValue($id, $pos, $subpop, 2*$ind-2), $ref_genotypes->getValue($id, $pos, $subpop, 2*$ind-1) );
									}
								}

							my $minnumber=$ref_haplotize->[$set] ? scalar(@subsample) : 2*scalar(@subsample);
							unless (@alleles==$minnumber){ warn "$id:$pos - pop $subpop: Invalid size of allele array!\n"; next GROUP }
							for my $allele (@alleles){ ++$nder[$set][$metapop] if ($allele ne $ancestral) }
							$nchrom[$set][$metapop]+=@alleles;
							if (++$nder[$set][$metapop]==$minnumber){ $fixedder[$set][$metapop]=1 }
							elsif (++$nder[$set][$metapop]){ $polym[$set][$metapop]=1 }
							else { $ancestral[$set][$metapop]=1 }
							}
						}
					$onevalidgroup=1;
					}
				unless ($onevalidgroup){ ++$filtered; next POSITION }	# at least one set needs to be covered.
				my ($trchrom, $trpos);
				if (defined $chrommap){
					($trchrom, $trpos)=$chrommap->findPosonChrom($id, $pos);
					}
				++$trpos if (defined $trpos);	# change to 1-based indexing. 
				print $output "$id\t", $pos+1, "\t", $trchrom // 'NA', "\t", $trpos // 'NA',
					"\t$type_of_change[$cat]\t", $ref_annotation->{_Annotations}{$id}{$pos} // 'NA';
				for my $group ( 0..$#{$groups} ){
					for my $set (0..$nsets-1){
						if ($validgroups[$group]){
							print $output "\t", sprintf("%.4f", $nder[$set][$group]/$nchrom[$set][$group]), "\t$nder[$set][$group]\t", $nchrom[$set][$group]/2;
							}
						else { print $output "\tNA\tNA\tNA" }
						}
					}
				print $output "\n";	# scaffold    position    chromosome    position    type_of_change    derived_allele_frequency    derived_allele_count    number_of_individuals_covered
				} continue { $prevpos=$pos }
			print "$id:$start-$end; total number of SNPs: $nsnps_tot, hardmasked: $hardmasked, not biallelic: $notbiallelic, undefined ancestral state: $undefanc, undefined category: $undefcat, filtered: $filtered\n";
			}
		}
	close $output;
	return;
	}


sub printSNPDifferentials{
	my ($self, $outfile, $regions, $ref_annotation, $refseq, $ancref, $chrommap, $ref_genotypes, $ref_popsizes, $ref_groups, $ref_mingroupind, $mincov, $haplotize, $ref_regionstomask, $onlyannotated)=@_;

	die "Script only supports a single pair of groups!\n" unless (@$ref_groups==2);
	open my $output, ">>", $outfile or die "Could not open output file $outfile!$!\n";
	my @type_of_change=('W2W', 'W2S', 'S2W', 'S2S');

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my ($startchrom, $endchrom, $trstart, $trend);
			if (defined $chrommap){ 
				($startchrom, $trstart)=$chrommap->findPosonChrom($id, $start);
				($endchrom, $trend)=$chrommap->findPosonChrom($id, $end);
				if (defined $startchrom && defined $endchrom && $startchrom ne $endchrom){ die "Start and end position of window $id:$start-$end are not on the same chromosome!\n" }
				}
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $undefanc, $undefcat, $filtered)=(0, 0, 0, 0, 0);
			my $nsnps_tot=0;
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				my $gap=$pos-$prevpos;
				substr($$ref_seqstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				substr($$ref_ancstring, 0, $gap, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1) );
				if ($gap>1){ $hardmasked+=$gap-1 }
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral reference sequence is undefined.

				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				next POSITION if ($onlyannotated && !exists $ref_annotation->{_Annotations}{$id}{$pos});	# only consider exonic sites.
				++$nsnps_tot;

				my $ancallele=$ref_genotypes->getAlleleforBase($id, $pos, $ancbase);
				my ($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($id, $pos, $ancallele);
				unless (defined $ancestral){ warn "$id:$pos - no ancestral allele for ancestral base $ancbase!\n"; ++$undefanc; next POSITION }
				$derived[0]=!$ancestral unless (@derived);
				if (@derived>1){ ++$notbiallelic; next POSITION }	# maximum one derived allele permitted.
				my $derbase=$ref_genotypes->getBase($id, $pos, $derived[0]);
				my $cat=Misc::getWeakStrong($ancbase, $derbase);
				unless (defined $cat){ ++$undefcat; next POSITION }

				my @nchrom=(0) x @$ref_groups;
				my @nder=(0) x @$ref_groups;

				GROUP: for my $metapop ( 0..$#{$ref_groups} ){
					my @validalleles;
					for my $subpop ( @{$ref_groups->[$metapop]} ){
						IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($haplotize){
								my $index=2*$ind+int(rand(2));
								my $allele=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele && $allele ne '.'){ next IND }
								push @validalleles, $allele;
								}
							else {
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								++$index;
								my $allele2=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								push @validalleles, ($allele1, $allele2);
								}
							}
						}

					my $minnumber=$haplotize ? $ref_mingroupind->[$metapop] : 2*$ref_mingroupind->[$metapop];
					unless (@validalleles>=$minnumber){ ++$filtered; next POSITION }
					for my $allele (@validalleles){ ++$nder[$metapop] if ($allele ne $ancestral) }
					$nchrom[$metapop]=@validalleles;
					}

				my ($trchrom, $trpos);
				if (defined $chrommap){
					($trchrom, $trpos)=$chrommap->findPosonChrom($id, $pos);
					}
				++$trpos if (defined $trpos);	# change to 1-based indexing.
				my $higher_group=($nder[0]/$nchrom[0]>$nder[1]/$nchrom[1]) ? 0 : 1;
				$higher_group="-" if ($nder[0]/$nchrom[0]==$nder[1]/$nchrom[1]);
				print $output "$id\t", $pos+1, "\t", $trchrom // 'NA', "\t", $trpos // 'NA',
					"\t$ancbase\t$derbase\t$type_of_change[$cat]\t", $ref_annotation->{_Annotations}{$id}{$pos} // 'noncoding',
					"\t", sprintf("%.4f", abs($nder[0]/$nchrom[0]-$nder[1]/$nchrom[1]) ), "\t$higher_group";
				for my $group ( 0..$#{$ref_groups} ){
					print $output "\t", sprintf("%.4f", $nder[$group]/$nchrom[$group]), "\t$nder[$group]\t", $nchrom[$group]/2;
					}
				print $output "\n";	# scaffold    position    chromosome    position    type_of_change     allele_frequency_difference   group_index_with_higher_derived_frequency {derived_allele_frequency    derived_allele_count   number_of_individuals_covered}*2
				} continue { $prevpos=$pos }
			print "$id:$start-$end; total number of SNPs: $nsnps_tot, hardmasked: $hardmasked, not biallelic: $notbiallelic, undefined ancestral state: $undefanc, undefined category: $undefcat, filtered: $filtered\n";
			}
		}
	close $output;
	return;
	}



sub generateJointSFS{
	my ($self, $regions, $outfile, $genotypes, $ref_missing, $groups, $mincov, $minind, $haplotize)=@_;

	print "Generating SFS..\n";
	my $npops=$ref_missing->getNPop();
	my @sfs;	# $sfs[poppair][nder1][nder2];
	my (@popsizes, @popsum, @groupsizes);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}
	for my $metapop ( 0..$#{$groups} ){
		push @groupsizes, Misc::sum(@{$minind}[ @{$groups->[$metapop]} ]);
		}

	my $poppair=0; my @validpos; my @undefanc; my @notbiallelic;
	for my $pop1 ( 0..($#{$groups}-1) ){
		for my $pop2 ( ($pop1+1)..$#{$groups} ){
			print "Population pair $poppair; total size metapopulation 1: $groupsizes[$pop1], metapopulation 2: $groupsizes[$pop2]\n";
			my @arraysizes;
			if ($haplotize){ @arraysizes=($groupsizes[$pop1], $groupsizes[$pop2]) }
			else { @arraysizes=(2*$groupsizes[$pop1], 2*$groupsizes[$pop2]) }
			for my $i (0..$arraysizes[0]){ $sfs[$poppair][$i]=[ (0) x ($arraysizes[1]+1) ] }
			$validpos[$poppair]=0; $undefanc[$poppair]=0; $notbiallelic[$poppair]=0;

			CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
				LOCUS: for my $locus ( $regions->allKeys($chrom) ){
					my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
					my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
					my $locuslength=$end-$start+1;
#					my $iterator=$regions->getreftoIndividuals('locus', $chrom, $locus);
					print "$chrom:$start-$end\n";

					POSITION: for my $pos ($start..$end){
						my $offset=$pos-$start;
						my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
						my $ancestral;
						if ($snp){
							my $nalleles=$genotypes->NAlleles($chrom, $pos);
							if ($nalleles>2){++$notbiallelic[$poppair]; next POSITION}
							elsif ($nalleles<2){$snp=0}
							$ancestral=$genotypes->getAncestral($chrom, $pos);
							unless (defined $ancestral && $ancestral<=3){++$undefanc[$poppair]; next POSITION}
							}
						my @nder=( (0) x scalar(@$groups) );
						my $validindgroup=0;
						for my $metapop ($pop1, $pop2){
							for my $subpop (@{$groups->[$metapop]}){
								my @validalleles;
#								INDEX: for my $index (@{$iterator->[$subpop]}){
								IND: for my $ind (0..$popsizes[$subpop]-1){
									my $allind=$ind+$popsum[$subpop];
									my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
									unless (defined $coverage && $coverage>=$mincov){next IND}
									if ($snp){
										if ($haplotize){
											my $index=2*$ind+int(rand(2));
											my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele && $allele ne '.'){next IND}
											push @validalleles, $allele;
											}
										else {
											my $index=2*$ind;
											my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele1 && $allele1 ne '.'){next IND}
											++$index;
											my $allele2=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele2 && $allele2 ne '.'){next IND}
											push @validalleles, ($allele1, $allele2);
											}
										}
									elsif ($haplotize){push @validalleles, 0}
									else {push @validalleles, (0, 0)}
									}
								my $minnumber=$haplotize ? $minind->[$subpop] : 2*$minind->[$subpop];
								unless (@validalleles>=$minnumber){next POSITION}
								if ($snp){
									my @suballeles=( Misc::shuffle(@validalleles) )[0..$minnumber-1];
									for my $allele (@suballeles){ if ($allele ne $ancestral){++$nder[$metapop]} }
									}
								}
							}
#						print "$chrom:$pos - nder1: $nder[$pop1], nder2: $nder[$pop2].\n";						
						++$sfs[$poppair][ $nder[$pop1] ][ $nder[$pop2] ];
						++$validpos[$poppair];
						}
					}
				}
			$sfs[$poppair][0][0]+=$sfs[$poppair][ $arraysizes[0] ][ $arraysizes[1] ]; # the top right bin needs to be combined with the lower left bin.
			$sfs[$poppair][ $arraysizes[0] ][ $arraysizes[1] ]=0;

			my $output;
			if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", "${outfile}_${pop1}_${pop2}.txt" or die "Could not open output file!\n"} 
			print $output "Population pair: $poppair\tvalid positions: $validpos[$poppair]\tundefined ancestral state: $undefanc[$poppair]\tnot-biallelic: $notbiallelic[$poppair]\n";
			for my $nder1 (0..$arraysizes[0]){print $output "\td${pop1}_${nder1}"}
			print $output "\n";
			for my $nder2 (0..$arraysizes[1]){
				print $output "d${pop2}_${nder2}";
				for my $nder1 ( 0..$arraysizes[0]){print $output "\t$sfs[$poppair][$nder1][$nder2]"}
				print $output "\n";
				}
			if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
			++$poppair;
			}
		}
	return;
	}


sub storeSFS{
	my ($self, $regions, $refseq, $ancref, $ref_genotypes, $ref_missing, $npops, $ref_popsizes, $groups, $mincov, $minind, $haplotize, $snpsonly)=@_;

	print "Generating SFS ...\n";
	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_popsizes->[$pop];
		$rt+=$temp;
		}
	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;

	my @validpos;
	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		for my $pop (0..$#{$groups}){
			my $arraysize;
			if ($haplotize){ $arraysize=$groupsizes[$pop] }
			else { $arraysize=2*$groupsizes[$pop] }
			$self->{_SFS}[$pop]=[ (0) x ($arraysize+1) ];
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($chrom, $start, $end, 1);
			my ($undefanc, $notbiallelic, $hardmasked)=(0,0,0);
			POSITION: for my $pos ($start..$end){
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1, "") );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral sequence is not defined.
				my $offset=$pos-$start;
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				next POSITION if ($snpsonly && !$snp);
				my $ancestral;
				if ($snp){
					my $nalleles=$ref_genotypes->NAlleles($chrom, $pos);
					if ($nalleles>2){ ++$notbiallelic; next POSITION }
					$ancestral=$ref_genotypes->getAncestral($chrom, $pos);
					unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POSITION }
					}
				elsif ($refbase ne $ancbase){	# position fixed for reference base, but not the same as ancestral base, should not occur in our case.
					$ancestral=1;
					warn "$chrom:$pos - position not polymorphic, but reference base $refbase not the same as ancestral base $ancbase!\n";
					}
				else { $ancestral=0 }

				GROUP: for my $metapop ( 0..$#{$groups} ){
					my $nder=0;
					for my $subpop (@{$groups->[$metapop]}){
						my @validalleles;
						IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=defined $ref_missing ? $ref_missing->getValue($chrom, $start, $allind, $offset) : $ref_genotypes->getCoverage($chrom, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){next IND}
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){next IND}
									push @validalleles, $allele;
									}
								else {
									my $index=2*$ind;
									my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){next IND}
									my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.'){next IND}
									push @validalleles, ($allele1, $allele2);
									}
								}
							elsif ($haplotize){push @validalleles, 0}
							else {push @validalleles, (0, 0)}
							}
						my $minnumber=$haplotize ? $minind->[$subpop] : 2*$minind->[$subpop];
						if (@validalleles<$minnumber){ next GROUP }
						if ($snp){
							my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $minnumber);
							for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
							}
						elsif ($ancestral==1){ $nder=$minnumber }
						}
					++$self->{_SFS}[$metapop][$nder];
					++$validpos[$metapop];
					}
				}
			print STDERR "$chrom:$start-$end; hardmasked: $hardmasked, undefined ancestral state: $undefanc, not biallelic: $notbiallelic\n";
			}
		}
	return;
	}


sub storeJointSFS{
	my ($self, $posfile, $regions, $refseq, $ancref, $genotypes, $ref_missing, $groups, $mincov, $minind, $folded, $haplotize, $ref_regionstomask)=@_;

	print "Generating Joint SFS ...\n";
	$ancref=$refseq unless defined $ancref;
	my $npops=$ref_missing->getNPop();
	my $popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		$rt+=$temp;
		}
	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;
	my @arraysizes;
	for my $group ( 0..$#{$groups} ){
		if ($haplotize && $folded){ $arraysizes[$group]=int($groupsizes[$group]/2) }
		elsif ($haplotize || $folded){ $arraysizes[$group]=$groupsizes[$group] }
		else { $arraysizes[$group]=2*$groupsizes[$group] }
		}

	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		my $poppair=0;
		for my $pop1 ( 0..($#{$groups}-1) ){
			for my $pop2 ( ($pop1+1)..$#{$groups} ){
				for my $i (0..$arraysizes[$pop1]){ $self->{_SFS}[$poppair][$i]=[ (0) x ($arraysizes[$pop2]+1) ] }
				++$poppair;
				}
			}
		}

	my $output;
	open $output, ">>", $posfile or die "Could not open output file!\n" if (defined $posfile);

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		my ($end, $prevend, $firststart);
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			$end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_ancstring=$ancref->randomAccess($chrom, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($chrom) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($chrom) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $undefanc, $notbiallelic)=(0, 0, 0);
			print "$chrom, $locus, $firststart, $ref_range->[0]-$ref_range->[-1]\n";
			my $prevpos=$start;
			POSITION: for my $pos ( @$ref_range ){
				my $offset=$pos-$firststart;
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				substr($$ref_ancstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				my $ancbase=uc( substr($$ref_ancstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				if ($ancbase eq 'N' ){ ++$undefanc; next POSITION }	# discard if base in ancestral sequence is not defined.
				my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
				my ($ancestral, $derived);
				if ($snp){
					my $nalleles=$genotypes->NAlleles($chrom, $pos);
					if ($nalleles>2){ ++$notbiallelic; next POSITION }	# maybe only apply on a per poppair basis?
					if (!$folded){
						$ancestral=$genotypes->getAncestral($chrom, $pos);
						$derived=$ancestral==0 ? 1 : 0;
						unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POSITION }
						my $ancgenobase=$genotypes->getBase($chrom, $pos, $ancestral);
						warn "$chrom:$pos - ancestral reference base $ancbase not the same as ancestral genotype base $ancgenobase!\n" unless ($ancbase eq $ancgenobase);
						}
					else { $ancestral=0 }
					}
				elsif ($refbase ne $ancbase){	# position fixed for reference base, but not the same as ancestral base, should not occur in our case.
					$ancestral=1;
					warn "$chrom:$pos - position not polymorphic, but reference base $refbase not the same as ancestral base $ancbase!\n";
					}
				else { $ancestral=0 }
				my @validgroup=( (1) x scalar(@$groups) );
				my @nder=( (0) x scalar(@$groups) );
				GROUP: for my $metapop ( 0..$#{$groups} ){
					for my $subpop ( @{$groups->[$metapop]} ){
						my @validalleles;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @validalleles, $allele;
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									++$index;
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @validalleles, ($allele1, $allele2);
									}
								}
							elsif ($haplotize){ push @validalleles, $ancestral }
							else { push @validalleles, ($ancestral, $ancestral) }
							}
						my $minnumber=$haplotize ? $minind->[$subpop] : 2*$minind->[$subpop];
						if (@validalleles<$minnumber){ $validgroup[$metapop]=0; next GROUP }
						if ($snp){
							my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $minnumber);
							for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder[$metapop] } }
							}
						elsif ($ancestral==1){ $nder[$metapop]=$minnumber }
						}
					if ($folded && $nder[$metapop]>$arraysizes[$metapop] ){ $nder[$metapop]=2*$arraysizes[$metapop]-$nder[$metapop] }
					}
				my $poppair=0;
				POP1: for my $pop1 ( 0..($#{$groups}-1) ){
					next POP1 unless ($validgroup[$pop1]);
					POP2: for my $pop2 ( ($pop1+1)..$#{$groups} ){
						next POP2 unless ($validgroup[$pop2]);
						++$self->{_SFS}[$poppair][ $nder[$pop1] ][ $nder[$pop2] ];
#						print "Valid position $chrom:$pos for group $poppair.\n";
						} continue { ++$poppair }
					}
				if (defined $posfile){
					print $output "$chrom\t$pos";
					print $output "\t", ($validgroup[$_]) ? $nder[$_] : "N" foreach ( 0..($#{$groups}) );
					print $output "\n";
					}
				} continue { $prevpos=$pos }
			print "$chrom:$start-$end; hardmasked: $hardmasked, undefined ancestral state: $undefanc, not biallelic: $notbiallelic\n";
			} continue { $prevend=$end }
		}
	close $output if (defined $posfile);
	return;
	}


sub fixedsharedStats{	# Double-check before using. Bug corrected in jumping over invariant sites and refseq chewing.
	my ($self, $regions, $refseq, $ref_genotypes, $npops, $popsizes, $groups, $mincov, $ref_minind, $haplotize, $ref_regionstomask)=@_;

	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$popsizes->[$pop];
		$rt+=$temp;
		}
	my @popminsizes=map { $haplotize ? $_ : 2*$_ } @$ref_minind;
	my @groupminsizes=map { Misc::sum(@popminsizes[ @$_ ]) } @$groups;
	my $totsize=Misc::sum(@groupminsizes);

	unless ( exists $self->{_Fixed} ){	# to prevent overwriting of previously stored counts.
		my $nsets=2**scalar(@$groups);
		$self->{_Fixed}=[ (0) x $nsets ];
		$self->{_Shared}=[ (0) x $nsets ];
		}

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($id, $start, $end, 1);
			my $ref_range;
			if ( defined $ref_regionstomask && $ref_regionstomask->getValue($id) ){ $ref_range=Misc::maskInterval( $start, $end, $ref_regionstomask->getValue($id) ) }
			else { $ref_range=[ ($start..$end) ] }
			my ($hardmasked, $notbiallelic, $undefanc, $filtered)=(0, 0, 0, 0);
			my ($nsnps_tot, $nsnps_all, $nsnps_ancdef)=(0, 0, 0);
			my $prevpos=$start;

			POSITION: for my $pos ( @$ref_range ){
				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );	# only consider polymorphic sites.
				++$nsnps_tot;
#				substr($$ref_seqstring, 0, $pos-$prevpos, "");	# The chewing appeared twice. Must have messed up the refbase extraction!
				my $refbase=uc( substr($$ref_seqstring, 0, 1) );
				if ($refbase eq 'N' ){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my ($ancestral, $derived);
				my $nalleles=$ref_genotypes->NAlleles($id, $pos);
				if ($nalleles>2){ ++$notbiallelic; next POSITION };
				$ancestral=$ref_genotypes->getAncestral($id, $pos);
				my $ancdef=1;
				unless (defined $ancestral && $ancestral<=3){ ++$undefanc; $ancestral=0; $ancdef=0 }
				$derived=$ancestral==0 ? 1 : 0;
				my ($fixder, $polym)=(0, 0);
				my $ndertot=0;
				GROUP: for my $metapop ( 0..$#{$groups} ){
					my $nder=0;
					for my $subpop ( @{$groups->[$metapop]} ){
						my @validalleles;
						IND: for my $ind (0..$popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $subpop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($haplotize){
								my $index=2*$ind+int(rand(2));
								my $allele=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele && $allele ne '.'){ next IND }
								push @validalleles, $allele;
								}
							else {
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								++$index;
								my $allele2=$ref_genotypes->getValue($id, $pos, $subpop, $index);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								push @validalleles, ($allele1, $allele2);
								}
							}
						unless (@validalleles>=$popminsizes[$subpop]){ ++$filtered; next POSITION }
						my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, $popminsizes[$subpop]);
						for my $allele (@$ref_suballeles){ ++$nder if ($allele ne $ancestral) }
						}
					if ( $nder==$groupminsizes[$metapop] ){ $fixder |= (1 << $metapop) }
					elsif ($nder>0){ $polym |= (1 << $metapop) }
					$ndertot+=$nder;
					}

#				print "$id:$pos - $fixder, $polym\n";
				if ($ancdef){	# only score if ancestral state is defined.
					my $polder=$fixder | $polym;
					++$self->{_Fixed}[$fixder];
					++$self->{_Derived}[$polder];
					++$nsnps_ancdef if ($polder);	# only count snps that are variable within focal species.
					}
				++$self->{_Shared}[$polym];	# also score if ancestral state is undefined.
				++$nsnps_all if ($ndertot>0 && $ndertot<$totsize);	# only count snps that are variable within focal species.
				} continue { $prevpos=$pos }
			print "$id:$start-$end; number of valid SNPs: $nsnps_tot, hardmasked: $hardmasked, not biallelic: $notbiallelic, undefined ancestral state: $undefanc, filtered: $filtered\n";
			$self->{_nsnps_all}+=$nsnps_all;
			$self->{_nsnps_ancdef}+=$nsnps_ancdef;
			}
		}
	return;
	}


sub storeSFSIntron{
	my ($self, $regions, $refseq, $genotypes, $ref_missing, $groups, $mincov, $minind, $folded, $haplotize, $filtercpg)=@_;

	print "Generating SFS ...\n";
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum, @groupsizes);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	for my $ref_minset (@$minind){
		my @minsizes;
		for my $ref_group (@$groups){
			push @minsizes, Misc::sum(@{$ref_minset}[ @$ref_group ]);
			}
		push @groupsizes, \@minsizes;
		}	

	my (@validpos, $undefanc, $notbiallelic, $cpgfiltered, $totalsites);
	my (@arraysizes, @minnumbers);
	for my $minset ( 0..$#{$minind} ){
		for my $group ( 0..$#{$groups} ){
			if ($haplotize && $folded){ $arraysizes[$minset][$group]=int($groupsizes[$minset][$group]/2) }
			elsif ($haplotize || $folded){ $arraysizes[$minset][$group]=$groupsizes[$minset][$group] }
			else { $arraysizes[$minset][$group]=2*$groupsizes[$minset][$group] }
			}
		for my $pop ( 0..$#{$minind->[$minset]} ){
			if ($haplotize){ push @{$minnumbers[$pop]}, $minind->[$minset][$pop] }
			else { push @{$minnumbers[$pop]}, 2*$minind->[$minset][$pop] }
			}
		}

	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		for my $minset ( 0..$#{$minind} ){
			for my $group ( 0..$#{$groups} ){
				$self->{_SFS}[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
				}
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$refseq->randomAccess($chrom, $start, $end, 1);
			my $ref_cpgstring=Misc::getCpGString($ref_seqstring, $genotypes, $chrom, $start, $end);

			POSITION: for my $pos ($start..$end){
				++$totalsites;
				my $offset=$pos-$start;
				if ($filtercpg){ unless ( substr($$ref_cpgstring, 0, 1, "") ){ ++$cpgfiltered; next POSITION } }
				my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
				if ($snp){
					if ( $genotypes->NAlleles($chrom, $pos)>2 ){ ++$notbiallelic; next POSITION }
					}
				my $ancestral;
				if (!$folded && $snp){
					$ancestral=$genotypes->getAncestral($chrom, $pos);
					unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POSITION }
					}
				else { $ancestral=0 }
				GROUP: for my $metapop ( 0..$#{$groups} ){
					my $nder=0; my $sfsarray=0;
					my @validalleles;
					for my $subpop ( @{$groups->[$metapop]} ){
						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @{$validalleles[$subpop]}, $allele;
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									++$index;
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @{$validalleles[$subpop]}, ($allele1, $allele2);
									}
								}
							elsif ($haplotize){ push @{$validalleles[$subpop]}, 0 }
							else { push @{$validalleles[$subpop]}, (0, 0) }
							}
						my $valid=0;
						unless ( defined $validalleles[$subpop] ){ next GROUP }
						SET: for my $minset ( 0..$#{$minnumbers[$subpop]} ){
							if ( @{$validalleles[$subpop]}>=$minnumbers[$subpop][$minset] ){
								$sfsarray=$minset if $minset>$sfsarray;
								$valid=1; last SET;
								}
							}
						unless ($valid){ next GROUP }
						}
					if ($snp){
						for my $subpop ( @{$groups->[$metapop]} ){
							my @suballeles=( Misc::shuffle( @{$validalleles[$subpop]} ) )[ 0..$minnumbers[$subpop][$sfsarray]-1 ];
							for my $allele (@suballeles){ if ($allele ne $ancestral){ ++$nder } }
							}
						}
					if ($folded && $nder>$arraysizes[$sfsarray][$metapop] ){ $nder=2*$arraysizes[$sfsarray][$metapop]-$nder }
					++$self->{_SFS}[$metapop][$sfsarray][$nder];
					++$validpos[$metapop];
					}
				}
			}
		}
	print "Total number of assessed sites: $totalsites, validsites: ", join(',', @validpos), ", non-biallelic: $notbiallelic, undefined ancestral state: $undefanc, cpg filtered: $cpgfiltered\n";
	return;
	}


sub storeSFSbyCodon{
	my ($self, $regions, $refseq, $genotypes, $ref_missing, $groups, $mincov, $minind, $folded, $haplotize, $filtercpg)=@_;

	print "Generating SFS ...\n";
	my $ref_deg=Misc::getDegRef();
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum, @groupsizes);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	for my $ref_minset (@$minind){
		my @minsizes;
		for my $ref_group (@$groups){
			push @minsizes, Misc::sum(@{$ref_minset}[ @$ref_group ]);
			}
		push @groupsizes, \@minsizes;
		}	

	my ($undefanc, $codonpattern, $cpgfiltered, $totalsites)=(0, 0, 0, 0);
	my (@validpos, @arraysizes, @minnumbers);
	for my $minset ( 0..$#{$minind} ){
		for my $group ( 0..$#{$groups} ){
			if ($haplotize && $folded){ $arraysizes[$minset][$group]=int($groupsizes[$minset][$group]/2) }
			elsif ($haplotize || $folded){ $arraysizes[$minset][$group]=$groupsizes[$minset][$group] }
			else { $arraysizes[$minset][$group]=2*$groupsizes[$minset][$group] }
			}
		for my $pop ( 0..$#{$minind->[$minset]} ){
			if ($haplotize){ push @{$minnumbers[$pop]}, $minind->[$minset][$pop] }
			else { push @{$minnumbers[$pop]}, 2*$minind->[$minset][$pop] }
			}
		}

	unless ( exists $self->{_Selected} ){	# to prevent overwriting of previously stored counts.
		for my $minset ( 0..$#{$minind} ){
			for my $group ( 0..$#{$groups} ){
				$self->{_Selected}[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
				$self->{_Neutral}[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
				}
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $frame=$regions->{_Loci}{$chrom}[$locus]{'frame'};
			my $dir=$regions->{_Loci}{$chrom}[$locus]{'dir'};
			my $locuslength=$end-$start+1;
			my $shift=$frame;
			if ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			my $ref_seqstring=$refseq->randomAccess($chrom, $start+$shift, $end, 1);
			my $ref_cpgstring=Misc::getCpGString($ref_seqstring, $genotypes, $chrom, $start+$shift, $end);
#			print "$chrom:$start-$end, $frame, $dir, $shift\n$$ref_seqstring\n";

			CODON: for (my $codstart=$start+$shift; $codstart+2<=$end; $codstart+=3){
				$totalsites+=3;
				my $codoffset=$codstart-$start;
				my $codon=substr($$ref_seqstring, 0, 3, "");	# chew up reference string by next three bases.
				my $cpg_mask=substr($$ref_cpgstring, 0, 3, "");	# chew up cpgmask string by next three bases.
#				my $codon=substr($$ref_seqstring, $codoffset, 3);
				$codon=uc($codon);
				if ($codon=~/N/){
#					print "Codon $chrom:$codstart-", $codstart+3, " has a missing base ($codon)!\n";
					++$codonpattern;
					next CODON;
					}
				if ($dir eq '-'){ $codon=reverse($codon) }
				my $altcodon=$codon;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=$dir eq '-' ? $codstart+(2-$codonpos) : $codstart+$codonpos;
					if ( defined $genotypes->getValue($chrom, $pos) ){
						if ($poly){ warn "Warning: Codon $chrom:$codstart-", $codstart+3, " has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if ( $genotypes->NAlleles($chrom, $pos)>2 ){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						substr($altcodon, $codonpos, 1)=$genotypes->getBase($chrom, $pos, 1);
						$poly=1;
						}
					}
				$codon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $deg=$ref_deg->{$codon};
#				print "$chrom:$codstart - $codon\n";
				if (!defined $deg){ warn "Warning, degeneracy pattern for codon $codon at position $codstart, offset $codoffset is not defined! $frame, $dir\n"; ++$codonpattern; next CODON }
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at position $codstart, offset $codoffset is not defined!\n"; ++$codonpattern; next CODON }
				if ($deg eq 'stop' || $altdeg eq 'stop'){
					if ($dir eq '+' && $end-$codstart>2){ warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-", $codstart+2, "! CDS end at $end.\n" }
					elsif ($dir eq '-' && $codstart>$start+$shift){ warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-", $codstart+2, "! CDS end at ", $start+$shift, ".\n" }
					next CODON;
					}

				POS: for my $codonpos (0..2){
					my $degstate=substr($deg, $codonpos, 1);
					unless ($degstate eq '1' || $degstate eq '4'){ ++$codonpattern; next POS }
					if ($poly && $deg ne $altdeg){ next POS if (substr($altdeg, $codonpos, 1) ne $degstate) }
					my $pos=$codstart+$codonpos;
					my $offset=$pos-$start;
					if ($filtercpg){ unless ( substr($cpg_mask, $codonpos, 1) ){ ++$cpgfiltered; next POS } }
					my $snp=0; my $ancestral;
					$snp=1 if ( $poly && defined $genotypes->getValue($chrom, $pos) );

					if (!$folded && $snp){
						$ancestral=$genotypes->getAncestral($chrom, $pos);
						unless (defined $ancestral && $ancestral<=3){ ++$undefanc; next POS }
						}
					else { $ancestral=0 }

					GROUP: for my $metapop ( 0..$#{$groups} ){
						my $nder=0; my $sfsarray=0;
						my @validalleles;
						for my $subpop ( @{$groups->[$metapop]} ){
							IND: for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
								unless (defined $coverage && $coverage>=$mincov){ next IND }
								if ($snp){
									if ($haplotize){
										my $index=2*$ind+int(rand(2));
										my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
										unless (defined $allele && $allele ne '.'){ next IND }
										push @{$validalleles[$subpop]}, $allele;
										}
									else {
										my $index=2*$ind;
										my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
										unless (defined $allele1 && $allele1 ne '.'){ next IND }
										++$index;
										my $allele2=$genotypes->getValue($chrom, $pos, $subpop, $index);
										unless (defined $allele2 && $allele2 ne '.'){ next IND }
										push @{$validalleles[$subpop]}, ($allele1, $allele2);
										}
									}
								elsif ($haplotize){ push @{$validalleles[$subpop]}, 0 }
								else { push @{$validalleles[$subpop]}, (0, 0) }
								}
							my $valid=0;
							unless ( defined $validalleles[$subpop] ){ next GROUP }
							SET: for my $minset ( 0..$#{$minnumbers[$subpop]} ){
								if ( @{$validalleles[$subpop]}>=$minnumbers[$subpop][$minset] ){
									$sfsarray=$minset if $minset>$sfsarray;
									$valid=1; last SET;
									}
								}
							unless ($valid){ next GROUP }
							}

						if ($snp){
							for my $subpop ( @{$groups->[$metapop]} ){
								my @suballeles=( Misc::shuffle( @{$validalleles[$subpop]} ) )[ 0..$minnumbers[$subpop][$sfsarray]-1 ];
								for my $allele (@suballeles){ if ($allele ne $ancestral){ ++$nder } }
								}
							}

						if ($folded && $nder>$arraysizes[$sfsarray][$metapop] ){ $nder=2*$arraysizes[$sfsarray][$metapop]-$nder }
						if ($degstate eq '1'){ ++$self->{_Selected}[$metapop][$sfsarray][$nder] }
						elsif ($degstate eq '4'){ ++$self->{_Neutral}[$metapop][$sfsarray][$nder] }
						++$validpos[$metapop];
						}
					}
				}
			}
		}
	print "Total number of assessed sites: $totalsites, validsites: ", join(',', @validpos), ", invalid codon pattern: $codonpattern, undefined ancestral state: $undefanc, cpg filtered: $cpgfiltered\n";
	return;
	}

sub storeSFSbyCodon2{	# built after storeSFSbyCodonCat2
	my ($self, $regions, $ref_genotypes, $ref_missing, $groups, $mincov, $minind, $haplotize, $folded, $ref_recbins, $useoutgroups)=@_;

	print "Generating SFS ...\n";
	my $ref_deg=Misc::getDegRef();
	my $maxbin=( sort { $b<=>$a } @$ref_recbins )[0];	# take largest value of bin in selection.
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum, @groupsizes);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	for my $ref_minset (@$minind){
		my @minsizes;
		for my $ref_group (@$groups){
			push @minsizes, Misc::sum(@{$ref_minset}[ @$ref_group ]);
			}
		push @groupsizes, \@minsizes;
		}	

	my ($undefanc, $notbiallelic, $codonpattern, $totalcodons)=(0, 0, 0, 0);
	my (@validpostot, @nodata, @arraysizes, @minnumbers);

	for my $minset ( 0..$#{$minind} ){
		for my $group ( 0..$#{$groups} ){
			if ($haplotize && $folded){ $arraysizes[$minset][$group]=int($groupsizes[$minset][$group]/2) }
			elsif ($haplotize || $folded){ $arraysizes[$minset][$group]=$groupsizes[$minset][$group] }
			else { $arraysizes[$minset][$group]=2*$groupsizes[$minset][$group] }
			}
		for my $pop ( 0..$#{$minind->[$minset]} ){
			if ($haplotize){ push @{$minnumbers[$pop]}, $minind->[$minset][$pop] }
			else { push @{$minnumbers[$pop]}, 2*$minind->[$minset][$pop] }
			}
		}

	unless ( exists $self->{_Selected} || exists $self->{_Neutral} ){	# to prevent overwriting of previously stored counts.
		for my $group ( 0..$#{$groups} ){
			for my $minset ( 0..$#{$minind} ){
				$self->{_ValidPos}[$group][$minset]=0;
				$self->{_Selected}[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
				$self->{_Neutral}[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
				}
			}
		}


	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $bin=$regions->getValue($chrom, $locus, 'bin');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ancref=$folded ? $refseq : $regions->getValue($chrom, $locus, 'ancref');	# if folded SFS is selected, set ancestral reference sequence to reference sequence.
			my $ref_positions=$regions->getValue($chrom, $locus, 'positions');
			my @positions;

			my (@selectedsfs, @neutralsfs, @validpos);	# use temporary arrays per locus in order to be able to discard counts if stop codon is encountered.
			for my $minset ( 0..$#{$minind} ){
				for my $group ( 0..$#{ $groups } ){
					$selectedsfs[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
					$neutralsfs[$group][$minset]=[ (0) x ($arraysizes[$minset][$group]+1) ];
					$validpos[$group][$minset]=0;
					}
				}

			if ($dir eq '-'){
				@positions=reverse(@$ref_positions);	# copy and reverse the array of positions.
				$refseq=uc( reverse($refseq) );	# copy and reverse the reference string.
				$ancref=uc( reverse($ancref) );	# reverse the ancestral reference string.
				}
			else {
				@positions=@$ref_positions;
				$refseq=uc($refseq);
				$ancref=uc($ancref);
				}

			my $locuslength=length($refseq);
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			if ($shift){
				substr($refseq, 0, $shift, "");
				substr($ancref, 0, $shift, "");
				splice(@positions, 0, $shift);
				}
			if ($clip){
				substr($refseq, -$clip, $clip, "");
				substr($ancref, -$clip, $clip, "");
				splice(@positions, -$clip);
				}
			if (length($refseq) % 3){ warn "Reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (length($ancref) % 3){ warn "Ancestral reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (@positions % 3){ warn "Position array has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			print "$chrom:$start-$end, ", length($refseq), ", ", length($ancref), ", ", scalar(@positions), ", $shift, $clip\n";
			unless ( scalar(@positions)==length($refseq) && length($refseq)==length($ancref) ){ die "Lengths of reference sequence, ancestral reference and position array do not match!\n" }

			CODON: while ( length($refseq)>=3 ){
				++$totalcodons;
				my $codon_plusstrand=substr($refseq, 0, 3, "");	# chew up reference string by next three bases.
				my $anccodon_plusstrand=substr($ancref, 0, 3, "");	# chew up ancestral reference string by next three bases.
				my @triplet=splice(@positions, 0, 3);	# chew up position array by next three bases.
				my $codstart=($dir eq '-') ? $triplet[-1] : $triplet[0];
				my $codend=($dir eq '-') ? $triplet[0] : $triplet[-1];
				if ($codon_plusstrand=~/N/){	# if one of the codon positions is hardmasked.
					++$codonpattern;
					next CODON;
					}
				my $basecodon=$codon_plusstrand;
				my $altcodon=$codon_plusstrand;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=$triplet[$codonpos];
					if ( defined $ref_genotypes->getValue($chrom, $pos) ){
						my @alleles=$ref_genotypes->getAlleles($chrom, $pos);
						if ($poly && @alleles>1){ warn "Warning: Codon $chrom:$codstart-$codend has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if (@alleles>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						if (@alleles==1 && $alleles[0]>0){	# if all samples in the set have a fixed difference to the reference.
							my $altbase=$ref_genotypes->getBase($chrom, $pos, $alleles[0]);
							substr($basecodon, $codonpos, 1)=$altbase;
							substr($altcodon, $codonpos, 1)=$altbase;
							}
						elsif (@alleles==2){	# if the position is polymorphic in the sample.
							substr($basecodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[0]) unless ($alleles[0]==0);	# only change base codon if first allele nonreference.
							substr($altcodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[1]);
							$poly=1;
							}
						}
					}

				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $deg=$ref_deg->{$basecodon};
				if (!defined $deg){ warn "Warning, degeneracy pattern for codon $basecodon at position at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				if ($deg eq 'stop' || $altdeg eq 'stop'){
					if (@positions){
						warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-$codend! CDS end at $end. Skipping transcript.\n";
						$regions->setExcluded($chrom, $locus, 1);
						next LOCUS;	# discard the whole transcript.
						}
					else { next CODON }	# regular stop codon, proceed to addition of counts.
					}

				POS: for my $codonpos (0..2){
					my $degstate=substr($deg, $codonpos, 1);
					unless ($degstate eq '1' || $degstate eq '4'){ next POS }
					if ($poly && $deg ne $altdeg){ next POS if (substr($altdeg, $codonpos, 1) ne $degstate) }
					my $pos=$triplet[$codonpos];
					my $ancrefbase=substr($anccodon_plusstrand, $codonpos, 1);
					unless (defined $ancrefbase && $ancrefbase ne 'N'){ ++$undefanc; next POS }

					my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
					my $ancestral;
					if ($snp){
						my @derived;
						if (!$folded && $useoutgroups){
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos);
							my $ancbase=$ref_genotypes->getBase($chrom, $pos, $ancestral);
							unless ($ancbase eq $ancrefbase){
								warn "Base from ancestral reference $ancrefbase at $chrom: $pos does not correspond to ancestral state from outgroups $ancbase!\n";
								++$undefanc;
								next POS;
								}
							}
						else {
							my $ancallele=$folded ? 0 : $ref_genotypes->getAlleleforBase($chrom, $pos, $ancrefbase);
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos, $ancallele);
							}
						unless (defined $ancestral){ ++$undefanc; next POS }
						if (@derived==0){ $snp=0 }	# position is monomorphic reference.
						elsif (@derived>1){ ++$notbiallelic; next POS }	# maximum one derived allele permitted.
						}

					GROUP: for my $metapop ( 0..$#{$groups} ){
						my @validalleles;
						my $sfsarray=0;
						for my $subpop ( @{$groups->[$metapop]} ){
							IND: for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValuebyPos($chrom, $pos, $allind);
								unless (defined $coverage && $coverage>=$mincov){ next IND }
								if ($snp && $haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @{$validalleles[$subpop]}, $allele;
									}
								elsif ($snp){
									my $index=2*$ind;
									my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @{$validalleles[$subpop]}, ($allele1, $allele2);
									}
								elsif ($haplotize) { push @{$validalleles[$subpop]}, 0 }
								else { push @{$validalleles[$subpop]}, (0, 0) }
								}
							unless ( defined $validalleles[$subpop] ){ ++$nodata[$metapop]; next GROUP }
							my $valid=0;
							SET: for my $minset ( 0..$#{$minnumbers[$subpop]} ){
								if ( @{$validalleles[$subpop]}>=$minnumbers[$subpop][$minset] ){
									$sfsarray=$minset if $minset>$sfsarray;
									$valid=1; last SET;
									}
								}
							unless ($valid){ ++$nodata[$metapop]; next GROUP }
							}

						my $nder=0;
						if ($snp){
							for my $subpop ( @{$groups->[$metapop]} ){
								my $ref_suballeles=Misc::samplewithoutRep($validalleles[$subpop], $minnumbers[$subpop][$sfsarray]);
								for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
								}
							}

						if ($folded && $nder>$arraysizes[$sfsarray][$metapop] ){ $nder=2*$arraysizes[$sfsarray][$metapop]-$nder }
						if ($degstate eq '1'){ ++$selectedsfs[$metapop][$sfsarray][$nder]	}
						elsif ($degstate eq '4'){ ++$neutralsfs[$metapop][$sfsarray][$nder] }
						++$validpos[$metapop][$sfsarray];
						++$validpostot[$metapop];
						}
					}
				}

			for my $group ( 0..$#{ $groups } ){
				for my $minset ( 0..$#{$minind} ){
					$self->{_ValidPos}[$group][$minset]+=$validpos[$group][$minset];
					$self->{_Selected}[$group][$minset][$_]+=$selectedsfs[$group][$minset][$_] foreach (0..$arraysizes[$group]);
					$self->{_Neutral}[$group][$minset][$_]+=$neutralsfs[$group][$minset][$_] foreach (0..$arraysizes[$group]);
					}
				}
			}
		}
	print "Total number of assessed codons: $totalcodons, invalid codon pattern: $codonpattern, validsites: ", join(',', @validpostot), ", undefined ancestral state: $undefanc, ",
		"not biallelic: $notbiallelic, no data: ", join(',', @nodata), "\n";
	return;
	}


sub storeSFSbyCodonCat{
	my ($self, $regions, $refseq, $ancref, $genotypes, $ref_missing, $groups, $mincov, $minind, $haplotize, $ref_recbins)=@_;

	print "Generating SFS ...\n";
	my $ref_deg=Misc::getDegRef();
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	my ($undefanc, $notbiallelic, $codonpattern, $notpoly, $totalsites)=(0, 0, 0, 0, 0);
	my (@validpos, @nodata, @arraysizes, @minnumbers);
	for my $ref_group (@$groups){
		if ($haplotize){ push @arraysizes, Misc::sum(@$minind[ @$ref_group ]) }
		else { push @arraysizes, 2*Misc::sum(@$minind[ @$ref_group ]) }
		}

	for my $pop ( 0..$#{ $minind } ){
		if ($haplotize){ push @minnumbers, $minind->[$pop] }
		else { push @minnumbers, 2*$minind->[$pop] }
		}

	unless ( exists $self->{_Selected} ){	# to prevent overwriting of previously stored counts.
		for my $group ( 0..$#{$groups} ){
			$self->{_Selected}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
			$self->{_Neutral}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
			for my $bin ( 0..$#{ $ref_recbins } ){
				$self->{_SelectedSegs}[$group][$bin]=[ (0) x 5 ];
				$self->{_NeutralSegs}[$group][$bin]=[ (0) x 5 ];
				}
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $bin=$regions->getValue($chrom, $locus, 'bin');
			unless (defined $bin){ warn "Bin is not defined for $chrom: $start-$end, skipping locus!\n"; next LOCUS }
			my $locuslength=$end-$start+1;
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			my $refseq_string=$refseq->randomAccess($chrom, $start+$shift, $end-$clip, 1);
			my $ancref_string=$ancref->randomAccess($chrom, $start+$shift, $end-$clip, 1);
			if (length($$refseq_string) % 3){ warn "$chrom:", $start+$shift, "-", $end-$clip, " is not a multiple of 3!\n" }

			my @codstarts;
			for (my $codstart=$start+$shift; $codstart+2<=$end; $codstart+=3){ push(@codstarts, $codstart) }
			if ($dir eq '-'){
				@codstarts=reverse(@codstarts);	# reverse the array of codon starting positions.
				$$refseq_string=reverse($$refseq_string);	# reverse the reference string.
				$$ancref_string=reverse($$ancref_string);	# reverse the ancestral reference string.
				}

			CODON: for my $codstart(@codstarts){	# loop through the array of codon start positions
				$totalsites+=3;
				my $codoffset=$codstart-$start;
				my $refcodon=substr($$refseq_string, 0, 3, "");	# chew up reference string by next three bases.
				my $anccodon=substr($$ancref_string, 0, 3, "");	# chew up ancestral reference string by next three bases.
				$refcodon=uc($refcodon);
				$anccodon=uc($anccodon);
				if ($refcodon=~/N/){	# if one of the codon positions is hardmasked.
					++$codonpattern;
					next CODON;
					}
				my $basecodon=$refcodon;
				my $altcodon=$refcodon;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=($dir eq '-') ? $codstart+(2-$codonpos) : $codstart+$codonpos;	# to make sure that reverse codons are processed backwards.
					if ( defined $genotypes->getValue($chrom, $pos) ){
						my @alleles=$genotypes->getAlleles($chrom, $pos);
						if ($poly && @alleles>1){ warn "Warning: Codon $chrom:$codstart-", $codstart+3, " has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if (@alleles>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						if (@alleles==1 && $alleles[0]>0){	# if all samples in the set have a fixed difference to the reference.
							my $altbase=$genotypes->getBase($chrom, $pos, $alleles[0]);
							substr($basecodon, $codonpos, 1)=$altbase;
							substr($altcodon, $codonpos, 1)=$altbase;
							}
						elsif (@alleles==2){	# if the position is polymorphic in the sample.
							substr($basecodon, $codonpos, 1)=$genotypes->getBase($chrom, $pos, $alleles[0]) unless ($alleles[0]==0);	# only change base codon if first allele nonreference.
							substr($altcodon, $codonpos, 1)=$genotypes->getBase($chrom, $pos, $alleles[1]);
							$poly=1;
							}
						}
					}
				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $deg=$ref_deg->{$basecodon};
				if (!defined $deg){ warn "Warning, degeneracy pattern for codon $basecodon at position $codstart, offset $codoffset is not defined! $frame, $dir\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at position $codstart, offset $codoffset is not defined!\n"; ++$codonpattern; next CODON }
				if ($deg eq 'stop' || $altdeg eq 'stop'){
					unless ($codstart==$codstarts[-1]){ warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-", $codstart+2, "! CDS end at $codstarts[-1].\n" }
					next LOCUS;
					}

				POS: for my $codonpos (0..2){
					my $degstate=substr($deg, $codonpos, 1);
					unless ($degstate eq '1' || $degstate eq '4'){ ++$codonpattern; next POS }
					if ($poly && $deg ne $altdeg){ next POS if (substr($altdeg, $codonpos, 1) ne $degstate) }
					my $pos=($dir eq '-') ? $codstart+(2-$codonpos) : $codstart+$codonpos;
					my $offset=$pos-$start;
					unless (defined $genotypes->getValue($chrom, $pos) && $genotypes->NAlleles($chrom, $pos)>1){ ++$notpoly; next POS };
					my ($ancestral, @derived)=$genotypes->getPolarizedAlleles($chrom, $pos);
					unless (defined $ancestral){ ++$undefanc; next POS }
					if (@derived>1){ ++$undefanc; next POS }	# maximum one derived allele permitted.
					my $ancbase=$genotypes->getBase($chrom, $pos, $ancestral);
					my $derbase=$genotypes->getBase($chrom, $pos, $derived[0]);
					my $ancrefbase=uc( substr($anccodon, $codonpos, 1) );
					unless ($ancbase eq $ancrefbase){
						warn "Base from ancestral reference $ancrefbase at $chrom: $pos does not correspond to ancestral state from outgroups $ancbase!\n";
						}
					my $cat=Misc::getWeakStrong($ancbase, $derbase);

					GROUP: for my $metapop ( 0..$#{$groups} ){
						my $nder=0; my $sfsarray=0;
						my @validalleles;
						for my $subpop ( @{$groups->[$metapop]} ){
							IND: for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
								unless (defined $coverage && $coverage>=$mincov){ next IND }
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @{$validalleles[$subpop]}, $allele;
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @{$validalleles[$subpop]}, ($allele1, $allele2);
									}
								}
							unless ( defined $validalleles[$subpop] && @{ $validalleles[$subpop] }>=$minnumbers[$subpop] ){ ++$nodata[$metapop]; next GROUP }
							}

						for my $subpop ( @{$groups->[$metapop]} ){
							my $ref_suballeles=Misc::samplewithoutRep($validalleles[$subpop], $minnumbers[$subpop]);
							for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
							}

						if ($degstate eq '1'){
							++$self->{_SelectedSFS}[$metapop][$cat][$nder];
							++$self->{_SelectedSegs}[$metapop][$bin][$cat] if ($nder>0 && $nder<$arraysizes[$metapop]);
							}
						elsif ($degstate eq '4'){
							++$self->{_NeutralSFS}[$metapop][$cat][$nder];
							++$self->{_NeutralSegs}[$metapop][$bin][$cat] if ($nder>0 && $nder<$arraysizes[$metapop]);
							}
						++$validpos[$metapop];
						if ($ancrefbase eq 'A' || $ancrefbase eq 'T'){ ++$self->{_Basecounts}[$metapop][$bin]{'AT'} }
						elsif ($ancrefbase eq 'G' || $ancrefbase eq 'C'){ ++$self->{_Basecounts}[$metapop][$bin]{'GC'} }
						else { ++$self->{_Basecounts}[$metapop][$bin]{'N'} }
						}
					}
				}
			}
		}
	print "Total number of assessed sites: $totalsites, validsites: ", join(',', @validpos), ", invalid codon pattern: $codonpattern, undefined ancestral state: $undefanc,
		not polymorphic: $notpoly, no data: ", join(',', @nodata), "\n";
	return;
	}


sub storeSFSbyCodonCat2{
	my ($self, $regions, $ref_genotypes, $ref_missing, $groups, $mincov, $minind, $haplotize, $ref_recbins, $useoutgroups, $printpergene, $outfile)=@_;

	my $output;
	if ($printpergene){
		die "Per gene output selected, but no outfile specified!\n" unless (defined $outfile);
		if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }
#		print $output "GeneID\tTranscriptID";
#		print $output "\tValidsites\tSel_W->W\tSel_W->S\tSel_S->W\tSel_S->S\tSel_undef\tSel_AT\tSel_GC\tSel_N\tSel_\%GC\tNeut_W->W\tNeut_W->S\tNeut_S->W\tNeut_S->S\tNeut_undef\tNeut_AT\tNeut_GC\tNeut_N\tNeut_\%GC" foreach (0..$#{$groups});
#		print $output "\n";
		}

	print "Generating SFS ...\n";
	my $ref_deg=Misc::getDegRef();
	my $maxbin=( sort { $b<=>$a } @$ref_recbins )[0];	# take largest value of bin in selection.
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	my ($undefanc, $notbiallelic, $codonpattern, $totalcodons)=(0, 0, 0, 0);
	my (@validpostot, @nodata, @arraysizes, @minnumbers);
	for my $ref_group (@$groups){
		if ($haplotize){ push @arraysizes, Misc::sum(@$minind[ @$ref_group ]) }
		else { push @arraysizes, 2*Misc::sum(@$minind[ @$ref_group ]) }
		}

	for my $pop ( 0..$#{ $minind } ){
		if ($haplotize){ push @minnumbers, $minind->[$pop] }
		else { push @minnumbers, 2*$minind->[$pop] }
		}

	unless ( exists $self->{_SelectedSFS} ){	# to prevent overwriting of previously stored counts.
		for my $group ( 0..$#{ $groups } ){
			$self->{_ValidPos}[$group]=[ (0) x ($maxbin+1) ];
#			$self->{_SelectedSFS}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
#			$self->{_NeutralSFS}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
			for my $bin (0..$maxbin){
				$self->{_SelectedSegs}[$group][$bin]=[ (0) x 5 ];
				$self->{_NeutralSegs}[$group][$bin]=[ (0) x 5 ];
				$self->{_Basecounts}[$group][$bin]={ 'selAT'=>0, 'selGC'=>0, 'selN'=>0, 'neutAT'=>0, 'neutGC'=>0, 'neutN'=>0 };
				$self->{_SelectedSFS}[$group][$bin][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$self->{_NeutralSFS}[$group][$bin][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				}
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $bin=$regions->getValue($chrom, $locus, 'bin');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ancref=$regions->getValue($chrom, $locus, 'ancref');
			my $ref_positions=$regions->getValue($chrom, $locus, 'positions');
			unless (defined $bin && $bin<=$maxbin){ warn "$geneid, $transcriptid has non-valid recombination bin $bin, skipping transcript!\n"; next LOCUS }
			my @positions;

			my (@selectedsfs, @neutralsfs, @selectedsegs, @neutralsegs, @basecounts, @validpos);	# use temporary arrays per locus in order to be able to discard counts if stop codon is encountered.
			for my $group ( 0..$#{ $groups } ){
				$selectedsfs[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$neutralsfs[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$selectedsegs[$group]=[ (0) x 5 ];
				$neutralsegs[$group]=[ (0) x 5 ];
				$basecounts[$group]={ 'selAT'=>0, 'selGC'=>0, 'selN'=>0, 'neutAT'=>0, 'neutGC'=>0, 'neutN'=>0 };
				$validpos[$group]=0;
				}

			if ($dir eq '-'){
				@positions=reverse(@$ref_positions);	# copy and reverse the array of positions.
				$refseq=uc( reverse($refseq) );	# copy and reverse the reference string.
				$ancref=uc( reverse($ancref) );	# reverse the ancestral reference string.
				}
			else {
				@positions=@$ref_positions;
				$refseq=uc($refseq);
				$ancref=uc($ancref);
				}

			my $locuslength=length($refseq);
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
#			print "$chrom:$start-$end, ", length($refseq), ", ", length($ancref), ", ", scalar(@positions), ", $shift, $clip\n";
			if ($shift){
				substr($refseq, 0, $shift, "");
				substr($ancref, 0, $shift, "");
				splice(@positions, 0, $shift);
				}
			if ($clip){
				substr($refseq, -$clip, $clip, "");
				substr($ancref, -$clip, $clip, "");
				splice(@positions, -$clip);
				}
			if (length($refseq) % 3){ warn "Reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (length($ancref) % 3){ warn "Ancestral reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (@positions % 3){ warn "Position array has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			print "$chrom:$start-$end, ", length($refseq), ", ", length($ancref), ", ", scalar(@positions), ", $shift, $clip\n";
			unless ( scalar(@positions)==length($refseq) && length($refseq)==length($ancref) ){ die "Lengths of reference sequence, ancestral reference and position array do not match!\n" }

			CODON: while ( length($refseq)>=3 ){
				++$totalcodons;
				my $codon_plusstrand=substr($refseq, 0, 3, "");	# chew up reference string by next three bases.
				my $anccodon_plusstrand=substr($ancref, 0, 3, "");	# chew up ancestral reference string by next three bases.
				my @triplet=splice(@positions, 0, 3);	# chew up position array by next three bases.
				my $codstart=($dir eq '-') ? $triplet[-1] : $triplet[0];
				my $codend=($dir eq '-') ? $triplet[0] : $triplet[-1];
				if ($codon_plusstrand=~/N/){	# if one of the codon positions is hardmasked.
					++$codonpattern;
					next CODON;
					}
				my $basecodon=$codon_plusstrand;
				my $altcodon=$codon_plusstrand;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=$triplet[$codonpos];
					if ( defined $ref_genotypes->getValue($chrom, $pos) ){
						my @alleles=$ref_genotypes->getAlleles($chrom, $pos);
						if ($poly && @alleles>1){ warn "Warning: Codon $chrom:$codstart-$codend has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if (@alleles>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						if (@alleles==1 && $alleles[0]>0){	# if all samples in the set have a fixed difference to the reference.
							my $altbase=$ref_genotypes->getBase($chrom, $pos, $alleles[0]);
							substr($basecodon, $codonpos, 1)=$altbase;
							substr($altcodon, $codonpos, 1)=$altbase;
							}
						elsif (@alleles==2){	# if the position is polymorphic in the sample.
							substr($basecodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[0]) unless ($alleles[0]==0);	# only change base codon if first allele nonreference.
							substr($altcodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[1]);
							$poly=1;
							}
						}
					}

				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $deg=$ref_deg->{$basecodon};
				if (!defined $deg){ warn "Warning, degeneracy pattern for codon $basecodon at position at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				if ($deg eq 'stop' || $altdeg eq 'stop'){
					if (@positions){
						warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-$codend! CDS end at $end. Skipping transcript.\n";
						$regions->setExcluded($chrom, $locus, 1);
						next LOCUS;	# discard the whole transcript.
						}
					else { next CODON }	# regular stop codon, proceed to addition of counts.
					}

				POS: for my $codonpos (0..2){
					my $degstate=substr($deg, $codonpos, 1);
					unless ($degstate eq '1' || $degstate eq '4'){ next POS }
					if ($poly && $deg ne $altdeg){ next POS if (substr($altdeg, $codonpos, 1) ne $degstate) }
					my $pos=$triplet[$codonpos];
					my $ancrefbase=substr($anccodon_plusstrand, $codonpos, 1);
					unless (defined $ancrefbase && $ancrefbase ne 'N'){ ++$undefanc; next POS }

					my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
					my $ancestral;
					my $cat=4;	# 4 is undefined category state.
					if ($snp){
						my @derived;
						if ($useoutgroups){
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos);
							my $ancbase=$ref_genotypes->getBase($chrom, $pos, $ancestral);
							unless ($ancbase eq $ancrefbase){
								warn "Base from ancestral reference $ancrefbase at $chrom: $pos does not correspond to ancestral state from outgroups $ancbase!\n";
								++$undefanc;
								next POS;
								}
							}
						else {
							my $ancallele=$ref_genotypes->getAlleleforBase($chrom, $pos, $ancrefbase);
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos, $ancallele);
							}
						unless (defined $ancestral){ ++$undefanc; next POS }
						if (@derived==0){ $snp=0 }	# position is monomorphic reference.
						elsif (@derived>1){ ++$notbiallelic; next POS }	# maximum one derived allele permitted.
						else {
							my $derbase=$ref_genotypes->getBase($chrom, $pos, $derived[0]);
							$cat=Misc::getWeakStrong($ancrefbase, $derbase);
							}
						}

					GROUP: for my $metapop ( 0..$#{$groups} ){
						my @validalleles;
						for my $subpop ( @{$groups->[$metapop]} ){
							IND: for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValuebyPos($chrom, $pos, $allind);
								unless (defined $coverage && $coverage>=$mincov){ next IND }
								if ($snp && $haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.'){ next IND }
									push @{$validalleles[$subpop]}, $allele;
									}
								elsif ($snp){
									my $index=2*$ind;
									my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.'){ next IND }
									my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.'){ next IND }
									push @{$validalleles[$subpop]}, ($allele1, $allele2);
									}
								elsif ($haplotize) { push @{$validalleles[$subpop]}, 0 }
								else { push @{$validalleles[$subpop]}, (0, 0) }
								}
							unless ( defined $validalleles[$subpop] && @{ $validalleles[$subpop] }>=$minnumbers[$subpop] ){ ++$nodata[$metapop]; next GROUP }
							}

						my $nder=0;
						if ($snp){
							for my $subpop ( @{$groups->[$metapop]} ){
								my $ref_suballeles=Misc::samplewithoutRep($validalleles[$subpop], $minnumbers[$subpop]);
								for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
								}
							}

						if ($nder>0){
							if ($degstate eq '1'){
								++$selectedsfs[$metapop][$cat][$nder];
								++$selectedsegs[$metapop][$cat] if ($nder<$arraysizes[$metapop]);
								}
							elsif ($degstate eq '4'){
								++$neutralsfs[$metapop][$cat][$nder];
								++$neutralsegs[$metapop][$cat] if ($nder<$arraysizes[$metapop]);
								}
							}

						++$validpostot[$metapop];
						++$validpos[$metapop];
						if ($degstate eq '1'){
							if ($ancrefbase eq 'A' || $ancrefbase eq 'T'){ ++$basecounts[$metapop]{'selAT'} }
							elsif ($ancrefbase eq 'G' || $ancrefbase eq 'C'){ ++$basecounts[$metapop]{'selGC'} }
							else { ++$basecounts[$metapop]{'selN'} }
							}
						elsif ($degstate eq '4'){
							if ($ancrefbase eq 'A' || $ancrefbase eq 'T'){ ++$basecounts[$metapop]{'neutAT'} }
							elsif ($ancrefbase eq 'G' || $ancrefbase eq 'C'){ ++$basecounts[$metapop]{'neutGC'} }
							else { ++$basecounts[$metapop]{'neutN'} }
							}
						}
					}
				}

			if ($printpergene){
				print $output "$geneid\t$transcriptid";
				for my $group (0..$#{$groups}){
					print $output "\t$validpos[$group]";
					for my $cat (0..4){
#						print $output "\t$selectedsegs[$group][$cat]";
						print $output "\t", join(',', @{ $selectedsfs[$group][$cat] });
						}
					print $output "\t$basecounts[$group]{'selAT'}\t$basecounts[$group]{'selGC'}\t$basecounts[$group]{'selN'}";
					my $validbases=$basecounts[$group]{'selAT'}+$basecounts[$group]{'selGC'};
					my $gccont=$validbases>0 ? $basecounts[$group]{'selGC'}/$validbases : 0;
					print $output "\t", sprintf("%.2f", $gccont*100);
					for my $cat (0..4){
#						print $output "\t$neutralsegs[$group][$cat]";
						print $output "\t", join(',', @{ $neutralsfs[$group][$cat] });
						}
					print $output "\t$basecounts[$group]{'neutAT'}\t$basecounts[$group]{'neutGC'}\t$basecounts[$group]{'neutN'}";
					$validbases=$basecounts[$group]{'neutAT'}+$basecounts[$group]{'neutGC'};
					$gccont=$validbases>0 ? $basecounts[$group]{'neutGC'}/$validbases : 0;
					print $output "\t", sprintf("%.2f", $gccont*100);
					}
				print $output "\n";
				}

			for my $group ( 0..$#{ $groups } ){
#				for my $cat (0..4){
#					$self->{_SelectedSFS}[$group][$cat][$_]+=$selectedsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
#					$self->{_NeutralSFS}[$group][$cat][$_]+=$neutralsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
#					}
				$self->{_ValidPos}[$group][$bin]+=$validpos[$group];
				$self->{_Basecounts}[$group][$bin]{$_}+=$basecounts[$group]{$_} foreach ('selAT', 'selGC', 'selN', 'neutAT', 'neutGC', 'neutN');
				for my $cat (0..4){
					$self->{_SelectedSegs}[$group][$bin][$cat]+=$selectedsegs[$group][$cat];
					$self->{_NeutralSegs}[$group][$bin][$cat]+=$neutralsegs[$group][$cat];
					$self->{_SelectedSFS}[$group][$bin][$cat][$_]+=$selectedsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
					$self->{_NeutralSFS}[$group][$bin][$cat][$_]+=$neutralsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
					}
				}
			}
		}

	if ($printpergene && $outfile ne "-" && $outfile ne "STDOUT"){ close $output }

	print "Total number of assessed codons: $totalcodons, invalid codon pattern: $codonpattern, validsites: ", join(',', @validpostot), ", undefined ancestral state: $undefanc, ",
		"not biallelic: $notbiallelic, no data: ", join(',', @nodata), "\n";
	return;
	}

sub storeSFSbyCodonCatZchrom{	# storeSFSbyCodonCat2 adapted for Z-chromosomal data
	my ($self, $regions, $ref_genotypes, $ref_missing, $groups, $ref_mincov, $ref_minalleles, $ref_sexes, $haplotize, $ref_recbins, $useoutgroups, $printpergene, $outfile)=@_;

	my $output;
	if ($printpergene){
		die "Per gene output selected, but no outfile specified!\n" unless (defined $outfile);
		if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }
#		print $output "GeneID\tTranscriptID";
#		print $output "\tValidsites\tSel_W->W\tSel_W->S\tSel_S->W\tSel_S->S\tSel_undef\tSel_AT\tSel_GC\tSel_N\tSel_\%GC\tNeut_W->W\tNeut_W->S\tNeut_S->W\tNeut_S->S\tNeut_undef\tNeut_AT\tNeut_GC\tNeut_N\tNeut_\%GC" foreach (0..$#{$groups});
#		print $output "\n";
		}

	print "Generating SFS ...\n";
	my $ref_deg=Misc::getDegRef();
	my $maxbin=( sort { $b<=>$a } @$ref_recbins )[0];	# take largest value of bin in selection.
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	my ($undefanc, $notbiallelic, $codonpattern, $totalcodons)=(0, 0, 0, 0);
	my (@validpostot, @nodata, @arraysizes, @minnumbers);
	for my $ref_group (@$groups){ push @arraysizes, Misc::sum(@$ref_minalleles[ @$ref_group ]) }
	for my $pop ( 0..$#{ $ref_minalleles } ){ push @minnumbers, $ref_minalleles->[$pop] }

	unless ( exists $self->{_SelectedSFS} ){	# to prevent overwriting of previously stored counts.
		for my $group ( 0..$#{ $groups } ){
			$self->{_ValidPos}[$group]=[ (0) x ($maxbin+1) ];
#			$self->{_SelectedSFS}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
#			$self->{_NeutralSFS}[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
			for my $bin (0..$maxbin){
				$self->{_SelectedSegs}[$group][$bin]=[ (0) x 5 ];
				$self->{_NeutralSegs}[$group][$bin]=[ (0) x 5 ];
				$self->{_Basecounts}[$group][$bin]={ 'selAT'=>0, 'selGC'=>0, 'selN'=>0, 'neutAT'=>0, 'neutGC'=>0, 'neutN'=>0 };
				$self->{_SelectedSFS}[$group][$bin][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$self->{_NeutralSFS}[$group][$bin][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				}
			}
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $bin=$regions->getValue($chrom, $locus, 'bin');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ancref=$regions->getValue($chrom, $locus, 'ancref');
			my $ref_positions=$regions->getValue($chrom, $locus, 'positions');
			unless (defined $bin && $bin<=$maxbin){ warn "$geneid, $transcriptid has non-valid recombination bin $bin, skipping transcript!\n"; next LOCUS }
			my @positions;

			my (@selectedsfs, @neutralsfs, @selectedsegs, @neutralsegs, @basecounts, @validpos);	# use temporary arrays per locus in order to be able to discard counts if stop codon is encountered.
			for my $group ( 0..$#{ $groups } ){
				$selectedsfs[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$neutralsfs[$group][$_]=[ (0) x ($arraysizes[$group]+1) ] foreach (0..4);
				$selectedsegs[$group]=[ (0) x 5 ];
				$neutralsegs[$group]=[ (0) x 5 ];
				$basecounts[$group]={ 'selAT'=>0, 'selGC'=>0, 'selN'=>0, 'neutAT'=>0, 'neutGC'=>0, 'neutN'=>0 };
				$validpos[$group]=0;
				}

			if ($dir eq '-'){
				@positions=reverse(@$ref_positions);	# copy and reverse the array of positions.
				$refseq=uc( reverse($refseq) );	# copy and reverse the reference string.
				$ancref=uc( reverse($ancref) );	# reverse the ancestral reference string.
				}
			else {
				@positions=@$ref_positions;
				$refseq=uc($refseq);
				$ancref=uc($ancref);
				}

			my $locuslength=length($refseq);
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
#			print "$chrom:$start-$end, ", length($refseq), ", ", length($ancref), ", ", scalar(@positions), ", $shift, $clip\n";
			if ($shift){
				substr($refseq, 0, $shift, "");
				substr($ancref, 0, $shift, "");
				splice(@positions, 0, $shift);
				}
			if ($clip){
				substr($refseq, -$clip, $clip, "");
				substr($ancref, -$clip, $clip, "");
				splice(@positions, -$clip);
				}
			if (length($refseq) % 3){ warn "Reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (length($ancref) % 3){ warn "Ancestral reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (@positions % 3){ warn "Position array has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			print "$chrom:$start-$end, ", length($refseq), ", ", length($ancref), ", ", scalar(@positions), ", $shift, $clip\n";
			unless ( scalar(@positions)==length($refseq) && length($refseq)==length($ancref) ){ die "Lengths of reference sequence, ancestral reference and position array do not match!\n" }

			CODON: while ( length($refseq)>=3 ){
				++$totalcodons;
				my $codon_plusstrand=substr($refseq, 0, 3, "");	# chew up reference string by next three bases.
				my $anccodon_plusstrand=substr($ancref, 0, 3, "");	# chew up ancestral reference string by next three bases.
				my @triplet=splice(@positions, 0, 3);	# chew up position array by next three bases.
				my $codstart=($dir eq '-') ? $triplet[-1] : $triplet[0];
				my $codend=($dir eq '-') ? $triplet[0] : $triplet[-1];
				if ($codon_plusstrand=~/N/){	# if one of the codon positions is hardmasked.
					++$codonpattern;
					next CODON;
					}
				my $basecodon=$codon_plusstrand;
				my $altcodon=$codon_plusstrand;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=$triplet[$codonpos];
					if ( defined $ref_genotypes->getValue($chrom, $pos) ){
						my @alleles=$ref_genotypes->getAlleles($chrom, $pos);
						if ($poly && @alleles>1){ warn "Warning: Codon $chrom:$codstart-$codend has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if (@alleles>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						if (@alleles==1 && $alleles[0]>0){	# if all samples in the set have a fixed difference to the reference.
							my $altbase=$ref_genotypes->getBase($chrom, $pos, $alleles[0]);
							substr($basecodon, $codonpos, 1)=$altbase;
							substr($altcodon, $codonpos, 1)=$altbase;
							}
						elsif (@alleles==2){	# if the position is polymorphic in the sample.
							substr($basecodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[0]) unless ($alleles[0]==0);	# only change base codon if first allele nonreference.
							substr($altcodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, $alleles[1]);
							$poly=1;
							}
						}
					}

				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $deg=$ref_deg->{$basecodon};
				if (!defined $deg){ warn "Warning, degeneracy pattern for codon $basecodon at position at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				if ($deg eq 'stop' || $altdeg eq 'stop'){
					if (@positions){
						warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-$codend! CDS end at $end. Skipping transcript.\n";
						$regions->setExcluded($chrom, $locus, 1);
						next LOCUS;	# discard the whole transcript.
						}
					else { next CODON }	# regular stop codon, proceed to addition of counts.
					}

				POS: for my $codonpos (0..2){
					my $degstate=substr($deg, $codonpos, 1);
					unless ($degstate eq '1' || $degstate eq '4'){ next POS }
					if ($poly && $deg ne $altdeg){ next POS if (substr($altdeg, $codonpos, 1) ne $degstate) }
					my $pos=$triplet[$codonpos];
					my $ancrefbase=substr($anccodon_plusstrand, $codonpos, 1);
					unless (defined $ancrefbase && $ancrefbase ne 'N'){ ++$undefanc; next POS }

					my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
					my $ancestral;
					my $cat=4;	# 4 is undefined category state.
					if ($snp){
						my @derived;
						if ($useoutgroups){
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos);
							my $ancbase=$ref_genotypes->getBase($chrom, $pos, $ancestral);
							unless ($ancbase eq $ancrefbase){
								warn "Base from ancestral reference $ancrefbase at $chrom: $pos does not correspond to ancestral state from outgroups $ancbase!\n";
								++$undefanc;
								next POS;
								}
							}
						else {
							my $ancallele=$ref_genotypes->getAlleleforBase($chrom, $pos, $ancrefbase);
							($ancestral, @derived)=$ref_genotypes->getPolarizedAlleles($chrom, $pos, $ancallele);
							}
						unless (defined $ancestral){ ++$undefanc; next POS }
						if (@derived==0){ $snp=0 }	# position is monomorphic reference.
						elsif (@derived>1){ ++$notbiallelic; next POS }	# maximum one derived allele permitted.
						else {
							my $derbase=$ref_genotypes->getBase($chrom, $pos, $derived[0]);
							$cat=Misc::getWeakStrong($ancrefbase, $derbase);
							}
						}

					GROUP: for my $metapop ( 0..$#{$groups} ){
						my @validalleles;
						for my $subpop ( @{$groups->[$metapop]} ){
							IND: for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValuebyPos($chrom, $pos, $allind);
								my $ishemi=$ref_sexes->[$allind];
								my $mincov=$ishemi ? $ref_mincov->[0] : $ref_mincov->[1];
								unless (defined $coverage && $coverage>=$mincov){ next IND }
								if ($snp){
									if ($ishemi){
										my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, 2*$ind+1);
										unless ($allele2 eq '-' || $allele2 eq '.'){ warn "$chrom:$pos, population $subpop, individual $ind is supposed to be hemizygous, but got an allele call for the second allele: $allele2!\n"; next IND }
										my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, 2*$ind);
										unless (defined $allele1 && $allele1 ne '.'){ next IND }
										push @{$validalleles[$subpop]}, $allele1;
										}
									elsif ($haplotize){
										my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, 2*$ind+int(rand(2)) );
										unless (defined $allele1 && $allele1 ne '.'){ next IND }
										push @{$validalleles[$subpop]}, $allele1;
										}
									else {
										my $index=2*$ind;
										my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
										unless (defined $allele1 && $allele1 ne '.'){ next IND }
										my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
										unless (defined $allele2 && $allele2 ne '.'){ next IND }
										push @{$validalleles[$subpop]}, ($allele1, $allele2);
										}
									}
								elsif ($haplotize || $ishemi){ push @{$validalleles[$subpop]}, 0 }
								else { push @{$validalleles[$subpop]}, (0, 0) }
								}
							my $validal=(defined $validalleles[$subpop]) ? scalar(@{ $validalleles[$subpop] }) : 0;
							unless ($validal>=$minnumbers[$subpop]){ ++$nodata[$metapop]; warn "$chrom:$pos is polymorphic, but population $subpop has only $validal (min: $minnumbers[$subpop]) valid alleles!\n" if ($snp); next GROUP }
							}

						my $nder=0;
						if ($snp){
							for my $subpop ( @{$groups->[$metapop]} ){
								my $ref_suballeles=Misc::samplewithoutRep($validalleles[$subpop], $minnumbers[$subpop]);
								for my $allele (@$ref_suballeles){ if ($allele ne $ancestral){ ++$nder } }
								}
							print STDERR "$chrom:$pos, group $metapop has $nder derived alleles.\n";
							}

						if ($nder>0){
							if ($degstate eq '1'){
								++$selectedsfs[$metapop][$cat][$nder];
								++$selectedsegs[$metapop][$cat] if ($nder<$arraysizes[$metapop]);
								}
							elsif ($degstate eq '4'){
								++$neutralsfs[$metapop][$cat][$nder];
								++$neutralsegs[$metapop][$cat] if ($nder<$arraysizes[$metapop]);
								}
							}

						++$validpostot[$metapop];
						++$validpos[$metapop];
						if ($degstate eq '1'){
							if ($ancrefbase eq 'A' || $ancrefbase eq 'T'){ ++$basecounts[$metapop]{'selAT'} }
							elsif ($ancrefbase eq 'G' || $ancrefbase eq 'C'){ ++$basecounts[$metapop]{'selGC'} }
							else { ++$basecounts[$metapop]{'selN'} }
							}
						elsif ($degstate eq '4'){
							if ($ancrefbase eq 'A' || $ancrefbase eq 'T'){ ++$basecounts[$metapop]{'neutAT'} }
							elsif ($ancrefbase eq 'G' || $ancrefbase eq 'C'){ ++$basecounts[$metapop]{'neutGC'} }
							else { ++$basecounts[$metapop]{'neutN'} }
							}
						}
					}
				}

			if ($printpergene){
				print $output "$geneid\t$transcriptid";
				for my $group (0..$#{$groups}){
					print $output "\t$validpos[$group]";
					for my $cat (0..4){
#						print $output "\t$selectedsegs[$group][$cat]";
						print $output "\t", join(',', @{ $selectedsfs[$group][$cat] });
						}
					print $output "\t$basecounts[$group]{'selAT'}\t$basecounts[$group]{'selGC'}\t$basecounts[$group]{'selN'}";
					my $validbases=$basecounts[$group]{'selAT'}+$basecounts[$group]{'selGC'};
					my $gccont=$validbases>0 ? $basecounts[$group]{'selGC'}/$validbases : 0;
					print $output "\t", sprintf("%.2f", $gccont*100);
					for my $cat (0..4){
#						print $output "\t$neutralsegs[$group][$cat]";
						print $output "\t", join(',', @{ $neutralsfs[$group][$cat] });
						}
					print $output "\t$basecounts[$group]{'neutAT'}\t$basecounts[$group]{'neutGC'}\t$basecounts[$group]{'neutN'}";
					$validbases=$basecounts[$group]{'neutAT'}+$basecounts[$group]{'neutGC'};
					$gccont=$validbases>0 ? $basecounts[$group]{'neutGC'}/$validbases : 0;
					print $output "\t", sprintf("%.2f", $gccont*100);
					}
				print $output "\n";
				}

			for my $group ( 0..$#{ $groups } ){
#				for my $cat (0..4){
#					$self->{_SelectedSFS}[$group][$cat][$_]+=$selectedsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
#					$self->{_NeutralSFS}[$group][$cat][$_]+=$neutralsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
#					}
				$self->{_ValidPos}[$group][$bin]+=$validpos[$group];
				$self->{_Basecounts}[$group][$bin]{$_}+=$basecounts[$group]{$_} foreach ('selAT', 'selGC', 'selN', 'neutAT', 'neutGC', 'neutN');
				for my $cat (0..4){
					$self->{_SelectedSegs}[$group][$bin][$cat]+=$selectedsegs[$group][$cat];
					$self->{_NeutralSegs}[$group][$bin][$cat]+=$neutralsegs[$group][$cat];
					$self->{_SelectedSFS}[$group][$bin][$cat][$_]+=$selectedsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
					$self->{_NeutralSFS}[$group][$bin][$cat][$_]+=$neutralsfs[$group][$cat][$_] foreach (0..$arraysizes[$group]);
					}
				}
			}
		}

	if ($printpergene && $outfile ne "-" && $outfile ne "STDOUT"){ close $output }

	print "Total number of assessed codons: $totalcodons, invalid codon pattern: $codonpattern, validsites: ", join(',', @validpostot), ", undefined ancestral state: $undefanc, ",
		"not biallelic: $notbiallelic, no data: ", join(',', @nodata), "\n";
	return;
	}


sub countCodons{
	my ($self, $regions, $refseq, $ref_missing, $groups, $mincov, $minind)=@_;

	print "Counting codons...\n";
	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum, @groupsizes);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	my ($validcodons, $NinCodon, $filtered, $ancstop, $totalcodons)=(0, 0, 0, 0, 0);

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $frame=$regions->{_Loci}{$chrom}[$locus]{'frame'};
			my $dir=$regions->{_Loci}{$chrom}[$locus]{'dir'};
			my $locuslength=$end-$start+1;
			my $shift=$frame;
			if ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			my $ref_seqstring=$refseq->randomAccess($chrom, $start+$shift, $end, 1);

			CODON: for (my $codstart=$start+$shift; $codstart+2<=$end; $codstart+=3){
				++$totalcodons;
				my $codoffset=$codstart-$start;
				my $codon=substr($$ref_seqstring, 0, 3, "");	# chew up reference string by next three bases.
				$codon=uc($codon);
				if ($codon=~/N/){
					++$NinCodon;
					next CODON;
					}
				if ($dir eq '-'){ $codon=reverse($codon); $codon=~tr/ATGC/TACG/ }

				POS: for my $codonpos (0..2){
					my $pos=$codstart+$codonpos;
					my $offset=$pos-$start;

					GROUP: for my $metapop ( 0..$#{$groups} ){
						for my $subpop ( @{$groups->[$metapop]} ){
							my $validind=0;
							for my $ind (0..$popsizes[$subpop]-1){
								my $allind=$ind+$popsum[$subpop];
								my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
								++$validind if (defined $coverage && $coverage>=$mincov);
								}
							unless ($validind>=$minind->[$subpop]){
								++$filtered;
								next CODON;
								}
							}
						}
					}
				++$self->{_CodonCounts}{$codon};
				++$validcodons;
				}
			}
		}
	print "Total number of codons assessed: $totalcodons, validcodons: $validcodons, codons with Ns: $NinCodon, codons with insufficient coverage: $filtered\n";
	$self->{_TotalCodons}+=$totalcodons;
	$self->{_ValidCodons}+=$validcodons;
	$self->{_MissingCodons}+=$NinCodon;
	$self->{_FilteredCodons}+=$filtered;
	return;
	}


sub printSFS{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	for my $pop (0..$#{ $self->{_SFS} } ){
		print $output "Population: $pop\n";
		for my $nder (0..$#{ $self->{_SFS}[$pop] } ){
			print $output "\td${pop}_${nder}";
			}
		print $output "\n";
		for my $nder (0..$#{ $self->{_SFS}[$pop] } ){
			print $output "\t$self->{_SFS}[$pop][$nder]";
			}
		print $output "\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printJointSFSonefile{
	my ($self, $groups, $outfile)=@_;

	open my $output, ">", $outfile or die "Could not open output file $outfile!\n";

	my $poppair=0;
	POP1: for my $pop1 ( 0..($#{$groups}-1) ){
		POP2: for my $pop2 ( ($pop1+1)..$#{$groups} ){
			my @arraysize=( $#{ $self->{_SFS}[$poppair] }, $#{ $self->{_SFS}[$poppair][0] } );
			$self->{_SFS}[$poppair][0][0]+=$self->{_SFS}[$poppair][-1][-1];	# the top right bin needs to be combined with the lower left bin.
			$self->{_SFS}[$poppair][-1][-1]=0;
			for my $nder1 (0..$arraysize[0]){ print $output "\td${pop1}_${nder1}" }
			print $output "\n";
			for my $nder2 (0..$arraysize[1]){
				print $output "d${pop2}_${nder2}";
				for my $nder1 (0..$arraysize[0]){print $output "\t$self->{_SFS}[$poppair][$nder1][$nder2]"}
				print $output "\n";
				}
			} continue { ++$poppair }
		}
	close $output;
	return;
	}


sub printJointSFS{
	my ($self, $groups, $outfolder, $outfile_prefix)=@_;

	my $poppair=0;
	POP1: for my $pop1 ( 0..($#{$groups}-1) ){
		POP2: for my $pop2 ( ($pop1+1)..$#{$groups} ){
			my $outfile="$outfile_prefix" . "_jointDAFpop" . $pop2 . "_" . $pop1 . ".obs";
			open my $output, ">", "$outfolder/$outfile" or die "Could not open output file $outfile!\n";
			print $output "1 observations\n";
			my @arraysize=( $#{ $self->{_SFS}[$poppair] }, $#{ $self->{_SFS}[$poppair][0] } );
			$self->{_SFS}[$poppair][0][0]+=$self->{_SFS}[$poppair][-1][-1];	# the top right bin needs to be combined with the lower left bin.
			$self->{_SFS}[$poppair][-1][-1]=0;
			for my $nder1 (0..$arraysize[0]){ print $output "\td${pop1}_${nder1}" }
			print $output "\n";
			for my $nder2 (0..$arraysize[1]){
				print $output "d${pop2}_${nder2}";
				for my $nder1 (0..$arraysize[0]){print $output "\t$self->{_SFS}[$poppair][$nder1][$nder2]"}
				print $output "\n";
				}
			close $output;
			} continue { ++$poppair }
		}
	return;
	}


sub printSFSmultiN{
	my ($self, $outfile)=@_;
	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	for my $pop (0..$#{ $self->{_SFS} } ){
		print $output "Population: $pop\n";
		for my $set (0..$#{ $self->{_SFS}[$pop] } ){
			print $output "SFS $set, minimal sampling size ", scalar( @{ $self->{_SFS}[$pop][$set] } ), "\n";
			for my $nder (0..$#{ $self->{_SFS}[$pop][$set] } ){
				print $output "\td${pop}_${nder}";
				}
			print $output "\n";
			for my $nder (0..$#{ $self->{_SFS}[$pop][$set] } ){
				print $output "\t$self->{_SFS}[$pop][$set][$nder]";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printSFSforDFEalpha{
	my ($self, $outfile, $folded)=@_;
	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	my @nalleles;
	print $output scalar @{ $self->{_Selected} }, "\n", scalar @{ $self->{_Selected}[0] }, "\n";
	for my $set (0..$#{ $self->{_Selected}[0] } ){
		push @nalleles, $folded ? (scalar @{ $self->{_Selected}[0][$set] }-1)*2 : scalar @{ $self->{_Selected}[0][$set] }-1;
		print $output "$nalleles[-1]\n";
		}

	for my $pop (0..$#{ $self->{_Selected} } ){
		print $output $pop+1, "\n???\n???\n???\n???\n";
		for my $set (0..$#{ $self->{_Selected}[$pop] } ){
			for my $nder (0..$nalleles[$set]){
				my $count=defined $self->{_Selected}[$pop][$set][$nder] ? $self->{_Selected}[$pop][$set][$nder] : 0;	# this needs to be the selected SFS.
				print $output "$count ";
				}
			print $output "\n";
			for my $nder (0..$nalleles[$set]){
				my $count=defined $self->{_Neutral}[$pop][$set][$nder] ? $self->{_Neutral}[$pop][$set][$nder] : 0;	# this needs to be the neutral SFS.
				print $output "$count ";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printSFSforStairway{
	my ($self, $outprefix, $ref_poplabels)=@_;

	for my $pop (0..$#{ $self->{_SFS} } ){
		my $poplabel=$ref_poplabels->[$pop];
		my $outfile=$outprefix . "_" . $poplabel . ".txt";
		my $samplesize=$#{ $self->{_SFS}[$pop] };
		my $seqlength=Misc::sum(@{ $self->{_SFS}[$pop] });

		open my $output, ">", $outfile or die "Could not open output file $outfile!\n$!\n";
		print $output "$poplabel\t$samplesize\t$seqlength\t1\t", $samplesize-1, "\n";
		print $output join("\t", @{ $self->{_SFS}[$pop] }[1..$samplesize-1]);
		close $output;
		}

	return;
	}


sub printSegsperCat{
	my ($self, $outfile)=@_;
	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file!\n" }

	print $output "RecBin";
	print $output "\tValidsites\tSel_W->W\tSel_W->S\tSel_S->W\tSel_S->S\tSel_undef\tSel_AT\tSel_GC\tSel_N\tSel_\%GC\tNeut_W->W\tNeut_W->S\tNeut_S->W\tNeut_S->S\tNeut_undef\tNeut_AT\tNeut_GC\tNeut_N\tNeut_\%GC" foreach (0..$#{ $self->{_SelectedSegs} } );
	print $output "\n";
	for my $bin (0..$#{ $self->{_SelectedSegs}[0] } ){
		print $output "$bin";
		for my $group (0..$#{ $self->{_SelectedSegs} } ){
			print $output "\t$self->{_ValidPos}[$group][$bin]";
			for my $cat (0..4){
				print $output "\t$self->{_SelectedSegs}[$group][$bin][$cat]";
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'selAT'}\t$self->{_Basecounts}[$group][$bin]{'selGC'}\t$self->{_Basecounts}[$group][$bin]{'selN'}";
			my $validbases=$self->{_Basecounts}[$group][$bin]{'selAT'}+$self->{_Basecounts}[$group][$bin]{'selGC'};
			my $gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'selGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);
			for my $cat (0..4){
				print $output "\t$self->{_NeutralSegs}[$group][$bin][$cat]";
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'neutAT'}\t$self->{_Basecounts}[$group][$bin]{'neutGC'}\t$self->{_Basecounts}[$group][$bin]{'neutN'}";
			$validbases=$self->{_Basecounts}[$group][$bin]{'neutAT'}+$self->{_Basecounts}[$group][$bin]{'neutGC'};
			$gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'neutGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);
			print $output "\n";
			}
		}

	print $output "\n";
	my @cats=("W->W", "W->S", "S->W", "S->S", "undef");
	for my $group (0..$#{ $self->{_SelectedSegs} } ){
		for my $cat (0..4){
			print $output "Group $group, $cats[$cat]\n";
			print $output join("\t", @{ $self->{_SelectedSFS}[$group][$cat] }), "\n";
			print $output join("\t", @{ $self->{_NeutralSFS}[$group][$cat] }), "\n\n";
			}
		print $output "\n\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}

sub printSFSperCat{
	my ($self, $outfile)=@_;
	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file!\n" }

	print $output "RecBin";
	print $output "\tValidsites\tSel_W->W\tSel_W->S\tSel_S->W\tSel_S->S\tSel_undef\tSel_AT\tSel_GC\tSel_N\tSel_\%GC\tNeut_W->W\tNeut_W->S\tNeut_S->W\tNeut_S->S\tNeut_undef\tNeut_AT\tNeut_GC\tNeut_N\tNeut_\%GC" foreach (0..$#{ $self->{_SelectedSFS} } );
	print $output "\n";
	for my $bin (0..$#{ $self->{_SelectedSFS}[0] } ){
		print $output "$bin";
		for my $group (0..$#{ $self->{_SelectedSFS} } ){
			print $output "\t$self->{_ValidPos}[$group][$bin]";
			for my $cat (0..4){
				print $output "\t", join(',', @{ $self->{_SelectedSFS}[$group][$bin][$cat] });
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'selAT'}\t$self->{_Basecounts}[$group][$bin]{'selGC'}\t$self->{_Basecounts}[$group][$bin]{'selN'}";
			my $validbases=$self->{_Basecounts}[$group][$bin]{'selAT'}+$self->{_Basecounts}[$group][$bin]{'selGC'};
			my $gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'selGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);
			for my $cat (0..4){
				print $output "\t", join(',', @{ $self->{_NeutralSFS}[$group][$bin][$cat] });
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'neutAT'}\t$self->{_Basecounts}[$group][$bin]{'neutGC'}\t$self->{_Basecounts}[$group][$bin]{'neutN'}";
			$validbases=$self->{_Basecounts}[$group][$bin]{'neutAT'}+$self->{_Basecounts}[$group][$bin]{'neutGC'};
			$gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'neutGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);
			}
		print $output "\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}

sub readSFSperCat{
	my ($self, $infile)=@_;
	my $input;
	if ($infile eq "-" || $infile eq "STDIN"){ $input=*STDIN } else { open $input, "<", $infile or die "Could not open input file!\n" }

	my $ngroups;

	LINE: while (<$input>){
		next LINE if (/^RecBin/);
		if (/^\d+\s\d+/){
			chomp;
			my @line=split('\t', $_);
			my $bin=shift(@line);
			my $ngroups=(@line-1)/15 unless defined $ngroups;
			die "Rows in input file differ in the number of fields!\n" if ($ngroups!=(@line-1)/15);
			for my $group (0..$ngroups){
				$self->{_ValidPos}[$group][$bin]=shift(@line);
				@{ $self->{_SelectedSFS}[$group][$bin][$_] }=split(',', shift(@line) ) foreach (0..3);
				$self->{_Basecounts}[$group][$bin]{'selAT'}=shift(@line);
				$self->{_Basecounts}[$group][$bin]{'selGC'}=shift(@line);
				shift(@line);
				@{ $self->{_NeutralSFS}[$group][$bin][$_] }=split(',', shift(@line) ) foreach (0..3);
				$self->{_Basecounts}[$group][$bin]{'neutAT'}=shift(@line);
				$self->{_Basecounts}[$group][$bin]{'neutGC'}=shift(@line);
				print STDERR "Reading recombination bin $bin ...\n";
				}
			}
		}

	if ($infile ne "-" && $infile ne "STDIN"){ close $input }
	return;
	}

sub perfomThetaWBootstrap{
	my ($self, $nbootstraps, $nchrom, $outfile)=@_;

	my $denom=0;
	$denom+=(1/$_) foreach (1..$nchrom-1);

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file!\n" }

	print $output "RecBin";
	print $output "\tValidsites\tSel_W->W\tSel_W->S\tSel_S->W\tSel_S->S\tSel_All\tSel_AT\tSel_GC\tSel_\%GC\tNeut_W->W\tNeut_W->S\tNeut_S->W\tNeut_S->S\tNeut_All\tNeut_AT\tNeut_GC\tNeut_\%GC" foreach (0..$#{ $self->{_SelectedSFS} } );
	print $output "\n";
	for my $bin (0..$#{ $self->{_SelectedSFS}[0] } ){
		print STDERR "Working on recombination bin $bin ...\n";
		print $output "$bin";
		for my $group (0..$#{ $self->{_SelectedSFS} } ){
			print $output "\t$self->{_ValidPos}[$group][$bin]";
			my $validbases=$self->{_Basecounts}[$group][$bin]{'selAT'}+$self->{_Basecounts}[$group][$bin]{'selGC'};
			my $nsnps=0;
			my @sampling=(0) x 4;;
			for my $cat (0..3){
				my @sfs=@{ $self->{_SelectedSFS}[$group][$bin][$cat] };
				$self->{_SelectedBootstrap}[$group][$bin][$cat]=[ (0) x $nbootstraps ];
				$nsnps+=Misc::sum(@sfs[1..$#sfs-1]);
				$sampling[$cat]=$nsnps;
				}
			print STDERR "validbases: $validbases, nsnps: $nsnps\n";
			print STDERR join(',', @sampling), "\n";
			for my $rep (0..$nbootstraps-1){
				my @counts=(0) x 5;
				for my $base (0..$validbases-1){
					my $rand=int( rand($validbases) );
					if ($rand>=$nsnps){ ++$counts[4] }
					else {
						my $cat=0;
						++$cat while ($rand>=$sampling[$cat]);
						if ($cat>3){
							print STDERR "validbases: $validbases, nsnps: $nsnps\n";
							print STDERR join(',', @sampling), "\n";
							die "Wrong category assignment $cat!\n";
							}
						++$counts[$cat];
						}
					}
				$self->{_SelectedBootstrap}[$group][$bin][$_][$rep]=$counts[$_]/($self->{_Basecounts}[$group][$bin]{'selAT'}*$denom) foreach (0..1);
				$self->{_SelectedBootstrap}[$group][$bin][$_][$rep]=$counts[$_]/($self->{_Basecounts}[$group][$bin]{'selGC'}*$denom) foreach (2..3);
				$self->{_SelectedBootstrap}[$group][$bin][4][$rep]=Misc::sum(@counts[0..3])/($validbases*$denom);
				}

			for my $cat (0..4){
				print $output "\t", join(',', @{ $self->{_SelectedBootstrap}[$group][$bin][$cat] });
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'selAT'}\t$self->{_Basecounts}[$group][$bin]{'selGC'}";
			my $gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'selGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);

			$validbases=$self->{_Basecounts}[$group][$bin]{'neutAT'}+$self->{_Basecounts}[$group][$bin]{'neutGC'};
			$nsnps=0;
			@sampling=(0) x 4;
			for my $cat (0..3){
				my @sfs=@{ $self->{_NeutralSFS}[$group][$bin][$cat] };
				$self->{_NeutralBootstrap}[$group][$bin][$cat]=[ (0) x $nbootstraps ];
				$nsnps+=Misc::sum(@sfs[1..$#sfs-1]);
				$sampling[$cat]=$nsnps;
				}
			print STDERR "validbases: $validbases, nsnps: $nsnps\n";
			print STDERR join(',', @sampling), "\n";
			for my $rep (0..$nbootstraps-1){
				my @counts=(0) x 5;
				for my $base (0..$validbases-1){
					my $rand=int( rand($validbases) );
					if ($rand>=$nsnps){ ++$counts[4] }
					else {
						my $cat=0;
						++$cat while ($rand>=$sampling[$cat]);
						if ($cat>3){
							print STDERR "validbases: $validbases, nsnps: $nsnps\n";
							print STDERR join(',', @sampling), "\n";
							die "Wrong category assignment $cat!\n";
							}
						++$counts[$cat];
						}
					}
				$self->{_NeutralBootstrap}[$group][$bin][$_][$rep]=$counts[$_]/($self->{_Basecounts}[$group][$bin]{'neutAT'}*$denom) foreach (0..1);
				$self->{_NeutralBootstrap}[$group][$bin][$_][$rep]=$counts[$_]/($self->{_Basecounts}[$group][$bin]{'neutGC'}*$denom) foreach (2..3);
				$self->{_NeutralBootstrap}[$group][$bin][4][$rep]=Misc::sum(@counts[0..3])/($validbases*$denom);
				}

			for my $cat (0..4){
				print $output "\t", join(',', @{ $self->{_NeutralBootstrap}[$group][$bin][$cat] });
				}
			print $output "\t$self->{_Basecounts}[$group][$bin]{'neutAT'}\t$self->{_Basecounts}[$group][$bin]{'neutGC'}";
			$gccont=$validbases>0 ? $self->{_Basecounts}[$group][$bin]{'neutGC'}/$validbases : 0;
			print $output "\t", sprintf("%.2f", $gccont*100);
			}
		print $output "\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub printCodonCounts{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	print $output "Total number of codons assessed: $self->{_TotalCodons}, valid codons: $self->{_ValidCodons}, codons with Ns: $self->{_MissingCodons}, codons with insufficient coverage: $self->{_FilteredCodons}\n";

	while ( my ($codon, $count)=each %{$self->{_CodonCounts} } ){
		print $output "$codon: $count\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printFixedShared{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file!\n" }

	print $output "$self->{_nsnps_ancdef}"; 
	for my $set (@{ $self->{_Fixed} }){
		print $output "\t$set"; 
		}
	print $output "\n";

	print $output "$self->{_nsnps_ancdef}"; 
	for my $set (@{ $self->{_Derived} }){
		print $output "\t$set"; 
		}
	print $output "\n";

	print $output "$self->{_nsnps_all}"; 
	for my $set (@{ $self->{_Shared} }){
		print $output "\t$set"; 
		}
	print $output "\n";

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output } 
	return;
	}


sub summarizeStats{
	my ($self, $folder, $outfile, $ref_header, $ref_sets, $prefix, $glob, @idlist)=@_;

	my @filelist;
	if ($glob){ @filelist=<$folder/${prefix}_*.txt> }
	else {
		for my $id (@idlist){
			my $filename="$folder" . "/" . "$prefix" . '_' . "$id" . '.txt';
			push @filelist, $filename;
			}
		}
	print join(', ', @filelist), "\n";

	my @sets;
	for my $filename (@filelist){
		print "Summarizing stats file $filename.\n";
		open my $fh, "<", $filename or die "Could not open $filename!\n";
		my $line=0;
		LINE: while (<$fh>){
			if (/^\d+/){
				chomp;
				my @set=split(/\s/, $_);
				$sets[$line][$_]+=$set[$_] foreach (0..$#set);
				++$line;
				}
			}
		close $fh;
		}
	unless (@sets==@$ref_sets){ die "Wrong numer of columns in stats files!\n" }
	unless (@{ $sets[0] }==@$ref_header){ die "Wrong numer of rows in stats files!\n" }

	my @areanames=("area1", "area2", "area3", "area4", "n12", "n13", "n14", "n23", "n24", "n34", "n123", "n124", "n134", "n234", "n1234");
	my @areas;
	$areas[$_]=[ (0) x scalar(@areanames) ] foreach (0..$#sets);
	my @inter=(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15);
	for my $line (0..$#sets){
		for my $bitset ( 1..scalar(@{ $sets[$line] })-2 ){
			for my $i (0..$#inter){
#				if (($bitset & $inter[$i])==$inter[$i]){ print "Line: $line, bitset: $bitset, inter: $inter[$i], ", $bitset & $inter[$i], " - Adding $sets[$line][$bitset+1]\n" }
				$areas[$line][$i]+=$sets[$line][$bitset+1] if ( ($bitset & $inter[$i])==$inter[$i] );
				}
			}
		}

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", "$folder/$outfile" or die "Could not open output file!\n" }
	print $output "\t", join( "\t", (@$ref_header, @areanames) ), "\n";
	for my $line (0..$#sets){
		print $output $ref_sets->[$line];
		for my $set ( @{ $sets[$line] } ){
			print $output "\t$set"; 
			}
		for my $area ( @{ $areas[$line] } ){
			print $output "\t$area"; 
			}
		print $output "\n";
		}

	for my $line (0..$#sets){
		print $output $ref_sets->[$line];
		my @counts=(@{ $sets[$line] }, @{ $areas[$line] });
		my $nsnps=shift @counts;
		print $output "\t$nsnps";
		for my $set (@counts){
			print $output "\t", sprintf("%.2f", $set/$nsnps*100); 
			}
		print $output "\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output } 
	return;
	}


sub summarizeSFS{
	my ($self, $groups, $minind, $folded, $haplotize, $folder, $prefix, $glob, @idlist)=@_;

	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;
	my @arraysizes;
	for my $group ( 0..$#{$groups} ){
		if ($haplotize && $folded){ $arraysizes[$group]=int($groupsizes[$group]/2) }
		elsif ($haplotize || $folded){ $arraysizes[$group]=$groupsizes[$group] }
		else { $arraysizes[$group]=2*$groupsizes[$group] }
		}

	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		for my $pop ( 0..($#{$groups}) ){
			$self->{_SFS}[$pop]=[ (0) x ($arraysizes[$pop]+1) ];
			}
		}

	my @filelist;
	if ($glob){ @filelist=<$folder/${prefix}_*.txt> }
	else {
		for my $id (@idlist){
			my $filename="$folder" . "/" . "$prefix" . '_' . "$id" . '.txt';
			push @filelist, $filename;
			}
		}
	print join(', ', @filelist), "\n";

	for my $filename (@filelist){
		print "Summarizing SFS file $filename.\n";
		open my $fh, "<", $filename or die "Could not open $filename!\n";
		my $pop;
		LINE: while (<$fh>){
			if (/Population: (\d+)/){ $pop=$1 }
			elsif (/^\td([0-9]+)_[0-9]+/){
				unless ($pop==$1){ die "Wrong format of SFS file $filename, line $.!\n" }
				my @line=split("\t", $_);
				shift(@line);
				unless (@line==$arraysizes[$pop]+1){ die "Wrong number of columns in SFS file $filename, ", scalar(@line), " vs. ", $arraysizes[$pop]+1, "!\n" }
				}
			elsif (/^\t\d+/){
				my @line=split("\t", $_);
				shift(@line);
				unless (@line==$arraysizes[$pop]+1){ die "Wrong number of columns in SFS file $filename, ", scalar(@line), " vs. ", $arraysizes[$pop]+1, "!\n" }
				for my $col (0..$#line){ $self->{_SFS}[$pop][$col]+=$line[$col] }
				}
			else { die "Wrong format of SFS file $filename, line $.!\n" }
			}
		close $fh;
		}

	return;
	}

sub summarizejointSFS{
	my ($self, $groups, $minind, $folded, $haplotize, $folder, $prefix, $glob, @idlist)=@_;

	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;
	my @arraysizes;
	for my $group ( 0..$#{$groups} ){
		if ($haplotize && $folded){ $arraysizes[$group]=int($groupsizes[$group]/2) }
		elsif ($haplotize || $folded){ $arraysizes[$group]=$groupsizes[$group] }
		else { $arraysizes[$group]=2*$groupsizes[$group] }
		}

	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		my $poppair=0;
		for my $pop1 ( 0..($#{$groups}-1) ){
			for my $pop2 ( ($pop1+1)..$#{$groups} ){
				for my $i (0..$arraysizes[$pop1]){ $self->{_SFS}[$poppair][$i]=[ (0) x ($arraysizes[$pop2]+1) ] }
				++$poppair;
				}
			}
		}

	my @filelist;
	if ($glob){ @filelist=<$folder/${prefix}_*.txt> }
	else {
		for my $id (@idlist){
			my $filename="$folder" . "/" . "$prefix" . '_' . "$id" . '.txt';
			push @filelist, $filename;
			}
		}
	print join(', ', @filelist), "\n";

	for my $filename (@filelist){
		print "Summarizing SFS file $filename.\n";
		open my $fh, "<", $filename or die "Could not open $filename!\n";
		my ($poppair, $ncols, $nrows, $pop1, $pop2);
		LINE: while (<$fh>){
			if (/^\td([0-9]+)_[0-9]+\td[0-9]+_[0-9]+/){
				if (defined $pop2 && $nrows != $arraysizes[$pop2]+1){ die "Wrong number of rows in SFS file for group $pop2, $nrows vs. ", $arraysizes[$pop2]+1, "!\n" }
				if (defined $poppair){ ++$poppair }
				else { $poppair=0 }
				my @line=split("\t", $_);
				$pop1=$1;
				$ncols=@line;
				unless ($ncols==$arraysizes[$pop1]+2){ die "Wrong number of columns in SFS file $filename, $ncols vs. ", $arraysizes[$pop1], "!\n" }
				$nrows=0;				
				}
			elsif (/^d([0-9]+)_([0-9]+)\t\d+/){
				unless (defined $poppair){ die "Wrong format of SFS file $filename, line $.!\n" }
				$pop2=$1;
				++$nrows;
				my $nder2=$2;
				unless ($nrows-1==$nder2){ die "Wrong number of rows in SFS file $filename, line $.!\n" }
				my @line=split("\t", $_);
				unless (@line==$ncols){ die "Wrong number of columns in SFS file $filename, line $.!\n" }
				for my $col (1..$#line){ $self->{_SFS}[$poppair][$col-1][$nder2]+=$line[$col] }
				}
			else { die "Wrong format of SFS file $filename, line $.!\n" }
			}
		close $fh;
		}

	return;
	}



1;

