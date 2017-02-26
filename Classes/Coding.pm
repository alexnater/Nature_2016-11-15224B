package Coding;

use strict;
use warnings;
use Classes::Regions;
use Classes::Misc;


sub new{
	my $class=shift;
	my $self={};
	bless $self, $class;
	return $self;
	}


sub readCDS{
	my ($self, $filename, $keepinfo, $ref_bins, $ref_ids)=@_;
	unless(@_>=3){ die "Wrong number of arguments!\n" }
	my %genes;	# %genes{'geneid'}{'transcriptid'}[exon]{'scaffold', 'start', 'end', 'dir', 'frame'}
	my %bins=map { $_=>1 } @$ref_bins if defined $ref_bins;
	my %ids=map { $_=>1 } @$ref_ids if defined $ref_ids;

	open my $input, "<", $filename or die "Could not open bed file $filename!\n";

	LINE: while (my $line=<$input>){
		chomp $line;
		if ($line=~/^\w+\t\w+\t\w+/){
			my @col=split(/\t/, $line);
			if ($col[2] eq 'CDS'){
				my ($gene_id)=$col[8]=~/gene_id\s+"(\w+)";/;
				my ($transcript_id)=$col[8]=~/transcript_id\s+"(\w+)";/;
				next LINE unless (defined $gene_id && defined $transcript_id);
				if (defined $ref_bins && $ref_bins->[0] ne 'all'){
					my $bin=$self->{_Bins}{$gene_id};
					next LINE unless ( defined $bin && $bins{$bin} );
					}
				if (defined $ref_ids && $ref_ids->[0] ne 'all'){
					next LINE unless ($ids{$col[0]});
					}
				next LINE if ( exists $self->{_Transcripts} && !$self->{_Transcripts}{$transcript_id} );
				push( @{ $genes{$gene_id}{$transcript_id} }, { 'id'=>$col[0], 'start'=>($col[3]-1), 'end'=>($col[4]-1), 'dir'=>$col[6], 'frame'=>$col[7] } );
				if ($keepinfo==1){ $genes{$gene_id}{$transcript_id}[-1]{'info'}=$col[8] }
				}
			}
		}
	close $input;
	$self->{_Genes}=\%genes;
	print STDERR "GTF file $filename successfully read in.\n";
	return;
	}


sub readBins{
	my ($self, @filenames)=@_;
	unless(@_>=2){die "Wrong number of arguments!\n"}

	for my $filename (@filenames){
		open my $input, "<", $filename or die "Could not open bin file $filename!\n";

		while (my $line=<$input>){
			chomp $line;
			if ($line=~/^"?(\w+)"?\s+"?(\w+)"?\s+(\d+)/){
				$self->{_Transcripts}{$2}=1;
				if (exists $self->{_Bins}{$1}){
					if ($self->{_Bins}{$1} ne $3){ warn "GeneID $1 has already been assigned to a recombination bin! ", $self->{_Bins}{$1}, " vs. $3\n" }
					}
				else { $self->{_Bins}{$1}=$3 }
				}
			}
		close $input;
		}
	return;
	}

sub getAllBins{
	my $self=shift;
	my %bins;
	++$bins{$_} foreach ( values %{ $self->{_Bins} } );
#	my %bins=map { $_=>1 } ( values %{ $self->{_Bins} } );
	my @allbins=sort {$a<=>$b} keys %bins;
	print "RecBin $_: $bins{$_}\n" foreach (@allbins);

	return(@allbins);
	}


sub addRefseq{
	my ($self, $fasta)=@_;
	GENE:	for my $ref_gene ( values %{ $self->{_Genes} } ){
		TRANS: for my $ref_transcript ( values %{ $ref_gene } ){
			EXON: for my $ref_exon (@$ref_transcript){
				my $ref_seqstring=$fasta->randomAccess($ref_exon->{'id'}, $ref_exon->{'start'}, $ref_exon->{'end'}, 1);
				$ref_exon->{'refseq'}=$$ref_seqstring;
				print "Added reference sequence to $ref_exon->{'id'}:$ref_exon->{'start'}-$ref_exon->{'end'}.\n";
				}
			}
		}
	return;
	}


sub concatenateTranscripts{
	my ($self, $fasta, $ancref_fasta, $keeplongest, $conflict_response)=@_;
	unless(@_>=2){ die "Wrong number of arguments!\n" }
	my %transcripts;	# %transcripts{scaffold}[locus]{'geneid', 'transcriptid', 'start', 'end', 'dir', 'refseq', 'positions'}

	GENE: for my $gene ( keys %{ $self->{_Genes} } ){
		my %temptranscripts;
		TRANS: for my $transcript ( keys %{ $self->{_Genes}{$gene} } ){
			my ($first_exon, $prev_exon, $clip, $refseq, $ancref);
			my $prev_clip=0;
			my @positions;
			EXON: for my $exon ( sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Genes}{$gene}{$transcript} } ){
				$first_exon=$exon unless (defined $first_exon);
				my $id=$exon->{'id'};
				my $start=$exon->{'start'};
				my $end=$exon->{'end'};
				if ( $id ne $first_exon->{'id'} ){
					warn "Gene $gene, transcript $transcript is not on the same scaffold for exon $first_exon->{'id'}:$first_exon->{'start'}-$first_exon->{'end'} and $id:$start-$end\n";
					next TRANS;
					}
				my $dir=$exon->{'dir'};
				if ( $dir ne $first_exon->{'dir'} ){
					warn "Gene $gene, transcript $transcript has not the same directions for exon $first_exon->{'id'}:$first_exon->{'start'}-$first_exon->{'end'} and $id:$start-$end\n";
					next TRANS;
					}
				my $frame=$exon->{'frame'};
				my $locuslength=$end-$start+1;
				my $leftover=($locuslength-$frame) % 3;
				my $shift=($dir eq '-') ? $leftover : $frame;
				$clip=($dir eq '-') ? $frame : $leftover;
				if ( defined $prev_exon && ($prev_clip+$shift) % 3 ){
					warn "Gene $gene, transcript $transcript has non-matching codon bounderies for exon ",
						"$prev_exon->{'id'}:$prev_exon->{'start'}-$prev_exon->{'end'} and $id:$start-$end! $dir, $frame, $prev_clip, $shift\n";
					unless ($locuslength-$frame){
						$exon=$prev_exon; $clip=$prev_clip;
						warn "Skipping codon, since it has a residual length of zero.\n";
						next EXON;
						}
					unless ($conflict_response){ die "Critical error since residual length is not zero!\n" }
					if ($conflict_response==1){
						warn "Residual length is not zero, skipping transcript!\n";
						next TRANS;
						}
					else {
						warn "Residual length is not zero, merging exons by removing overlaps!\n";
						$start+=$shift;
						$locuslength-=$shift;
						if ($prev_clip>0){
							$refseq=substr($refseq, 0, -$prev_clip);
							$ancref=substr($ancref, 0, -$prev_clip) if length($ancref);
							splice(@positions, -$prev_clip);
							}
						unless (@positions==length($refseq) ){
							die "Error: Unequal length of refeq string and positions array detected! ", length($refseq), " vs. ", scalar(@positions), "\n";
							}
						}
					}
				my $ref_exonseq=$fasta->randomAccess($id, $start, $end, 1);
				unless (length($$ref_exonseq)==$locuslength){ die "Wrong length of extracted reference sequence string!\n" }
				$refseq.=$$ref_exonseq;
				my $ref_exonancseq;
				if (defined $ancref_fasta){
					$ref_exonancseq=$ancref_fasta->randomAccess($id, $start, $end, 1);
					unless (length($$ref_exonancseq)==$locuslength){ die "Wrong length of extracted ancestral reference sequence string!\n" }
					$ancref.=$$ref_exonancseq;
					}
				push @positions, ($start..$end);
				} continue { $prev_exon=$exon; $prev_clip=$clip }

			my $transframe=($first_exon->{'dir'} eq '-') ? $prev_exon->{'frame'} : $first_exon->{'frame'};
			my $transleftover=(length($refseq)-$transframe) % 3;
			if (length($refseq) % 3){ warn "Gene $gene, transcript $transcript is not a multiple of 3!\n" }
			if ($transframe){ warn "Gene $gene, transcript $transcript does not start with a full codon!\n" }
			if ($transleftover){ warn "Gene $gene, transcript $transcript has $transleftover leftover bases at the end!\n" }
			if ($prev_clip){ warn "Gene $gene, transcript $transcript has $prev_clip leftover bases at the end, looking in forward direction!\n" }
			unless ( @positions==length($refseq) ){
				die "Gene $gene, transcript $transcript has not the same length for reference sequence and position array! ", length($refseq), " vs. ", scalar(@positions), "\n";
				}
			my $bin=defined $self->{_Bins}{$gene} ? $self->{_Bins}{$gene} : undef;
			unless (defined $bin){ warn "Gene $gene does not have a recombination rate bin defined!\n" }

			push @{ $temptranscripts{ $first_exon->{'id'} } }, { 'start'=>$first_exon->{'start'}, 'end'=>$prev_exon->{'end'},
											 'geneid'=>$gene, 'transcriptid'=>$transcript, 'dir'=>$first_exon->{'dir'}, 'frame'=>$transframe,
											 'refseq'=>$refseq, 'ancref'=>$ancref, 'positions'=>\@positions, 'bin'=>$bin };
			print "Gene $gene, transcript $transcript concatenated.\n",
				"Location: $first_exon->{'id'}:$first_exon->{'start'}-$prev_exon->{'end'}, number of codons: ", length($refseq)/3, " direction: $first_exon->{'dir'}, frame: $transframe\n";
			}

		my ($id_longest, $ref_longest);
		my $maxlength=0;
		for my $id (keys %temptranscripts){
			for my $ref_locus (@{ $temptranscripts{$id} }){
				if (length($ref_locus->{'refseq'})>$maxlength){
					$id_longest=$id;
					$ref_longest=$ref_locus;
					$maxlength=length($ref_locus->{'refseq'});
					}
				push @{ $transcripts{$id} }, $ref_locus unless ($keeplongest);
				}
			}
		if ($keeplongest && $maxlength){
			push @{ $transcripts{$id_longest} }, $ref_longest;
			print "Gene $gene: Keeping transcript $ref_longest->{'transcriptid'}, $id_longest:$ref_longest->{'start'}-$ref_longest->{'end'}, length: $maxlength.\n";
			}
		}

	my $transcripts_loci=Regions->new(\%transcripts);
	return($transcripts_loci);
	}


sub checkStopCodon{
	my ($self, $regions, $ref_genotypes)=@_;

	my $ref_deg=Misc::getDegRef();
	my ($codonpattern, $refstop, $fixstop, $polystop, $totalcodons)=(0, 0, 0, 0, 0);
	my %stop_codons;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ref_positions=$regions->getValue($chrom, $locus, 'positions');
			my @positions;

			if ($dir eq '-'){
				@positions=reverse(@$ref_positions);	# copy and reverse the array of positions.
				$refseq=uc( reverse($refseq) );	# copy and reverse the reference string.
				}
			else {
				@positions=@$ref_positions;
				$refseq=uc($refseq);
				}

			my $locuslength=length($refseq);
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			substr($refseq, 0, $shift, "") if $shift;
			substr($refseq, -$clip, $clip, "") if $clip;
			if (length($refseq) % 3){ warn "Reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			splice(@positions, 0, $shift) if $shift;
			splice(@positions, -$clip) if $clip;
			if (@positions % 3){ warn "Position array has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }

			CODON: while ( length($refseq)>=3 ){
				++$totalcodons;
				my $codon_plusstrand=substr($refseq, 0, 3, "");	# chew up reference string by next three bases.
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
				my $refcodon=$codon_plusstrand;
				$refcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $refdeg=$ref_deg->{$refcodon};
				if (!defined $refdeg){ warn "Warning, degeneracy pattern for reference codon $refcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $basedeg=$ref_deg->{$basecodon};
				if (!defined $basedeg){ warn "Warning, degeneracy pattern for base codon $basecodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				if (!defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				if ($refdeg eq 'stop' || $basedeg eq 'stop' || $altdeg eq 'stop'){
					if (@positions){	# still codons left in transcript.
						my $codonsleft=scalar(@positions)/3;	# number of full codons left in the transcript.
						print "Warning, premature stop codon ($refdeg/$basedeg/$altdeg) detected at $chrom:$codstart-$codend! CDS end at $positions[-1], ",
							"$codonsleft codons left in transcript.\n";
						if ($refdeg eq 'stop'){
							++$refstop;
							print "Stop codon $refcodon in reference sequence at position $chrom:$codstart-$codend!\n";
							}
						if ($basedeg eq 'stop'){
							++$fixstop;
							print "Fixed stop codon $basecodon at position $chrom:$codstart-$codend!\n";
							}
						if ($altdeg eq 'stop'){
							++$polystop;
							print "Polymorphic stop codon $altcodon at position $chrom:$codstart-$codend!\n";
							}
						push @{ $stop_codons{$chrom} }, { 'start'=>$codstart, 'end'=>$codend, 'geneid'=>$geneid, 'transcriptid'=>$transcriptid,
												'dir'=>$dir, 'refseq'=>$codon_plusstrand, 'positions'=>\@triplet };
						}
#					next LOCUS;
					}
				}
			}
		}
	print "Total number of assessed codons: $totalcodons, invalid codon pattern: $codonpattern, ",
		"premature referenece stop codons: $refstop, premature fixed stop codons: $fixstop, premature polymorphic stop codons: $polystop\n";
	my $loci=Regions->new(\%stop_codons);
	return($loci);
	}


sub annotateSNPs{
	my ($self, $regions, $ref_genotypes, $reportfile)=@_;

	my $output;
	if (defined $reportfile){ open $output, ">", $reportfile or die "Could not open output file $reportfile!$!\n" }

	print "Annotating polymorphic sites ...\n";
	my $ref_codon_hash=Misc::getCodonRef();
	my ($totalcodons, $codonpattern, $undefanc, $notbiallelic, $conflicting)=(0, 0, 0, 0, 0);

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ref_positions=$regions->getValue($chrom, $locus, 'positions');
			my @positions;

			if ($dir eq '-'){
				@positions=reverse(@$ref_positions);	# copy and reverse the array of positions.
				$refseq=uc( reverse($refseq) );	# copy and reverse the reference string.
				}
			else {
				@positions=@$ref_positions;
				$refseq=uc($refseq);
				}

			my $locuslength=length($refseq);
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			if ($shift){
				substr($refseq, 0, $shift, "");
				splice(@positions, 0, $shift);
				}
			if ($clip){
				substr($refseq, -$clip, $clip, "");
				splice(@positions, -$clip);
				}
			if (length($refseq) % 3){ warn "Reference sequence has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			if (@positions % 3){ warn "Position array has wrong length for gene $geneid, transcript $transcriptid, $chrom:$start-$end!\n" }
			print "$chrom:$start-$end, ", length($refseq), ", ", scalar(@positions), ", $shift, $clip\n";
			unless (scalar(@positions)==length($refseq) ){ die "Lengths of reference sequence and position array do not match!\n" }

			my ($change, $poly, $basecodon, $altcodon, $baseaa, $altaa, @triplet);
			CODON: while ( length($refseq)>=3 ){
				++$totalcodons;
				@triplet=();
				$change='undefined';
				$poly=0;

				my $codon_plusstrand=substr($refseq, 0, 3, "");	# chew up reference string by next three bases.
				@triplet=splice(@positions, 0, 3);	# chew up position array by next three bases.
				unless (@triplet==3){ die "Invalid extraction of positions from array!\n" }
				my $codstart=($dir eq '-') ? $triplet[-1] : $triplet[0];
				my $codend=($dir eq '-') ? $triplet[0] : $triplet[-1];
				if ($codon_plusstrand=~/N/){	# if one of the codon positions is hardmasked.
					++$codonpattern;
					next CODON;
					}
				$basecodon=$codon_plusstrand;
				$altcodon=$codon_plusstrand;
				POS: for my $codonpos (0..2){
					my $pos=$triplet[$codonpos];
					if ( defined $ref_genotypes->getValue($chrom, $pos) ){
						if ($poly){ warn "Warning: Codon $chrom:$codstart-$codend has more than one polymorphic position!\n"; ++$codonpattern; next CODON }
						if ($ref_genotypes->getBase($chrom, $pos)>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$codonpattern; next CODON }
						substr($altcodon, $codonpos, 1)=$ref_genotypes->getBase($chrom, $pos, 1);
						$poly=1;
						}
					}
				next CODON unless $poly;	# no polymorphic positions in codon.

				$basecodon=~tr/ATGC/TACG/ if ($dir eq '-');
				$baseaa=$ref_codon_hash->{$basecodon};
				if (!defined $baseaa){ warn "Warning, amino acid for codon $basecodon at position at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				$altaa=$ref_codon_hash->{$altcodon};
				if (!defined $altaa){ warn "Warning, amino acid for alternative codon $altcodon at $chrom:$codstart-$codend is not defined!\n"; ++$codonpattern; next CODON }
				$change=($baseaa eq $altaa) ? 'synonymous' : 'nonsynonymous';
				if ($baseaa eq '_' || $altaa eq '_'){
					$change='polystop' if ($baseaa ne $altaa);
					if (@positions){ warn "Warning, premature stop codon ($baseaa/$altaa) detected at $chrom:$codstart-$codend! CDS end at $end.\n" }
					else { next CODON }	# regular stop codon, proceed to recording.
					}
				}
			continue {
				if ($poly){
					POS: for my $codonpos (0..2){
						my $pos=$triplet[$codonpos];
						next POS unless (defined $ref_genotypes->getValue($chrom, $pos) );
						if (!exists $self->{_Annotations}{$chrom}{$pos}){ $self->{_Annotations}{$chrom}{$pos}=$change }
						else {
							unless ($change eq $self->{_Annotations}{$chrom}{$pos} || $change eq 'undefined'){
								$self->{_Annotations}{$chrom}{$pos}='conflicting';
								++$conflicting;
								}
							}
						print $output "$chrom\t", $pos+1, "\t$change\t${basecodon}/${altcodon}\t${baseaa}/${altaa}\t$geneid\t$transcriptid\t$start\t$end\t$dir\t$frame\n";
						}
					}
				}
			}
		}

	close $output if (defined $reportfile);
	print "Total number of assessed codons: $totalcodons, invalid codon pattern: $codonpattern, undefined ancestral state: $undefanc, not biallelic: $notbiallelic, conflicting annotations: $conflicting\n";
	return;
	}


sub AssessStopCodons{
	my ($self, $regions, $ref_missing, $ref_genotypes, $groups, $mincov, $minind)=@_;
	my $ref_deg=Misc::getDegRef();
	my %stop_positions;

	my @groupsizes=map { Misc::sum(@{$minind}[ @$_ ]) } @$groups;
	unless ( exists $self->{_SFS} ){	# to prevent overwriting of previously stored counts.
		$self->{_SFS}[$_]=[ (0) x (2*$groupsizes[$_]+1) ] foreach (0..$#groupsizes);
		}
	unless ( exists $self->{_jointSFS} ){	# to prevent overwriting of previously stored counts.
		my $poppair=0;
		for my $group1 ( 0..($#{$groups}-1) ){
			for my $group2 ( ($group1+1)..$#{$groups} ){
				$self->{_jointSFS}[$poppair][$_]=[ (0) x (2*$groupsizes[$group2]+1) ] foreach (0..2*$groupsizes[$group1]);
				++$poppair;
				}
			}
		}


	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $geneid=$regions->getValue($chrom, $locus, 'geneid');
			my $transcriptid=$regions->getValue($chrom, $locus, 'transcriptid');
			my $refseq=$regions->getValue($chrom, $locus, 'refseq');
			my $ref_triplet=$regions->getValue($chrom, $locus, 'positions');
			my @triplet=@$ref_triplet;
			my $codstart=($dir eq '-') ? $triplet[-1] : $triplet[0];
			my $codend=($dir eq '-') ? $triplet[0] : $triplet[-1];

			print "Assessing gene id $geneid, transcript $transcriptid, codon $chrom:$codstart-$codend for the occurrence of stop codons...\n";
			my $stopcodon=0;
			my $npops=$ref_missing->getNPop();
			my $allind=0;
			my @validgroup=( (1) x scalar(@$groups) );
			my @nder=( (0) x scalar(@$groups) );
			my @validspecal=( (0) x scalar(@$groups) );
			my @speciesfreq=( ("NA") x scalar(@$groups) );
			GROUP: for my $metapop ( 0..$#{$groups} ){
				my $speciesstop=0;
				SUBPOP: for my $subpop ( @{$groups->[$metapop]} ){
					my @validalleles;
					my $stopcodons=0;
					IND: for my $ind ( 0..$ref_missing->getNInd($subpop)-1 ){
						my $indcodon1=$refseq;
						my $indcodon2=$refseq;
						my $hetero=0;
						POS: for my $codonpos (0..2){
							my $pos=$triplet[$codonpos];
							my $coverage=$ref_missing->getValuebyPos($chrom, $pos, $allind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }	# skip individual if any codon position is not sufficiently covered.
							if ( defined $ref_genotypes->getValue($chrom, $pos) ){
								my $index=2*$ind;
								my $base1=$ref_genotypes->getBase($chrom, $pos, $subpop, $index);
								unless (defined $base1 && $base1 ne 'N'){ next IND }
								my $base2=$ref_genotypes->getBase($chrom, $pos, $subpop, ++$index);
								unless (defined $base2 && $base2 ne 'N'){ next IND }
								if ($base1 ne $base2){
									next IND if $hetero;	# Individual has already a heterozygous position in codon.
									$hetero=1;
									}
								substr($indcodon1, $codonpos, 1)=$base1;
								substr($indcodon2, $codonpos, 1)=$base2;										
								}
							}
						$indcodon1=~tr/ATGC/TACG/ if ($dir eq '-');
						if ($ref_deg->{$indcodon1} eq 'stop'){ ++$stopcodons; push(@validalleles, 1) }
						else { push(@validalleles, 0) }
						$indcodon2=~tr/ATGC/TACG/ if ($dir eq '-');
						if ($ref_deg->{$indcodon2} eq 'stop'){ ++$stopcodons; push(@validalleles, 1) }
						else { push(@validalleles, 0) }
						} continue { ++$allind }
					$speciesstop+=$stopcodons;
					$validspecal[$metapop]+=scalar(@validalleles);
					my $freq=@validalleles ? $stopcodons/scalar(@validalleles) : 0;
					print "Population $subpop: Stop codons occur at a frequency of ", sprintf("%.4f", $freq), " (", scalar(@validalleles), " alleles covered).\n";
					if (@validalleles<2*$minind->[$subpop]){ $validgroup[$metapop]=0 }
					elsif ($validgroup[$metapop] && $freq>0){
						my $ref_suballeles=Misc::samplewithoutRep(\@validalleles, 2*$minind->[$subpop]);
						for my $allele (@$ref_suballeles){ ++$nder[$metapop] if ($allele) }
						}
					}
				$speciesfreq[$metapop]=sprintf( "%.4f", $speciesstop/$validspecal[$metapop]) if ($validspecal[$metapop]>0);
				++$self->{_SFS}[$metapop][ $nder[$metapop] ] if $validgroup[$metapop];
				}

			my $poppair=0;
			POP1: for my $pop1 ( 0..($#{$groups}-1) ){
				POP2: for my $pop2 ( ($pop1+1)..$#{$groups} ){
					next POP2 unless ($validgroup[$pop1] && $validgroup[$pop2]);
					++$self->{_jointSFS}[$poppair][ $nder[$pop1] ][ $nder[$pop2] ];
					} continue { ++$poppair }
				}
			my $freqstring=join(';', @validspecal);
			$freqstring.="\t" . join(';', @speciesfreq);
			$regions->setValue($chrom, $locus, 'info', $freqstring);
			}
		}
	my $loci=Regions->new(\%stop_positions);
	return($loci);
	}


sub printDeg{
	my ($self, $outfile, $regions, $refseq, $ref_genotypes)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file $outfile!\n" }

	my $ref_deg=Misc::getDegRef();
	my ($validcodons, $NinCodon, $notbiallelic, $twopoly, $undefdeg, $prematstop, $totalcodons)=(0, 0, 0, 0, 0, 0, 0);

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->getValue($chrom, $locus, 'start');
			my $end=$regions->getValue($chrom, $locus, 'end');
			my $frame=$regions->getValue($chrom, $locus, 'frame');
			my $dir=$regions->getValue($chrom, $locus, 'dir');
			my $locuslength=$end-$start+1;
			my $shift=$frame;
			my $clip=$frame;
			if ($dir eq '+'){ $clip=($locuslength-$frame) % 3 }
			elsif ($dir eq '-'){ $shift=($locuslength-$frame) % 3 }
			if ($locuslength-$shift-$clip<3){ print "No complete codon in $chrom:$start-$end, skipping locus!\n"; next LOCUS }
			my $refseq_string=$refseq->randomAccess($chrom, $start+$shift, $end-$clip, 1);
			if (length($$refseq_string) % 3){ warn "$chrom:", $start+$shift, "-", $end-$clip, " is not a multiple of 3!\n" }
			my $valid=0;

			my @codstarts;
			for (my $codstart=$start+$shift; $codstart+2<=$end; $codstart+=3){ push(@codstarts, $codstart) }
			if ($dir eq '-'){
				@codstarts=reverse(@codstarts);	# reverse the array of codon starting positions.
				$$refseq_string=reverse($$refseq_string);	# reverse the reference string.
				}

			CODON: for my $codstart(@codstarts){	# loop through the array of codon start positions
				++$totalcodons;
				my $codoffset=$codstart-$start;
				my $refcodon=substr($$refseq_string, 0, 3, "");	# chew up reference string by next three bases.
				$refcodon=uc($refcodon);
				if ($refcodon=~/N/){	# if one of the codon positions is hardmasked.
					++$NinCodon;
					next CODON;
					}
				my $basecodon=$refcodon;
				my $altcodon=$refcodon;
				my $poly=0;
				POS: for my $codonpos (0..2){
					my $pos=($dir eq '-') ? $codstart+(2-$codonpos) : $codstart+$codonpos;	# to make sure that reverse codons are processed backwards.
					if ( defined $ref_genotypes->getValue($chrom, $pos) ){
						my @alleles=$ref_genotypes->getAlleles($chrom, $pos);
						if ($poly && @alleles>1){ warn "Warning: Codon $chrom:$codstart-", $codstart+3, " has more than one polymorphic position!\n"; ++$twopoly; next CODON }
						if (@alleles>2){ warn "Warning: $chrom:$pos is not biallelic!\n"; ++$notbiallelic; next CODON }
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
				unless (defined $deg){ warn "Warning, degeneracy pattern for codon $basecodon at position $codstart, offset $codoffset is not defined! $frame, $dir\n"; ++$undefdeg; next CODON }
				$altcodon=~tr/ATGC/TACG/ if ($dir eq '-');
				my $altdeg=$ref_deg->{$altcodon};
				unless (defined $altdeg){ warn "Warning, degeneracy pattern for alternative codon $altcodon at position $codstart, offset $codoffset is not defined!\n"; ++$undefdeg; next CODON }

				if ($deg eq 'stop' || $altdeg eq 'stop'){
					my $string=($basecodon eq $altcodon) ? "$basecodon" : "$basecodon/$altcodon";
					if ($codstart!=$codstarts[-1]){
						warn "Warning, premature stop codon ($deg/$altdeg) detected at $chrom:$codstart-", $codstart+2, "! CDS end at $codstarts[-1].\n";
						++$prematstop;
						print $output "$chrom\t", $codstart, "\t", $codstart+3, "\tpremature stop codon: $string\n";
						next LOCUS;
						}
					else {
						print $output "$chrom\t", $codstart, "\t", $codstart+3, "\tstop codon: $string\n";
						}
					}
				else {
					for my $codonpos (0..2){
						my $pos=($dir eq '-') ? $codstart+(2-$codonpos) : $codstart+$codonpos;	# to make sure that reverse codons are processed backwards.
						my $degstate=substr($deg, 0, 1, "");
						my $altdegstate=substr($altdeg, 0, 1, "");
						if ($degstate eq $altdegstate){ print $output "$chrom\t$pos\t", $pos+1, "\t$degstate\n" }
						else { print $output "$chrom\t$pos\t", $pos+1, "\t$degstate/$altdegstate\n" }
						}
					}
				$valid=1;
				++$validcodons;
				} continue { 
					print $output "$chrom\t", $codstart, "\t", $codstart+3, "\tundefined\n" unless $valid;
					}
			}
		}
	print "Total number of codons assessed: $totalcodons\n";
	print "Validcodons: $validcodons\n";
	print "Codons with Ns: $NinCodon\n";
	print "Codons with non-biallelic positions: $notbiallelic\n";
	print "Codons with two or more polymorphic positions: $twopoly\n";
	print "Codons with undefined degenarcy pattern: $undefdeg\n";
	print "Premature stop codons: $prematstop\n";
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub printTranscripts{
	my ($self, $regions, $outfile, $append, $removeexcluded)=@_;
	my $output;

	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT }
	elsif ($append){ open $output, ">>", $outfile or die "Could not open output file!\n" }
	else { open $output, ">", $outfile or die "Could not open output file!\n" }

	SCAFF: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} keys %{ $regions->{_Loci} } ){
		LOCUS: for my $locus ( sort {$a->{'start'}<=>$b->{'start'}} @{ $regions->{_Loci}{$chrom} } ){
			if ( $removeexcluded && $locus->{'excluded'} ){ next LOCUS }
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'geneid'}\t$locus->{'transcriptid'}\t", $locus->{'bin'} // "", "\t", $locus->{'info'} || "", "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

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
			my @arraysize=( $#{ $self->{_jointSFS}[$poppair] }, $#{ $self->{_jointSFS}[$poppair][0] } );
			$self->{_jointSFS}[$poppair][0][0]+=$self->{_jointSFS}[$poppair][-1][-1];	# the top right bin needs to be combined with the lower left bin.
			$self->{_jointSFS}[$poppair][-1][-1]=0;
			for my $nder1 (0..$arraysize[0]){ print $output "\td${pop1}_${nder1}" }
			print $output "\n";
			for my $nder2 (0..$arraysize[1]){
				print $output "d${pop2}_${nder2}";
				for my $nder1 (0..$arraysize[0]){ print $output "\t$self->{_jointSFS}[$poppair][$nder1][$nder2]" }
				print $output "\n";
				}
			} continue { ++$poppair }
		}
	close $output;
	return;
	}




1;

