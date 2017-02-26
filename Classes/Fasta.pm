package Fasta;

use strict;
use warnings;
use Classes::Regions;
use Classes::Misc;


sub new{
	my $class=shift;
	my $seqs=shift;
	my $self={};
	$self->{_Seqs}=defined $seqs ? $seqs : {};
	bless $self, $class;
	return $self;
	}


sub getSortedIDs{
	my $self=shift;
	my @ids=( exists $self->{_Offsets} ) ? sort { $self->{_Offsets}{$a}<=>$self->{_Offsets}{$b} } keys %{ $self->{_Seqs} } : sort { Misc::expand($a) cmp Misc::expand($b) } keys %{ $self->{_Seqs} };
	return(@ids);
	}

sub getSequences{
	my $self=shift;
	return($self->{_Seqs});
	}

sub setSequences{
	my ($self, $ref_seqs)=@_;
	$self->{_Seqs}=$ref_seqs if defined($ref_seqs);
	return($self->{_Seqs});
	}

sub getSequence{
	my ($self, $id)=@_;
	unless (@_==2){ die "Wrong number of arguments!\n" }
	return($self->{_Seqs}{$id});
	}

sub setSequence{
	my ($self, $id, $ref_seq)=@_;
	unless (@_==3){ die "Wrong number of arguments!\n" }
	$self->{_Seqs}{$id}=$$ref_seq;
	return;
	}

sub getLength{
	my ($self, $id)=@_;
	unless (@_==2){ die "Wrong number of arguments!\n" }
	return($self->{_Length}{$id});
	}


sub readIndex{
	my ($self, $index_file, $fasta_file, $getseqs)=@_;

	$self->{_Filename}=$fasta_file;
	my $fasta;
	open $fasta, "<", $self->{_Filename} or die "Could not open fasta file $self->{_Filename}!\n" if ($getseqs);
	open my $index, "<", $index_file or die "Could not open fasta index $index_file!\n";

	while (<$index>){
		chomp;
		if (/^\w+\s/){
			my @entry=split("\t", $_);
			my $scaffold=$entry[0];
			$self->{_Lengths}{$scaffold}=$entry[1];
			$self->{_Offsets}{$scaffold}=$entry[2];
			if (defined $self->{_Rowlength} && $self->{_Rowlength} ne $entry[3]){warn "Rowlength not identical for each line in fasta file! $scaffold: $self->{_Rowlength} vs. $entry[3]\n"}
			elsif ( !defined $self->{_Rowlength} ){ $self->{_Rowlength}=$entry[3] }
			if (defined $self->{_Rowbytes} && $self->{_Rowbytes} ne $entry[4]){warn "Bytes per line not identical for each line in fasta file! $scaffold: $self->{_Rowbytes} vs. $entry[4]\n"}
			elsif ( !defined $self->{_Rowbytes} ){ $self->{_Rowbytes}=$entry[4] }
			if ($getseqs){
				my $newlinebytes=$self->{_Rowbytes}-$self->{_Rowlength};	# number of bytes in excess of characters per line, usually 1 for the newline character.
				my $nrows=int( $self->{_Lengths}{$scaffold} / $self->{_Rowlength} );	# number of rows for scaffold.
				my $bytelen=$self->{_Lengths}{$scaffold} + ( $nrows * $newlinebytes );
				unless ( seek($fasta, $self->{_Offsets}{$scaffold}, 0) ){ die "Could not set position in fasta file to byte offset $self->{_Offsets}{$scaffold}!\n" };
				unless ( read($fasta, $self->{_Seqs}{$scaffold}, $bytelen) ){ die "Could not extract $scaffold!\n" }
				$self->{_Seqs}{$scaffold}=~tr/\n//d;
				unless ( length( $self->{_Seqs}{$scaffold} )==$self->{_Lengths}{$scaffold} ){
					die "Wrong length of extracted sequence string for scaffold $scaffold! (", length( $self->{_Seqs}{$scaffold} ), " vs. $self->{_Lengths}{$scaffold})\n";
					}
				}
			else { $self->{_Seqs}{$scaffold}=1 }
			}
		}
	close $fasta if $getseqs;
	close $index;
	print STDERR "Fasta index $index_file successfully read in.\n";
	return;
	}


sub locifromIndex{
	my ($self, $sizelimit, @chromosomes)=@_;
	if (@chromosomes){
		@chromosomes=keys %{ $self->{_Seqs} } if ($chromosomes[0] eq 'all');
		}
	else { @chromosomes=keys %{ $self->{_Seqs} } }
	my %loci;

	CHROM: for my $chrom (@chromosomes){
		if ( defined $loci{$chrom} ){ warn "Sequence ID $chrom appreared twice in fasta index! Skipping sequence.\n"; next CHROM }
		my $seqlength=$self->{_Lengths}{$chrom};
		unless (defined $seqlength){ warn "Sequence ID $chrom not found in fasta index! Skipping sequence.\n"; next CHROM }
		my $maxsize=(defined $sizelimit) ? $sizelimit : $seqlength;
		my $start;
		my $end=0;
		while ($end<$seqlength){
			$start=$end+1;
			$end+=$maxsize;
			$end=$seqlength if ($end>$seqlength);
			push @{$loci{$chrom}}, {'start'=>$start-1, 'end'=>$end-1};
			}
		}

	my $regions=Regions->new(\%loci);
	return($regions);
	}


sub randomAccess{
	unless (@_>=4){ die "Wrong number of arguments!\n" }
	my ($self, $id, $start, $end, $rmnewline)=@_;
	if ($end>$self->{_Lengths}{$id}-1){ warn "Invalid end position for scaffold $id!\n"; $end=$self->{_Lengths}{$id}-1 }
	my $nrowsstart=int($start/$self->{_Rowlength});
	my $nrowsend=int($end/$self->{_Rowlength});
	my $postart=$start % $self->{_Rowlength};
	my $posend=$end % $self->{_Rowlength};
	my $startbyte=$self->{_Offsets}{$id}+$nrowsstart*$self->{_Rowbytes}+$postart;
	my $endbyte=$self->{_Offsets}{$id}+$nrowsend*$self->{_Rowbytes}+$posend;
	my $bytelen=$endbyte-$startbyte+1;
	my $sequence;
	open my $fasta, "<", $self->{_Filename} or die "Could not open fasta file $self->{_Filename}!\n";
	unless ( seek($fasta, $startbyte, 0) ){ die "Could not set position in fasta file to $startbyte!\n" };
	unless ( read($fasta, $sequence, $bytelen) ){ die "Could not extract $id:$start-$end!\n" }
#	print STDERR "Sequence $id:$start-$end successfully extracted.\n";
	close $fasta;
	$sequence=~tr/\n//d if $rmnewline;
	return(\$sequence);
	}


sub readSequences{
	my ($self, $filename, $maxseq)=@_;

	$maxseq=100000 unless (defined $maxseq);

	open my $fasta, "<", $filename or die "Could not open fasta file $filename!\n";

	my %seqs;
	my $nseq=0;
	my $flag=0;
	my $seqtitle;

	READ: while (my $line= <$fasta>){
		chomp $line;

		if ($line=~/^\s*#/ || length($line)==0){next READ}

		elsif ($line=~/^\s*\>/){
			++$nseq;
			if ($nseq>$maxseq){last READ} else {$flag=1}
			if ($line=~/\>([^\n\s]+)/){$seqtitle=$1} else {$seqtitle="seq_".$nseq}
			$seqs{$seqtitle}="";
			next;
			}

		elsif ($flag==1){ $seqs{$seqtitle}.=$line }
		}

	close $fasta;
	$self->{_Seqs}=\%seqs;
	%{ $self->{_Lengths} }=map { $_ => length($seqs{$_}) } keys %seqs;
	return;
	}

sub readSequencesList{
	my ($self, $ref_filelist, $maxseq)=@_;

	die "Must provide reference to array of filenames!\n" unless ( ref($ref_filelist) eq "ARRAY" );
	$maxseq=100000 unless (defined $maxseq);
	my %seqs;
	my $nseq=0;

	for my $filename (@$ref_filelist){
		open my $fasta, "<", $filename or die "Could not open fasta file $filename!\n";
		my $flag=0;
		my $seqtitle;
		READ: while (my $line= <$fasta>){
			chomp $line;
			if ($line=~/^\s*#/ || length($line)==0){ next READ }
			elsif ($line=~/^\s*\>/){
				++$nseq;
				if ($nseq>$maxseq){ last READ } else { $flag=1 }
				if ($line=~/\>([^\n\s]+)/){ $seqtitle=$1 } else { $seqtitle="seq_".$nseq }
				$seqs{$seqtitle}="";
				next;
				}
			elsif ($flag==1){ $seqs{$seqtitle}.=$line }
			}
		close $fasta;
		}

	$self->{_Seqs}=\%seqs;
	return;
	}


sub concatSeq{
	my ($self, $seq_name, $ref_seq)=@_;
	$self->{_Seqs}{$seq_name}.=$$ref_seq;
	$self->{_Lengths}{$seq_name}=length($self->{_Seqs}{$seq_name});	# update sequence length information.
	return;
	}


sub extractSeqs{
	my ($self, $ref_names)=@_;
	unless (@_==2){ die "Wrong number of arguments!\n" }
	my $extracted=Fasta->new();

	for my $sel (@{$ref_names}){
		if (exists $self->{_Seqs}{$sel}){
			$extracted->{_Seqs}{$sel}=$self->{_Seqs}{$sel};
			$extracted->{_Lengths}{$sel}=$self->{_Lengths}{$sel};
			$extracted->{_Offsets}{$sel}=$self->{_Offsets}{$sel};
			}
		else { warn "Sequence $sel does not exist in fasta file!\n" }
		}

	$extracted->{_Rowlength}=$self->{_Rowlength};
	$extracted->{_Rowbytes}=$self->{_Rowbytes};
	$extracted->{_Filename}=$self->{_Filename};
	return($extracted);
	}


sub extractRegions{
	my ($self, $regions)=@_;

	my %selregs;

	for my $seq ($regions->allKeys()){
		for my $reg ($regions->allValues($seq)){
			my $lower=$reg->{'start'};
			my $upper=$reg->{'end'};
			my $length=$upper-$lower+1;
			my $sequence=substr($self->{_Seqs}{$seq}, $lower, $length);
			my $newname=$seq . "_start_" . $lower . "_end_" . $upper;
			$selregs{$newname}=$sequence;
#			push(@{$selregs{$seq}}, {'start'=>$lower, 'end'=>$upper, 'seq'=>$sequence});
			}
		}

	my $extracted=Fasta->new(\%selregs);
	return($extracted);
	}


sub stackSequences{
	my ($self, $min_length, $locus_length)=@_;

	$min_length=$min_length>$locus_length ? $min_length : $locus_length;
	my $totlength=0;
	my $start=0;
	my $end=0;
	my %borders;
	my @excluded;

	for my $key (sort keys %{$self->{_Seqs}}){
		my $seq_length;
		if (defined $self->{_Lengths}{$key}){$seq_length=$self->{_Lengths}{$key}}
		else {$seq_length=length($self->{_Seqs}{$key})}
		if ($seq_length==1){die "Wrong information of sequence length provided!\n"}
		elsif($seq_length>=$min_length){
			$totlength+=($seq_length-$locus_length);
			$end=($start+$seq_length-$locus_length-1);
			push @{$borders{$key}}, {'start'=>$start, 'end'=>$end, 'offset'=>0};
			$start=$end+1;
			}
		else{push(@excluded, $key)}
		}

	$self->{_StackLength}=$totlength;
	$self->{_Borders}=\%borders;
	$self->{_Excluded}=\@excluded;
	return;
	}


sub stackRegions{
	my ($self, $ref_regions, $min_length, $locus_length)=@_;

	$min_length=$min_length>$locus_length ? $min_length : $locus_length;
	my $totlength=0;
	my $start=0;
	my $end=0;
	my %borders;
	my @excluded;

	for my $key (sort $ref_regions->allKeys()){
		for my $reg ($ref_regions->allValues($key)){
			my $length=$reg->{'end'}-$reg->{'start'}+1;
			if($length>=$min_length){
				$totlength+=($length-$locus_length);
				$end=($start+$length-$locus_length-1);
				push @{$borders{$key}}, {'start'=>$start, 'end'=>$end, 'offset'=>$reg->{'start'}};
				$start=$end+1;
				}
			else{push(@excluded, {'id'=>$key, 'start'=>$reg->{'start'}, 'end'=>$reg->{'end'}})}
			}
		}
	$self->{_StackLength}=$totlength;
	$self->{_Borders}=\%borders;
	$self->{_Excluded}=\@excluded;
	return;
	}


sub printBorders{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	for my $key (sort keys %{$self->{_Borders}}){
		for my $reg (@{$self->{_Borders}{$key}}){
			print $output "$key:\t", length($self->{_Seqs}{$key}), "\tStart: $reg->{'start'}\tEnd: $reg->{'end'}\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub sampleLoci{
	my ($self, $nloci, $locus_length, $prop)=@_;

	my @loci;

	while (@loci<$nloci){
		my @random=sort {$a<=>$b} (map { int(rand($self->{_StackLength})) } (1..$nloci));
		my ($remove, $index);

		SEQ: for my $key (sort keys %{$self->{_Borders}}){
			REG: for my $reg (@{$self->{_Borders}{$key}}){
				$remove=0; $index=0;

				LOCUS: for my $pos (@random){
					if ($pos<$reg->{'start'}){ ++$index; next LOCUS }
					elsif ($pos>$reg->{'end'}){ next REG }
					else {
						my $start=$pos-$reg->{'start'}+$reg->{'offset'};
						my $end=$start+$locus_length-1;
						print STDERR "Position $pos assigned to scaffold $key: $start-$end.\n";
						my $seq;
						if ( defined $self->{_Filename} ){
							$seq=${ $self->randomAccess($key, $start, $end, 1) };
							unless (length($seq)==$locus_length){
								die "Wrong length of referenece sequence for $key:$start-$end (", length($seq), " vs. $locus_length)!\n";
								}
							}
						my $seqcopy=$seq;
						my $N=($seqcopy=~tr/Nn//);
						print STDERR "$N Ns in reference sequence, ", $locus_length*$prop, " maximum acceptable!\n";
						unless ($N>$locus_length*$prop){ push @loci, {'id'=>$key, 'start'=>$start, 'end'=>$end, 'refseq'=>$seq } }
						++$remove;
						}
					}
				} continue {
					if ($remove>0){ splice(@random, $index, $remove) }
					if (scalar(@random)==0){ last SEQ }
					}
			}
		}

	my $candidates=Candidates->new([ @loci[0..($nloci-1)] ]);
	return($candidates);
	}


sub maskFasta{
	my ($self, $ref_regions, $mode)=@_;
	unless ($mode eq 'hard' || $mode eq 'soft'){ die "Mode $mode not recognized!\n" }
	my $totmasked=0;
	my $totlength=0;

	CHROM: for my $chrom ( $ref_regions->allKeys() ){
		unless ( exists $self->{_Seqs}{$chrom} ){
			warn "Scaffold $chrom does not exist in reference genome, skipping all intervals for scaffold $chrom!\n";
			next CHROM;
			}
		my $chrommasked=0;
		LOCUS: for my $locus ( $ref_regions->allKeys($chrom) ){
			my $start=$ref_regions->getValue($chrom, $locus, 'start');
			my $end=$ref_regions->getValue($chrom, $locus, 'end');
			unless ( $end < $self->{_Lengths}{$chrom} ){
				if ( $start < $self->{_Lengths}{$chrom} ){
					warn "Interval $chrom:$start-$end is partially outside of reference sequence, adjusting end of interval!\n";
					$end=$self->{_Lengths}{$chrom}-1;
					}
				else {
					warn "Interval $chrom:$start-$end is completely outside of reference sequence, skipping interval!\n";
					next LOCUS;
					}
				}
			my $locuslength=$end-$start+1;

			if ($mode eq 'hard'){
				substr($self->{_Seqs}{$chrom}, $start, $locuslength, "N" x $locuslength);
				}
			else {
				my $replacement=lc( substr($self->{_Seqs}{$chrom}, $start, $locuslength) );
				substr($self->{_Seqs}{$chrom}, $start, $locuslength, $replacement);
				}
			$chrommasked+=$locuslength;
			}
		$totmasked+=$chrommasked;
		$totlength+=$self->{_Lengths}{$chrom};
		print STDERR "Masked $chrommasked of $self->{_Lengths}{$chrom} (", sprintf("%.2f", $chrommasked/$self->{_Lengths}{$chrom}*100), "%) sites for scaffold $chrom.\n";
		}
	print STDERR "Masked a total of $totmasked of $totlength (", sprintf("%.2f", $totmasked/$totlength*100), "%) sites for the entire reference genome.\n";
	return;
	}

sub maskCpG{
	my ($self, $mode, $includeN)=@_;
	unless ($mode eq 'hard' || $mode eq 'soft'){ die "Mode $mode not recognized!\n" }
	my $totmasked=0;
	my $totlength=0;

	ID: for my $id ( keys %{ $self->{_Seqs} } ){
		my $idmasked=0;
		if ($mode eq 'hard'){
			$idmasked+=2*($self->{_Seqs}{$id}=~s/CA/NN/ig);
			$idmasked+=2*($self->{_Seqs}{$id}=~s/TG/NN/ig);
			$idmasked+=2*($self->{_Seqs}{$id}=~s/CG/NN/ig);
			if ($includeN){
				$idmasked+=($self->{_Seqs}{$id}=~s/CN/NN/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/NA/NN/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/TN/NN/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/NG/NN/ig);
				}
			}
		else{
			$idmasked+=2*($self->{_Seqs}{$id}=~s/CA/ca/ig);
			$idmasked+=2*($self->{_Seqs}{$id}=~s/TG/tg/ig);
			$idmasked+=2*($self->{_Seqs}{$id}=~s/CG/cg/ig);
			if ($includeN){
				$idmasked+=($self->{_Seqs}{$id}=~s/CN/cN/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/NA/Na/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/TN/tN/ig);
				$idmasked+=($self->{_Seqs}{$id}=~s/NG/Ng/ig);
				}
			}

		$totmasked+=$idmasked;
		$totlength+=$self->{_Lengths}{$id};
		print STDERR "Masked $idmasked of $self->{_Lengths}{$id} (", sprintf("%.2f", $idmasked/$self->{_Lengths}{$id}*100), "%) sites for scaffold $id.\n";
		}
	print STDERR "Masked a total of $totmasked of $totlength (", sprintf("%.2f", $totmasked/$totlength*100), "%) sites for the entire reference genome.\n";
	return;
	}


sub generateAncRef{
	my ($self, $regions, $ref_genotypes, $ref_missing, $ref_minind, $groups, $mincov, $excluderefN)=@_;

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

	ID: for my $id ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		my $prevend;
		LOCUS: for my $locus ( $regions->allKeys($id) ){
			my $start=$regions->{_Loci}{$id}[$locus]{'start'};
			my $end=$regions->{_Loci}{$id}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$self->randomAccess($id, $start, $end, 1);
			if ( length($$ref_seqstring) != $locuslength ){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my ($totalpos, $hardmasked, $filtered, $undefanc)=(0, 0, 0, 0);
			if (defined $prevend && $start-$prevend != 1){
				die "Gap/overlap detected between two consecutive loci on the same scaffold/chromosome!\n";
				}
			my ($refbase, $ancbase);
			my $seq="";

			POSITION: for my $pos ($start..$end){
				$refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				$ancbase='N';
				my $offset=$pos-$start;
				if ($excluderefN && $refbase eq 'N'){ ++$hardmasked; next POSITION }	# discard if base in reference sequence is hardmasked.
				my $snp=defined $ref_genotypes->getValue($id, $pos) ? 1 : 0;
				my @alleles=(0) x scalar(@$groups);

				GROUP: for my $group ( 0..$#{$groups} ){
					my $groupgts="";
					for my $subpop ( @{$groups->[$group]} ){
						my $validind=0;
						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($id, $start, $allind, $offset);
							++$validind if (defined $coverage && $coverage>=$mincov);
							}

						if ($validind<$ref_minind->[$subpop]){ ++$filtered; next POSITION }	# not enough individuals callable in a population.

						if ($snp){
							my $gtstring=$ref_genotypes->getValue($id, $pos, $subpop);
							my $nalleles=length($gtstring);
							my $nmissing=($gtstring=~tr/\.//);
							if ( ( ($nalleles-$nmissing)/2 ) < $ref_minind->[$subpop] ){ ++$filtered; next POSITION }	# not enough individuals callable in a population.
							$groupgts.=$gtstring;
							}
						}
					if ($snp && $groupgts){
						for my $allele (0..3){ $alleles[$group] |= (1 << $allele) if ($groupgts=~/$allele/) }
						}
					}
				unless ($snp){ $ancbase=$refbase; next POSITION }	# Position not polymorphic but sufficently covered.

#				my $all=($alleles[0] | $alleles[1]) | $alleles[2];
#				my ($nalleles, $outstate)=Misc::countbits2($all);
#				if ($nalleles>2){ next POSITION }

				my ($nalleles01, $outstate01)=Misc::countbits2($alleles[0] | $alleles[1]);
				if ($nalleles01==1){ $ancbase=$ref_genotypes->getBase($id, $pos, $outstate01); next POSITION }	# Ingroups plus parva monomorphic.

				my ($nalleles12, $outstate12)=Misc::countbits2($alleles[1] | $alleles[2]);
				if ($nalleles12==1){ $ancbase=$ref_genotypes->getBase($id, $pos, $outstate12); next POSITION }	# Both outgroups monomorphic.

				my ($nalleles02, $outstate02)=Misc::countbits2($alleles[0] | $alleles[2]);
				if ($nalleles02==1){ $ancbase=$ref_genotypes->getBase($id, $pos, $outstate02); next POSITION }	# Ingroups plus HYP monomorphic.

				++$undefanc;	# no valid allele constellation for ancestral state reconstruction.
				} continue { $seq.=$ancbase; ++$totalpos }

			unless ( length($seq)==$locuslength ){ die "Wrong length of sequence string $id, ", length($seq), " vs. $locuslength!\n" }
			$self->{_Ancestral}{$id}.=$seq;
			$prevend=$end;
			print STDERR "$id:$start-$end, locus length: $locuslength, positions considered: $totalpos, hardmasked: $hardmasked, filtered: $filtered, undefined ancestral state: $undefanc.\n";
			}
		}
	return;
	}


sub fastaStats{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file $outfile!\n" }

	my @totcounts=(0) x 6;
	my $totlength=0;

	for my $title ( sort {Misc::expand($a) cmp Misc::expand($b)} keys %{ $self->{_Seqs} } ){
		my @counts=(0) x 6;
		my $seq=$self->{_Seqs}{$title};	# make hardcopy of sequence string.		
		my $seqlength=length($seq);
		$totlength+=$seqlength;
		$counts[0]=($seq=~tr/Gg//);
		$counts[1]=($seq=~tr/Cc//);
		$counts[2]=($seq=~tr/Aa//);
		$counts[3]=($seq=~tr/Tt//);
		$counts[4]=($seq=~tr/Nn//);
		$counts[5]=($seq=~tr/WwRrMmKkYySs//);
		unless (Misc::sum(@counts)==$seqlength){ warn "Unidentified bases in sequence string of scaffold $title!\n" }
		$totcounts[$_]+=$counts[$_] foreach (0..$#counts);
		my $validsites=Misc::sum(@counts[0..3]);

		print $output "Scaffold $title:\n",
			"Length: $seqlength bp\n",
			"G: ", sprintf("%.2f%%", $counts[0]/$seqlength*100), "\n",
			"C: ", sprintf("%.2f%%", $counts[1]/$seqlength*100), "\n",
			"A: ", sprintf("%.2f%%", $counts[2]/$seqlength*100), "\n",
			"T: ", sprintf("%.2f%%", $counts[3]/$seqlength*100), "\n",
			"N: ", sprintf("%.2f%%", $counts[4]/$seqlength*100), "\n";
		print $output "GC content: ", sprintf("%.2f%%", ($counts[0]+$counts[1])/$validsites*100), "\n" if ($validsites);
		print $output "Site heterozygosity: ", sprintf("%.6f", $counts[5]/$validsites), "\n" if ($validsites+$counts[5]);
		print $output "\n";
		}

	unless (Misc::sum(@totcounts)==$totlength){ warn "Total length of sequences does not match sum of base counts!\n" }
	my $totvalidsites=Misc::sum(@totcounts[0..3]);

	print $output "All sequences in fasta file:\n",
		"Length: $totlength bp\n",
		"G: ", sprintf("%.2f%%", $totcounts[0]/$totlength*100), "\n",
		"C: ", sprintf("%.2f%%", $totcounts[1]/$totlength*100), "\n",
		"A: ", sprintf("%.2f%%", $totcounts[2]/$totlength*100), "\n",
		"T: ", sprintf("%.2f%%", $totcounts[3]/$totlength*100), "\n",
		"N: ", sprintf("%.2f%%", $totcounts[4]/$totlength*100), "\n";
	print $output "GC content: ", sprintf("%.2f%%", ($totcounts[0]+$totcounts[1])/$totvalidsites*100), "\n" if ($totvalidsites);
	print $output "Site heterozygosity: ", sprintf("%.6f", $totcounts[5]/$totvalidsites), "\n" if ($totvalidsites+$totcounts[5]);
	print $output "\n";

	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub recodeVCF{
	my ($self, $vcffile, $outfile, $verbose)=@_;

	open my $input, "<", $vcffile or die "Could not open vcf file $vcffile!\n$!\n";
	open my $output, ">", $outfile or die "Could not open output vcf file $outfile!\n";
	print "Recoding vcf file $vcffile...\n";
	my ($prevpos, $prevscaffold);
	my $ref_ancseq;
	my %bitflags=(1 => "2MISSING", 2 => "1MISSING", 3 => "", 4 => "2ALLELES", 5 => "3+ALLELES", 6 => "1POLYM", 7 => "2POLYM", 8 => "NOTDEFINED");

	LINE: while (<$input>){
		if (/^#/){	# row belongs to the header.
			if(/^##FILTER/){
				print $output $_;
				print $output '##FILTER=<ID=NOTDEFINED,Description="Ancestral state not defined">', "\n";
				}
			elsif(/^##INFO=<ID=MLEAF/){
				print $output $_;
				print $output '##INFO=<ID=RECODING,Number=A,Type=String,Description=\"Recoding of reference allele">', "\n";
				}
			elsif(/^#/){ print $output $_ }
			}

		elsif(/^\w+\s/){	# SNP entries start here.
			chomp;
			my @line=split(/\s+/, $_);
			my $scaffold=$line[0];
			my $pos=$line[1];
			my @bases=($line[3], split(',' , $line[4]) );
			my $nalleles=scalar(@bases);

			if (!defined $prevscaffold || $prevscaffold ne $scaffold){
				print STDERR "$scaffold, $self->{_Lengths}{$scaffold}\n";
				$prevscaffold=$scaffold; $prevpos=1;
				unless ($self->{_Lengths}{$scaffold}){ warn "Scaffold $scaffold is not included in ancestral reference genome! Skipping scaffold.\n"; next LINE }
				$ref_ancseq=$self->randomAccess($scaffold, 0, $self->{_Lengths}{$scaffold}-1, 1);
				die "Wrong length of extracted sequence string!\n" unless ( length($$ref_ancseq)==$self->{_Lengths}{$scaffold} );
				}
			elsif (!$self->{_Lengths}{$scaffold}){ next LINE }

			substr($$ref_ancseq, 0, $pos-$prevpos, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
			my $ancbase=uc( substr($$ref_ancseq, 0, 1) );
			$prevpos=$pos;

			unless (defined $ancbase && $ancbase ne "N"){	# if outstate is not defined for site.
				$line[3]="N";	# set ancestral base to N.
				$line[4]=join(',', @bases);	# set derived alleles to the complete set of alleles.
				Misc::setFilter(\@line, 128, \%bitflags);
				Misc::setMissing(\@line);	# set all individual genotypes to missing.
				print $output join("\t", @line), "\n";
				warn "Outstate not defined for $scaffold:$pos, set all genotypes to missing!\n" if $verbose;
				next LINE;
				}

			my $outstate;
			for my $index (0..$#bases){
				if ($bases[$index] eq $ancbase){
					$outstate=$index;
					last;
					}
				}
			unless (defined $outstate){
				warn "Ancestral base $ancbase does not occur in list of bases at $scaffold:$pos!\n";
				push @bases, $ancbase;	# add ancestral base to array of bases at this SNP.
				$outstate=$#bases;	# after push, index of ancestral base is the index of last element in bases array.
				}

			my %translation;
			Misc::setFilter(\@line, 0, \%bitflags);	# set filter string to PASS.

			if ($outstate==0){ print $output join("\t", @line), "\n"; next LINE }	# reference base equals ancestral base, leave line unchanged.
			elsif ($outstate<3 && $nalleles<4){	# more than 3 alleles are currently not recoded.
				$line[3]=$ancbase;
				my @altalleles;
				for my $base (@bases){ push @altalleles, $base if ($base ne $ancbase) }
				$line[4]=join(',', @altalleles);
				my @fields=split(':', $line[8]);
				if (@fields<2){	# if no coverage at all for this site, finish line here already.
					if ($line[7] eq "."){ $line[7]="RECODING=" . $bases[0] . "->" . $ancbase }
					else { $line[7].=";RECODING=" . $bases[0] . "->" . $ancbase }
					print $output join("\t", @line), "\n";
					next LINE;
					}
				%translation=('0'=>$outstate, $outstate=>'0');
				}
			else {
				$line[3]="N";	# set ancestral base to N.
				$line[4]=join(',', @bases);	# set derived alleles to the complete set of alleles.
				Misc::setMissing(\@line);	# set all individual genotypes to missing.
				print $output join("\t", @line), "\n";
				next LINE;
				}

			my ($MLEAC1, $MLEAC2, $MLEAF1, $MLEAF2)=$line[7]=~/(?:MLE)?AC=(\d+),?(\d+)?;(?:MLE)?AF=(\d+\.?\d*),?(\d+\.?\d*)?/;
			my $nind=0;

			IND: for my $col (9..$#line){
				if ($line[$col]=~/([0-2])\/([0-2])/){
					my @call=split(':', $line[$col]);
					my @gt=split('/', $call[0]);
					my @pl=split(',', $call[3]);
					my $al1=exists $translation{$gt[0]} ? $translation{$gt[0]} : $gt[0];
					my $al2=exists $translation{$gt[1]} ? $translation{$gt[1]} : $gt[1];
					$call[0]=join('/', sort($al1, $al2) );	#	recode alleles with translation hash.
					if (@pl==3){
						$call[3]=join( ',', reverse(@pl) );	# flip order of genotype likelihood fields.
						}
					elsif (@pl==6){
						if ($outstate==1){
							$call[3]=join( ',', ( reverse(@pl[0..2]), $pl[4], $pl[3], $pl[5] ) );	# flip order of genotype likelihood fields.
							}
						elsif ($outstate==2){
							$call[3]=join( ',', ($pl[5], $pl[4], $pl[2], $pl[3], $pl[1], $pl[0]) );	# flip order of genotype likelihood fields.
							}
						else { die "Illegal outstate $outstate in $scaffold:$pos!\n" }
						}
					else { die "Wrongly formated genotype likelihood field for $scaffold:$pos, individual ", $col-8, "!\n" }

					$line[$col]=join(':', @call);	# reassemble genotype column.
					++$nind;
					}
				else { $line[$col]="./." }	# missing data in genotype field.
				}

			if ($nalleles==2){
				$line[7]="MLEAC=" . (2*$nind-$MLEAC1) . ";MLEAF=" . sprintf("%.3f", 1-$MLEAF1);	# recode allele frequency annotation.
				}
			elsif ($nalleles==3){
				if ($outstate==1){
					$line[7]="MLEAC=" . (2*$nind-$MLEAC1-$MLEAC2) . "," . $MLEAC2 . ";MLEAF=" . sprintf("%.3f", 1-$MLEAF1-$MLEAF2) . "," . $MLEAF2;	# recode allele frequency annotation.
					}
				elsif ($outstate==2){
					$line[7]="MLEAC=" . $MLEAC1 . "," . (2*$nind-$MLEAC1-$MLEAC2) . ";MLEAF=" . $MLEAF1 . "," . sprintf("%.3f", 1-$MLEAF1-$MLEAF2);	# recode allele frequency annotation.
					}
				else { die "Illegal outstate $outstate in $scaffold:$pos!\n" }
				}
			$line[7].=";RECODING=" . $bases[0] . "->" . $ancbase;
			print $output join("\t", @line), "\n";	# print modified line.
			}
		}
	close $input;
	close $output;
	return;
	}


sub printSeqs{
	my ($self, $outfile, $rowlength)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	for my $key (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Seqs}}){
		print $output ">$key\n";
		my $pos=0;
		while ($pos<length($self->{_Seqs}{$key})){
			print $output substr($self->{_Seqs}{$key}, $pos, $rowlength), "\n";
			$pos+=$rowlength;
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub printFasta{
	my ($self, $outfile, $rowlength, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file $outfile!\n" }

	my @titles=($sort && exists $self->{_Offsets}) ? sort { $self->{_Offsets}{$a}<=>$self->{_Offsets}{$b} } keys %{ $self->{_Seqs} } : sort {Misc::expand($a) cmp Misc::expand($b)} keys %{ $self->{_Seqs} };
	for my $title (@titles){
		print STDERR "Printing scaffold $title.\n";
		print $output ">$title\n";
		my $seq=$self->{_Seqs}{$title};	# make hardcopy of sequence string.
		unless ( length($seq)==$self->{_Lengths}{$title} ){ die "Length of sequence string $title does not match with fasta index anymore!\n" }
		while ( length($seq)>0 ){
			print $output substr($seq, 0, $rowlength, ""), "\n";	# chew up sequence string until empty.
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub printAncestral{
	my ($self, $outfile, $rowlength)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file $outfile!\n" }

	for my $title ( sort { $self->{_Offsets}{$a}<=>$self->{_Offsets}{$b} } keys %{ $self->{_Ancestral} } ){
		print STDERR "Printing scaffold $title.\n";
		print $output ">$title\n";
		my $seq=$self->{_Ancestral}{$title};	# make hardcopy of sequence string.
		unless ( length($seq)==$self->{_Lengths}{$title} ){ die "Length of sequence string $title does not match with fasta index anymore!\n" }
		while ( length($seq)>0 ){
			print $output substr($seq, 0, $rowlength, ""), "\n";	# chew up sequence string until empty.
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


1;

