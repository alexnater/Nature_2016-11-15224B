package Regions;

use strict;
use warnings;
use Classes::Misc;
use List::Util qw( min );


sub new{
	my $class=shift;
	my $loci=shift;
	my $self={};
	$self->{_Loci}=defined $loci ? $loci : {};
	bless $self, $class;
	return $self;
	}


sub getValue{
	my ($self, $chrom, $locus, $key)=@_;

	if (@_==1){return($self->{_Loci})}
	elsif (@_==2){return($self->{_Loci}{$chrom})}
	elsif (@_==3){return($self->{_Loci}{$chrom}[$locus])}
	elsif (@_==4){return($self->{_Loci}{$chrom}[$locus]{$key})}
	else {die "Invalid number of arguments!\n"}
	}


sub setValue{
	my $self=shift;
	unless (@_==4){die "Wrong number of arguments!\n"}
	my ($chrom, $locus, $key, $value)=@_;

	$self->{_Loci}{$chrom}[$locus]{$key}=$value;

	return($self->{_Loci}{$chrom}[$locus]{$key});
	}


sub getExcluded{
	my $self=shift;
	my ($chrom, $locus)=@_;

	if (@_==2){ return( $self->{_Loci}{$chrom}[$locus]{'excluded'} ) }
	else { die "Invalid number of arguments!\n" }
	}


sub setExcluded{
	my $self=shift;
	unless (@_==3){ die "Wrong number of arguments!\n" }
	my ($chrom, $locus, $value)=@_;
	unless ($value==0 || $value==1){ die "Invalid value, 0 or 1 required!\n" }

	$self->{_Loci}{$chrom}[$locus]{'excluded'}=$value;

	return( $self->{_Loci}{$chrom}[$locus]{'excluded'} );
	}


sub getSize{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(scalar(keys %{$self->{_Loci}}))}
	elsif (@_==2){return(scalar @{$self->{_Loci}{$chrom}})}
	elsif (@_==3){return(scalar(keys %{$self->{_Loci}{$chrom}[$locus]}))}
	else {die "Invalid number of arguments!\n"}
	}


sub allKeys{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(keys %{$self->{_Loci}})}
	elsif (@_==2){return((0..$#{$self->{_Loci}{$chrom}}))}
	elsif (@_==3){return(keys %{$self->{_Loci}{$chrom}[$locus]})}
	else {die "Invalid number of arguments!\n"}
	}


sub allValues{
	my ($self, $chrom, $locus)=@_;

	if (@_==1){return(values %{$self->{_Loci}})}
	elsif (@_==2){return(@{$self->{_Loci}{$chrom}})}
	elsif (@_==3){return(values %{$self->{_Loci}{$chrom}[$locus]})}
	else {die "Invalid number of arguments!\n"}
	}


sub findLocus{
	my ($self, $id, $pos)=@_;
	unless (@_==3){die "Wrong number of arguments!\n"}

	my $index;
	LOCUS: for my $locus (sort {$self->{_Loci}{$id}[$a]{'start'}<=>$self->{_Loci}{$id}[$b]{'start'}} $self->allKeys($id)){
		if ($pos>=$self->{_Loci}{$id}[$locus]{'start'} && $pos<=$self->{_Loci}{$id}[$locus]{'end'}){
			$index=$locus;
			last LOCUS;
			}
		elsif ($pos>$self->{_Loci}{$id}[$locus]{'end'}){next LOCUS}
		else {die "Could not find locus index for $id:$pos!\n"}
		}
	return $index;
	}


sub getreftoIndividuals{
	my ($self, $mode, $chrom, $id)=@_;
	if ($mode eq 'position'){return($self->{_PosSummary}{$chrom}{$id})}
	elsif ($mode eq 'locus'){return($self->{_IndSummary}{$chrom}[$id])}
	else {die "Mode not recognized!\n"}
	}


sub nextIterator{
	my $self=shift;
	return(each %{$self->{_Loci}});
	}


sub resetIterator{
	my $self=shift;
	keys %{$self->{_Loci}};
	return;
	}


sub addLocus{
	my ($self, %args)=@_;
	unless ( defined $args{'id'} && defined $args{'start'} && defined $args{'end'} ){ die "Locus location arguments missing!\n" }
	push @{ $self->{_Loci}{ $args{'id'} } }, { 'start'=>$args{'start'}, 'end'=>$args{'end'} };
	if ( defined $args{'refseq'} ){ $self->{_Loci}{ $args{'id'} }[-1]{'refseq'}=$args{'refseq'} }
	if ( defined $args{'info'} ){ $self->{_Loci}{ $args{'id'} }[-1]{'info'}=$args{'info'} }
	return;
	}


sub addRegions{
	my ($self, $ref_regions, $sort)=@_;
	my $lociadded=0;
	keys %{ $ref_regions->{_Loci} };
	CHROM: while ( my ($chrom, $ref_chrom)=each %{ $ref_regions->{_Loci} } ){
		if ( defined $self->{_Loci}{$chrom} ){ push @{ $self->{_Loci}{$chrom} }, @$ref_chrom }
		else { $self->{_Loci}{$chrom}=$ref_chrom }
		$lociadded+=scalar(@{ $ref_chrom });
		if ($sort){	# sort locus array by start position after adding new loci.
			my @sorted=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} };
			$self->{_Loci}{$chrom}=\@sorted;
			}
		}
	print STDERR "Added a total of $lociadded regions to the object.\n";
	return($lociadded);
	}


sub addRefseq{
	my ($self, $fasta)=@_;
	ID: while ( my ($id, $ref_id)=each %{ $self->{_Loci} } ){
		LOCUS: for my $ref_locus (@$ref_id){
			my $ref_seqstring=$fasta->randomAccess($id, $ref_locus->{'start'}, $ref_locus->{'end'}, 1);
			$ref_locus->{'refseq'}=$$ref_seqstring;
			print STDERR "Added reference sequence to $id:$ref_locus->{'start'}-$ref_locus->{'end'}.\n";
#			print STDERR "$ref_locus->{'refseq'}\n";
			}
		}
	return;
	}


sub readBED{
	my ($self, $filename, $format, $keepinfo)=@_;
	unless(@_==4){die "Wrong number of arguments!\n"}
	if ($filename eq "none"){ return }

	my %regions;

	open my $bed, "<", $filename or die "Could not open bed file $filename!\n";

	while (my $line= <$bed>){
		chomp $line;
		if ($line=~/^(\w+)\t(\d+)\t(\d+)\t?(.*)/){
			if ($2>$3){ warn "$1:$2-$3 has smaller end than start coordinate!\n"; next }
			if ($format eq 'bed'){push(@{$regions{$1}}, {'start'=>$2, 'end'=>($3-1)})}
			elsif ($format eq 'onebased'){push(@{$regions{$1}}, {'start'=>($2-1), 'end'=>($3-1)})}
			else {push(@{$regions{$1}}, {'start'=>$2, 'end'=>$3})}
			if ($keepinfo==1){$regions{$1}[-1]{'info'}=$4}
			}
		}
	close $bed;
	$self->{_Loci}=\%regions;
	return;
	}


sub readGTF{
	my ($self, $filename, $keepinfo, $type)=@_;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	my %regions;

	open my $input, "<", $filename or die "Could not open bed file $filename!\n";

	while (my $line=<$input>){
		chomp $line;
		if ($line=~/^\w+\t\w+\t\w+/){
			my @col=split(/\t/, $line);
			unless (defined $type && $col[2] ne $type){
				my ($gene_id)=$col[8]=~/gene_id\s+"(\w+)";/;
				my ($transcript_id)=$col[8]=~/transcript_id\s+"(\w+)";/;
#				print STDERR "$col[0]:$col[3]-$col[4], $col[6], $gene_id\n";
				push( @{$regions{ $col[0] }}, {'type'=>$col[2], 'start'=>($col[3]-1), 'end'=>($col[4]-1), 'dir'=>$col[6], 'frame'=>$col[7], 'geneid'=>$gene_id, 'transcriptid'=>$transcript_id} );
				if ($keepinfo==1){ $regions{ $col[0] }[-1]{'info'}=$col[8] }
				}
			}
		}
	close $input;
	$self->{_Loci}=\%regions;
	return;
	}


sub readBins{
	my ($self, @filenames)=@_;
	unless(@_>=2){die "Wrong number of arguments!\n"}

	for my $filename (@filenames){
		open my $input, "<", $filename or die "Could not open bin file $filename!\n";

		while (my $line=<$input>){
			chomp $line;
			if ($line=~/^"?(\w+)"?\s+(\d+)/){
				if (exists $self->{_Bins}{$1}){
					if ($self->{_Bins}{$1} ne $2){ warn "GeneID $1 has already been assigned to a recombination bin! ", $self->{_Bins}{$1}, " vs. $2\n" }
					}
				else { $self->{_Bins}{$1}=$2 }
				}
			}
		close $input;
		}
	return;
	}


sub checkRegions{
	my ($self, $minlength)=@_;
	unless(@_==2){die "Wrong number of arguments!\n"}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $self->allKeys() ){
		my ($prevstart, $prevend);
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			if ($locuslength<$minlength){ warn "Warning: locus $chrom:$start-$end has length below minimum locus length! $locuslength vs. $minlength\n" }
			if (defined $prevend && $start<=$prevend){ warn "Warning: overlap detected between $chrom:$prevstart-$prevend and $chrom:$start-$end!\n" }
			$prevstart=$start;
			$prevend=$end unless (defined $prevend && $prevend>$end);
			}
		}
	return;
	}

sub cleanRegions{
	my ($self, $minlength, $sort, $remove_overlap, $ref_bins)=@_;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	my %cleaned;
	my %bins;
	%bins=map { $_=>1 } @$ref_bins if defined $ref_bins;

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $self->allKeys() ){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} }
		else { @loci=@{$self->{_Loci}{$chrom}} }
		my $prevlocus;
		my $discard=0;
		LOCUS: for my $locus (@loci){
			$locus->{'bin'} //= $self->{_Bins}{ $locus->{'geneid'} };	# assign bin hash element if not already defined.
			if (defined $ref_bins && $ref_bins->[0] ne 'all'){
				next LOCUS unless ( defined $locus->{'bin'} && $bins{ $locus->{'bin'} } );
				}
			my $accepted=1;
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			if ($locuslength<$minlength){
				warn "Warning: locus $chrom:$start-$end has length below minimum locus length! $locuslength vs. $minlength\n";
				next LOCUS;
				}
			unless (defined $prevlocus){ $prevlocus=$locus; next LOCUS }
			if ( $start<=$prevlocus->{'end'} ){
				my $infostring="$chrom:$prevlocus->{'start'}-$prevlocus->{'end'} (gene_id: $prevlocus->{'geneid'}, transcript_id: $prevlocus->{'transcriptid'}) and $chrom:$start-$end (gene_id: $locus->{'geneid'}, transcript_id: $locus->{'transcriptid'})";
#				print "Warning: overlap detected between $infostring!\n";
				if ( $prevlocus->{'geneid'} eq $locus->{'geneid'} ){
					unless ( $prevlocus->{'dir'} eq $locus->{'dir'} ){
						print STDERR "$infostring\n";
						die "Overlapping intervals have the same gene id but are on different strands!\n";
						}
					if ( $start==$prevlocus->{'start'} && $end==$prevlocus->{'end'} ){ next LOCUS }
					my $prevshift=$prevlocus->{'frame'};
					if ($prevlocus->{'dir'} eq '-'){ $prevshift=( $prevlocus->{'end'} - $prevlocus->{'start'} + 1 - $prevlocus->{'frame'} ) % 3 }
					my $shift=$locus->{'frame'};
					if ($locus->{'dir'} eq '-'){ $shift=( $locuslength - $locus->{'frame'} ) % 3 }
					if ( ($prevlocus->{'start'} + $prevshift - $start - $shift) % 3 ){	# check if both overlapping intervals are in the same reading frame.
						print STDERR "$infostring\n";
						warn "Overlapping intervals have the same gene id but are not in the same reading frame! Discarding both intervals.\n";
						print STDERR "$chrom:$prevlocus->{'start'}-$prevlocus->{'end'}, direction: $prevlocus->{'dir'}, frame: $prevlocus->{'frame'}, shift: $prevshift, codon start: ", $prevlocus->{'start'} + $prevshift, "\n"; 
						print STDERR "$chrom:$start-$end, direction: $locus->{'dir'}, frame: $locus->{'frame'}, shift: $shift, codon start: ", $start + $shift, "\n";
						$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
						$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};	# fuse intervals to detect further overlaps with either interval.
						$discard=1;
						next LOCUS;
						}
					else {
						print STDERR "$infostring\n";
						$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
						$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};
						warn "Overlapping intervals have same gene id, fusing intervals to $chrom:$prevlocus->{'start'}-$prevlocus->{'end'}!\n";
						next LOCUS;
						}
					}
				elsif ($remove_overlap){
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids, removing both intervals!\n";
					print STDERR "Removing $prevlocus->{'geneid'}\t$prevlocus->{'transcriptid'}\n";
					print STDERR "Removing $locus->{'geneid'}\t$locus->{'transcriptid'}\n";
					$prevlocus={ %$prevlocus };	# make a copy of the prevlocus hash to avoid changing the original locus hash.
					$prevlocus->{'end'}=$end if $end > $prevlocus->{'end'};	# fuse intervals to detect further overlaps with either interval.
					$discard=1;
					next LOCUS;
					}
				elsif ($prevlocus->{'dir'} ne $locus->{'dir'}){
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids and different directions, keeping both intervals!\n";
					}
				else {
					print STDERR "$infostring\n";
					warn "Overlapping regions have different gene ids and same directions, keeping longer interval!\n";
					if ( $locuslength < ($prevlocus->{'end'}-$prevlocus->{'start'}+1) ){ print "Removing $chrom:$start-$end!\n"; next LOCUS }
					else { print STDERR "Removing $chrom:$prevlocus->{'start'}-$prevlocus->{'end'}!\n"; $accepted=0 }
					}
				}
			$accepted=0 if $discard;
			push( @{$cleaned{$chrom}}, { %$prevlocus } ) if $accepted;
			$prevlocus=$locus;
			$discard=0;
			}
		push( @{$cleaned{$chrom}}, $prevlocus ) if (defined $prevlocus && !$discard);
		}
	my $cleaned_genes=Regions->new(\%cleaned);
	return($cleaned_genes);
	}

sub chopRegions{
	my ($self, $rmstart, $rmend, $minlength)=@_;
	unless(@_==4){die "Wrong number of arguments!\n"}
	my $minsize=$rmstart+$rmend+$minlength;

	CHROM: for my $chrom (keys %{$self->{_Loci}}){
		my @chopped;
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}} ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			unless($locuslength){ warn "$chrom:$start-$end has smaller end than start coordinate!\n" }
			unless($locuslength>$minsize){ warn "$chrom:$start-$end is below minimum size after chopping, locus omitted!\n"; next LOCUS }
			$locus->{'start'}=$start+$rmstart;
			$locus->{'end'}=$end-$rmend;
			push @chopped, $locus;
			}
		$self->{_Loci}{$chrom}=\@chopped;
		}

	return;
	}


sub concatRegions{
	my ($self)=@_;
	unless(@_==1){die "Wrong number of arguments!\n"}

	my %concat;
	CHROM: for my $chrom (keys %{$self->{_Loci}}){
#		my $totsize=0;
		my ($ref_locus, $newstart, $prevend);
		LOCUS: for my $locus ( sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } ){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			unless (defined $newstart){ $newstart=$start }
#			$totsize=$end-$newstart+1;
			if ( defined $prevend && $start>$prevend+1 ){
				push @{$concat{$chrom}}, { %$ref_locus };
				$concat{$chrom}[-1]{'start'}=$newstart;
				$concat{$chrom}[-1]{'end'}=$prevend;
				$newstart=$start;
#				$totsize=0;
				}
			$prevend=$end unless (defined $prevend && $prevend>$end);
			$ref_locus=$locus;
			}
		push @{$concat{$chrom}}, { %$ref_locus };
		$concat{$chrom}[-1]{'start'}=$newstart;
		$concat{$chrom}[-1]{'end'}=$prevend;
		}
	my $newloci=Regions->new(\%concat);
	return($newloci);
	}


sub windowingLoci{
	my ($self, $minlength, $windowsize, $stepsize, $offset, @chromosomes)=@_;
	my $verbose=0;
	unless(@_>=3){die "Wrong number of arguments!\n"}
	$stepsize=$windowsize unless (defined $stepsize);
	$offset=0 unless (defined $offset);
	unless (@chromosomes && $chromosomes[0] ne 'all'){ @chromosomes=keys %{$self->{_Loci}} }

	my %windows;
	CHROM: for my $chrom (@chromosomes){
		LOCUS: for my $locus (sort { $a->{'start'}<=>$b->{'start'} } @{$self->{_Loci}{$chrom}}){
			my $start=$locus->{'start'};
			my $end=$locus->{'end'};
			my $locuslength=$end-$start+1;
			next LOCUS unless ($locuslength>=$minlength);
			my $windowstart=$start+$offset;
			print STDERR "$chrom, $start-$end\n" if $verbose;
			WINDOW: while ($windowstart<$end){
				my $windowend=$windowstart+$windowsize-1;
				if ($windowend>$end){
					$windowend=$end;
					last WINDOW if ($end-$windowstart+1<$minlength);	# Don't add incomplete windows if smaller than minlength.
					}
				push @{$windows{$chrom}}, {'start'=>$windowstart, 'end'=>$windowend};
#				print STDERR "$windowstart, $windowend\n";
				last WINDOW if $windowend==$end;
				$windowstart+=$stepsize;
				}
			}
		}
	my $newloci=Regions->new(\%windows);
	return($newloci);
	}


sub subsetLocibyID{
	my ($self, $ref_ids)=@_;
	unless(@_==2){die "Wrong number of arguments!\n"}
	my %subset;
	ID: for my $id (@$ref_ids){
		unless (defined $self->{_Loci}{$id} ){ warn "Could not find $id in regions object!\n"; next ID }
		push @{ $subset{$id} }, { %$_ } foreach (@{ $self->{_Loci}{$id} });
		}
	my $newloci=Regions->new(\%subset);
	return($newloci);
	}

sub subsetLoci{	# Attention: Uses hash references. Any change to the subsetted loci will also affect the full set of loci!
	my ($self, $maxsize, @chromosomes)=@_;
	unless(@_>=2){die "Wrong number of arguments!\n"}
	unless (@chromosomes && $chromosomes[0] ne 'all'){ @chromosomes=sort {Misc::expand($a) cmp Misc::expand($b)} keys %{ $self->{_Loci} } }
	my %subset;
	my $prevend=0;
	my $rt=0;
	CHROM: for my $chrom (@chromosomes){
		next CHROM unless (defined $self->{_Loci}{$chrom} );
		unless ( defined $self->{_LastIndex}{$chrom} ){
			@{$self->{_Loci}{$chrom} }=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} };	# make sure to have loci sorted by scaffold position.
			$self->{_LastIndex}{$chrom}=0;
			}
		my $endindex=$self->getSize($chrom)-1;
		next CHROM if $self->{_LastIndex}{$chrom}>$endindex;
		LOCUS: for my $locus ( $self->{_LastIndex}{$chrom}..$endindex){
#			if ($self->getExcluded($chrom, $locus)){ next LOCUS }
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			if ($start<$prevend){ $rt+=$end-$prevend }
			else { $rt+=$end-$start+1 }
			if ($rt<=$maxsize){
				push @{$subset{$chrom}}, $self->{_Loci}{$chrom}[$locus];	# makes sure to take over all the information in the hash.
#				$self->setExcluded($chrom, $locus, 1);
				}
			else {$self->{_LastIndex}{$chrom}=$locus; last CHROM}
			$prevend=$end;
			}
		$self->{_LastIndex}{$chrom}=$self->getSize($chrom);
		last CHROM if exists $subset{$chrom};	# only process one scaffold at a time to prevent disruption of scaffold order on chromosome (scaffolds are stored as hash elements)!
		}
	my $newloci=Regions->new(\%subset);
	return($newloci);
	}


sub excludeRegions{
	my ($self, $regions, $min_dist)=@_;

	my %reduced;	

	SCAFF: for my $chrom (sort keys %{$self->{_Loci}}){
		LOCUS: for my $locus ( @{ $self->{_Loci}{$chrom} } ){
			if ( defined $regions->getValue($chrom) ){
				REG: for my $reg ( sort { $a->{'start'}<=>$b->{'start'} } $regions->allValues($chrom) ){
					if ( ($reg->{'end'}+$min_dist)>=$locus->{'start'} && ($reg->{'start'}-$min_dist)<=$locus->{'end'} ){
						warn "Overlap detected between $chrom $reg->{'start'}-$reg->{'end'} and $locus->{'start'}-$locus->{'end'}\n";
						next LOCUS;
						}
					elsif ( $reg->{'start'}>$locus->{'end'} ){ last REG }
					}
				}
			push @{ $reduced{$chrom} }, $locus;
			}
		}

	my $newloci=Regions->new(\%reduced);
	return($newloci);
	}


sub testProximity{
	my ($self, $candidates, $min_dist)=@_;

	my @accepted;

	CAND: for my $cand ($candidates->allValues()){
		my $id=$cand->{'id'};
		my $start=$cand->{'start'};
		my $end=$cand->{'end'};

		if (defined $self->{_Loci}{$id}){
			REG: for my $reg (@{$self->{_Loci}{$id}}){
				if (($reg->{'end'}+$min_dist)>=$start && ($reg->{'start'}-$min_dist)<=$end){
					warn "Overlap detected between $id: $reg->{'start'}-$reg->{'end'} and $start-$end, locus rejected!\n";
					push @accepted, 0;
					next CAND;
					}
				}
			}
		print STDERR "No overlap detected for locus $id: $start-$end, locus accepted.\n";
		
		$self->addLocus('id'=>$id, 'start'=>$start, 'end'=>$end);	# uses now named arguments.
		push @accepted, 1;
		}

	return;
	}


sub samplefromFasta{
	my ($self, $fasta, $genes, $peaks, $chrommap, $min_length, $nloci, $locus_length)=@_;

	my $found=0;
	my @accepted;

	LOCUS: while ($found<$nloci){
		my $candidates=$fasta->sampleLoci(1, $locus_length, 0.1);
		my $trcandidates=$candidates->translateScafftoChrom($chrommap, 1);
		my $id=$trcandidates->getValue(0, 'id');
		if ($id eq "NA"){next LOCUS}
		my $start=$trcandidates->getValue(0, 'start');
		my $end=$trcandidates->getValue(0, 'end');

		@accepted=();
		@accepted=$trcandidates->testbyRegions($self, 100000);
		@accepted=$trcandidates->testbyRegions($genes, 1000) if $accepted[0];
		@accepted=$candidates->testbyRegions($peaks, 10000) if $accepted[0];
#		@accepted=$self->testProximity($trcandidates, 100000) if $accepted[0];

		if ($accepted[0]){
			$self->addLocus('id'=>$id, 'start'=>$start, 'end'=>$end);	# uses now named arguments.
			$found++;
			print STDERR "$id: $start-$end accepted!\n";
			print STDERR "Found a total of $found of $nloci loci!\n";
			}
		else {warn "$id: $start-$end excluded!\n"}
		}

#	$loci->printLoci("-");
	return;
	}


sub samplefromFasta_mod{
	my ($self, $fasta, $foundloci, $genes, $peaks, $chrommap, $nloci, $locus_length, $amongdist, $excludeX, $excludeZ, $excludeMT, $allsitesfolder, $allsiteslist, $vcffolder, $vcflist, $minind, $minvq, $minmq, $mincov, $testcov, $minpropcov, $testbiallelic, $maxnotbiallelic, $maxiter)=@_;

	my ($allsitesfiles, $allsitessamples)=Misc::openvcflist($allsitesfolder, $allsiteslist);
	my ($vcffiles, $vcfsamples)=Misc::openvcflist($vcffolder, $vcflist);
	my $found=0;
	$maxiter=100*$nloci unless (defined $maxiter);
	my @accepted;

	LOCUS: while ($found<$nloci && $maxiter--){
		my $candidates=$fasta->sampleLoci(1, $locus_length, 0.1);
		my $trcandidates;
		if ($chrommap ne "none"){$trcandidates=$candidates->translateScafftoChrom($chrommap, 1)} else {$trcandidates=$candidates}
		my $id=$trcandidates->getValue(0, 'id');
		if ($id eq "NA" || $id=~/random/ || $id=~/Un/){ next LOCUS }
		if ( ($excludeX && $id=~/X/) || ($excludeZ && $id=~/Z/) || ($excludeMT && $id=~/MT/) ){ next LOCUS }
		my $start=$trcandidates->getValue(0, 'start');
		my $end=$trcandidates->getValue(0, 'end');

		@accepted=();
		@accepted=$trcandidates->testbyRegions($self, $amongdist);
		if ($accepted[0] && $foundloci ne "none"){ @accepted=$trcandidates->testbyRegions($foundloci, $amongdist) }
		if ($accepted[0] && $genes ne "none"){ @accepted=$candidates->testbyRegions($genes, 10000) }
		if ($accepted[0] && $peaks ne "none"){ @accepted=$candidates->testbyRegions($peaks, 10000) }
		if ($accepted[0] && $testcov){ @accepted=$candidates->testCoverage($allsitesfiles, $allsitessamples, $minind, $minvq, $minmq, $mincov, $minpropcov) }
		if ($accepted[0] && $testbiallelic){ @accepted=$candidates->testBiallelic($vcffiles, $vcfsamples, $minvq, $mincov, $maxnotbiallelic) }

		if ($accepted[0]){
			my $refseq=$fasta->randomAccess( $candidates->getValue(0, 'id'), $candidates->getValue(0, 'start'), $candidates->getValue(0, 'end'), 1 );
			$self->addLocus('id'=>$id, 'start'=>$start, 'end'=>$end, 'refseq'=>$$refseq);	# uses now named arguments.
			++$found;
			print STDERR "$id: $start-$end accepted!\n";
			print STDERR "Found a total of $found of $nloci loci!\n";
			}
		else {warn "$id: $start-$end excluded!\n"}
		}

	unless ($found==$nloci){ print STDERR "Could only find $found out of the requested $nloci loci within the maximum number of search iterations!\n" }
	return;
	}


sub excludeCoveragebyPosition{
	my ($self, $ref_missing, $mincov, $minprop, $minind, $haplotize)=@_;
	#$type('individual', 'position');

	my %locsummary; # %locsummary{chrom}[locus]
#	my %indsummary; # %indsummary{chrom}[locus][pop][ind]=index
	my %possummary; # %possummary{chrom}{pos}[pop][ind]=(0/1)
	my %reduced;
	my $npops=$ref_missing->getNPop();

	for my $chrom (sort keys %{$self->{_Loci}}){
		for my $locus (0..$#{$self->{_Loci}{$chrom}}){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my %loccov;
			my $validpos=0;
			my $offset=0;

			POSITION: for my $pos ($start..$end){
				my $valid=1;
				my $allind=0;
				POP: for my $pop (0..($npops-1)){
					my $nind=$ref_missing->getNInd($pop);
					my @indcov;
					my $validind=0;
					IND: for my $ind (0..$nind-1){
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						if (defined $coverage && $coverage>=$mincov){
							push @indcov, $coverage;
							$validind++;
							}
						else {
							$ref_missing->setValue2($chrom, $locus, $allind, $offset, 0);
							push @indcov, 0;
							}
						$allind++;
						}
					my @rawindices=sort {$a<=>$b} ( (sort {$indcov[$b]<=>$indcov[$a]} (0..$#indcov))[ 0..$minind->[$pop]-1 ] );
					my @indices;
					if ($haplotize){@indices=map {2*$_+int(rand(2))} @rawindices}
					else {@indices=map {(2*$_, 2*$_+1)} @rawindices}
#					my @temp=(sort {$indcov[$b]<=>$indcov[$a]} (0..$#indcov))[0..$minind->[$pop]-1];
#					push @indices, ($_, $_+1) foreach (sort {$a<=>$b} @temp);
					$loccov{$pos}[$pop]=\@indices;
					if ($validind<$minind->[$pop]){$valid=0}
					}
				$validpos++ if $valid;
				$offset++;
				}

			my $propcov=$validpos/$locuslength;
			if ($propcov>=$minprop){
				push @{$reduced{$chrom}}, $self->{_Loci}{$chrom}[$locus];
				push @{$locsummary{$chrom}}, 1;
#				my @ref_pop;
#				for my $pop (0..($ref_missing->getNPop()-1)){push @ref_pop, [0..$ref_missing->getNInd($pop)-1]};
#				push @{$indsummary{$chrom}}, \@ref_pop;
				for my $pos ($start..$end){$possummary{$chrom}{$pos}=$loccov{$pos}}
				}
			else {push @{$locsummary{$chrom}}, 0}
			print STDERR "$chrom: $start-$end\t$propcov\n";
			}
		}

	$self->{_LocSummary}=\%locsummary;
	my $newloci=Regions->new(\%reduced);
#	$newloci->{_IndSummary}=\%indsummary;
	$newloci->{_PosSummary}=\%possummary;
	return($newloci);
	}


sub excludeCoveragebyLocus{
	my ($self, $ref_missing, $mincov, $minprop, $minind, $haplotize, $subind)=@_;

	my %indsummary; # %indsummary{chrom}[locus][pop][ind]=index
	my %reduced;
	my $npops=$ref_missing->getNPop();

	for my $chrom (sort keys %{$self->{_Loci}}){
		for my $locus (0..$#{$self->{_Loci}{$chrom}}){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $valid=1;
			my @ref_pop;
			my $allind=0;

			for my $pop (0..($npops-1)){
				my $nind=$ref_missing->getNInd($pop);
				my $validind=0;
				my @propcov;

				for my $ind (0..$nind-1){
					my $covered=0;
					for my $offset (0..$end-$start){
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						if (defined $coverage && $coverage>=$mincov){ ++$covered }
						else {
							$ref_missing->setValue($chrom, $start, $allind, $offset, 0);
							}
						}
					push @propcov, ($covered/$locuslength);
					if ($propcov[-1]>=$minprop){++$validind}
					++$allind;
#					print STDERR "$chrom: $start-$end\t$ind: $covered - $propcov[-1]\n";
					}
				print STDERR "Population $pop - $chrom: $start-$end: ", $validind, "\n";
				my @rawindices;
				if (defined $subind){ @rawindices=sort {$a<=>$b} ( (sort {$propcov[$b]<=>$propcov[$a]} (0..$#propcov))[ 0..$subind->[$pop]-1 ] ) }
				else { @rawindices=(0..$#propcov) }
				my @indices;
				if ($haplotize){ @indices=map { 2*$_+int(rand(2)) } @rawindices }
				else { @indices=map { (2*$_, 2*$_+1) } @rawindices }
				push @ref_pop, \@indices;
				$valid=0 unless ($validind>=$minind->[$pop]);
				}
			$self->{_IndSummary}{$chrom}[$locus]=\@ref_pop;

			if ($valid){
				$self->{_LocSummary}{$chrom}[$locus]=1;	# is this needed? Basically the same as {'excluded'}.
				$self->setExcluded($chrom, $locus, 0);
				push @{$reduced{$chrom}}, $self->{_Loci}{$chrom}[$locus];
				push @{$indsummary{$chrom}}, \@ref_pop;
				}
			else {
				$self->{_LocSummary}{$chrom}[$locus]=0;
				$self->setExcluded($chrom, $locus, 1);
				warn "Locus $chrom: $start-$end excluded!\n";
				}
			}
		}
	my $newloci=Regions->new(\%reduced);
	$newloci->{_IndSummary}=\%indsummary;
	return($newloci);
	}

sub excludeCoveragebyLocusBlocks{
	my ($self, $ref_missing, $mincov, $minprop, $minind, $haplotize, $subind)=@_;

#	$subind=$minind unless (defined $subind);
	my %indsummary; # %indsummary{chrom}[locus][pop][ind]=index
	my %reduced;
	my $npops=$ref_missing->getNPop();

	for my $chrom (sort keys %{$self->{_Loci}}){
		for my $locus (0..$#{$self->{_Loci}{$chrom}}){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $valid=1;
			my @ref_pop;
			my $allind=0;

			for my $pop (0..($npops-1)){
				my $nind=$ref_missing->getNInd($pop);
				my $validind=0;
				my @propcov;

				for my $ind (0..$nind-1){
					my $covered=0;
					for my $pos ($start..$end){
						my $coverage=$ref_missing->getValue3($chrom, $pos, $allind);
						if (defined $coverage && $coverage>=$mincov){ ++$covered }
						else {
							$ref_missing->setValue3($chrom, $pos, $allind, 0);
							}
						}
					push @propcov, ($covered/$locuslength);
					if ($propcov[-1]>=$minprop){++$validind}
					++$allind;
#					print STDERR "$chrom: $start-$end\t$ind: $covered - $propcov[-1]\n";
					}
				print STDERR "Population $pop - $chrom: $start-$end: ", $validind, "\n";
				my @rawindices;
				if (defined $subind){ @rawindices=sort {$a<=>$b} ( (sort {$propcov[$b]<=>$propcov[$a]} (0..$#propcov))[ 0..$subind->[$pop]-1 ] ) }
				else { @rawindices=(0..$#propcov) }
				my @indices;
				if ($haplotize){ @indices=map { 2*$_+int(rand(2)) } @rawindices }
				else { @indices=map { (2*$_, 2*$_+1) } @rawindices }
				push @ref_pop, \@indices;
				$valid=0 unless ($validind>=$minind->[$pop]);
				}
			$self->{_IndSummary}{$chrom}[$locus]=\@ref_pop;

			if ($valid){
				$self->{_LocSummary}{$chrom}[$locus]=1;	# is this needed? Basically the same as {'excluded'}.
				$self->setExcluded($chrom, $locus, 0);
				push @{$reduced{$chrom}}, $self->{_Loci}{$chrom}[$locus];
				push @{$indsummary{$chrom}}, \@ref_pop;
				}
			else {
				$self->{_LocSummary}{$chrom}[$locus]=0;
				$self->setExcluded($chrom, $locus, 1);
				warn "Locus $chrom: $start-$end excluded!\n";
				}
			}
		}
	my $newloci=Regions->new(\%reduced);
	$newloci->{_IndSummary}=\%indsummary;
	return($newloci);
	}


sub Haplotyze{
	my ($self, $mode)=@_;

	if ($mode eq 'position'){
		CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_PosSummary}}){
			POSITION: while (my ($pos, $ref_pos)=each %{$ref_chrom}){
				POP: for my $pop (0..$#{$ref_pos}){
					my @indices;
					IND: for (my $ind=0;$ind<$#{$ref_pos->[$pop]};$ind=$ind+2){
						push @indices, $ref_pos->[$pop][$ind]+int(rand(2));
						}
					$ref_pos->[$pop]=\@indices;
					}
				}
			}
		}
	elsif ($mode eq 'locus'){
		CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_IndSummary}}){ #[locus][pop][ind]
			LOCUS: for my $ref_locus (@{$ref_chrom}){
				POP: for my $pop (0..$#{$ref_locus}){
					my @indices;
					IND: for (my $ind=0;$ind<$#{$ref_locus->[$pop]};$ind=$ind+2){
						push @indices, $ref_locus->[$pop][$ind]+int(rand(2));
						}
					$ref_locus->[$pop]=\@indices;
					}
				}
			}
		}
	else {die "Mode not recognized!\n"}
	return;
	}


sub translateScafftoChrom{
	my ($self, $chrommap)=@_; # bed file format with chromosome locations

	my %translated;

	keys %{$self->{_Loci}};
	SCAFF: while (my ($scaff, $ref_scaff)=each %{$self->{_Loci}}){
		CHROM: for my $chrom (keys %{$chrommap->{_Chrom}}){
			if (defined $chrommap->getValue($chrom, $scaff)){
				my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
				if ($dir eq '-'){
					my $offset=$chrommap->getValue($chrom, $scaff, 'end');
					LOCUS: for my $locus (@{$ref_scaff}){
						push @{$translated{$chrom}}, {%{$locus}};
						$translated{$chrom}[-1]{'start'}=$offset-$locus->{'end'};
						$translated{$chrom}[-1]{'end'}=$offset-$locus->{'start'};
						print STDERR "Scaffold $scaff: $locus->{'start'}-$locus->{'end'} translated to chromosome $chrom: $translated{$chrom}[-1]{'start'}-$translated{$chrom}[-1]{'end'}.\n";
						}
					}
				else {
					my $offset=$chrommap->getValue($chrom, $scaff, 'start');
					LOCUS: for my $locus (@{$ref_scaff}){
						push @{$translated{$chrom}}, {%{$locus}};
						$translated{$chrom}[-1]{'start'}+=$offset;
						$translated{$chrom}[-1]{'end'}+=$offset;
						print STDERR "Scaffold $scaff: $locus->{'start'}-$locus->{'end'} translated to chromosome $chrom: $translated{$chrom}[-1]{'start'}-$translated{$chrom}[-1]{'end'}.\n";
						}
					}
				if ($dir ne '+' && $dir ne '-') {warn "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n"}
				next SCAFF;
				}
			else {next CHROM}
			}
		warn "Warning: scaffold $scaff has not been located on a chromosome! Loci on this scaffold are excluded in the output!\n";
		}

	my $tr_loci=Regions->new(\%translated);
	return($tr_loci);
	}


sub translateChromtoScaff{	# only works if whole region is located on a single scaffold!
	my ($self, $chrommap)=@_; # bed file format with chromosome locations

	my %translated;

	keys %{$self->{_Loci}};
	CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_Loci}}){
		LOCUS: for my $locus (@{$ref_chrom}){
			SCAFF: for my $scaff (keys %{$chrommap->{_Chrom}{$chrom}}){
				my $scaffstart=$chrommap->getValue($chrom, $scaff, 'start');
				my $scaffend=$chrommap->getValue($chrom, $scaff, 'end');
				if ($locus->{'start'}>=$scaffstart && $locus->{'end'}<=$scaffend){
					push @{$translated{$scaff}}, {%{$locus}};
					my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
					if ($dir eq '-'){
						$translated{$scaff}[-1]{'start'}=$scaffend-$locus->{'end'};
						$translated{$scaff}[-1]{'end'}=$scaffend-$locus->{'start'};
						}
					else {
						$translated{$scaff}[-1]{'start'}-=$scaffstart;
						$translated{$scaff}[-1]{'end'}-=$scaffstart;
						}
					if ($dir ne '+' && $dir ne '-') {warn "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n"}
					print STDERR "Chromosome $chrom: $locus->{'start'}-$locus->{'end'} translated to scaffold $scaff: $translated{$scaff}[-1]{'start'}-$translated{$scaff}[-1]{'end'}.\n";
					next LOCUS;
					}
				}
			warn "Warning: Chromosome $chrom: $locus->{'start'}-$locus->{'end'} has not been located on a scaffold! The locus is not included in the output!\n";
			}
		}

	my $tr_loci=Regions->new(\%translated);
	return($tr_loci);
	}


sub translateChromtoScaffmod{
	my ($self, $chrommap)=@_; # bed file format with chromosome locations

	my %translated;

	keys %{$self->{_Loci}};
	CHROM: while (my ($chrom, $ref_chrom)=each %{$self->{_Loci}}){
		my @scafflist=sort {$chrommap->{_Chrom}{$chrom}{$a}{'start'}<=>$chrommap->{_Chrom}{$chrom}{$b}{'start'} } ( $chrommap->allKeys($chrom) );
		LOCUS: for my $locus (@{$ref_chrom}){
			my @foundscaffs;
			SCAFF: for my $scaff (@scafflist){
				my $scaffstart=$chrommap->getValue($chrom, $scaff, 'start');
				my $scaffend=$chrommap->getValue($chrom, $scaff, 'end');
				if ($scaffend<$locus->{'start'}){ next SCAFF }
				elsif ($scaffstart>$locus->{'end'}){ last SCAFF }
				else {
					my $frontoverhang=$locus->{'start'}-$scaffstart;
					my $backoverhang=$scaffend-$locus->{'end'};
					my $dir=$chrommap->getValue($chrom, $scaff, 'dir');
					push @{$translated{$scaff}}, {%{$locus}};
					my ($trstart, $trend);
					if ($dir eq '-'){
						$trstart=$backoverhang>0 ? $backoverhang : 0;
						$trend=$frontoverhang>0 ? $scaffend-$locus->{'start'} : $scaffend-$scaffstart;
						}
					else {
						$trstart=$frontoverhang>0 ? $frontoverhang : 0;
						$trend=$backoverhang>0 ? $locus->{'end'}-$scaffstart : $scaffend-$scaffstart;
						}
					if ($dir ne '+' && $dir ne '-') {warn "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n"}
					$translated{$scaff}[-1]{'start'}=$trstart;
					$translated{$scaff}[-1]{'end'}=$trend;
					push @foundscaffs, [ $scaff, $trstart, $trend ];
					}
				}
			if (@foundscaffs){
				print STDERR "Chromosome $chrom: $locus->{'start'}-$locus->{'end'} has been translated to the following scaffolds:\n";
				print STDERR "$_->[0]:$_->[1]-$_->[2]\n" for (@foundscaffs);
				}
			else { warn "Warning: Chromosome $chrom: $locus->{'start'}-$locus->{'end'} has not been located on a scaffold! The locus is not included in the output!\n" }
			}
		}

	my $tr_loci=Regions->new(\%translated);
	return($tr_loci);
	}


sub calculatePi{	# number of pairwise differences within populations, calculated per site
	my ($self, $genotypes, $ref_missing, $mincov, $minprop, $minind)=@_;

	my $npops=$ref_missing->getNPop();

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		LOCUS: for my $locus (0..$#{$self->{_Loci}{$chrom}}){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $rt=0;

			POP: for my $pop (0..($npops-1)){
				my $nind=$ref_missing->getNInd($pop);
				my $validpos=0;
				my $pi=0;
				POSITION: for my $pos ($start..$end){
					my $offset=$pos-$start;
					my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
					my $validind=0;
					my @alleles=(0,0,0,0,0);
					IND: for my $index (0..2*$nind-1){
						my $allind=int($index/2)+$rt;
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						if (defined $coverage && $coverage>=$mincov){
							if ($snp){
								my $call=$genotypes->getValue($chrom, $pos, $pop, $index);
								if ($call eq 0 || $call eq 1 || $call eq 2 || $call eq 3){$alleles[$call]++; $validind++}
								else {$alleles[4]++}
								}
							else {$validind++}
							}
						}
					if ($validind>=2*$minind->[$pop]){
						$validpos++;
						if ($snp){
							my $total=0;
							my $denom=$validind*($validind-1);
							for my $call (0..3){
								if (defined $alleles[$call] && $alleles[$call]>0){$total+=($alleles[$call]*($alleles[$call]-1))/$denom}
								}
							$pi+=(1-$total);
#							print STDERR "$chrom: $pos - ", join(":", @alleles), "\t$total\n";
							}
						}
					}
				my $propcov=$validpos/$locuslength;
				if ($propcov>=$minprop){
					my $pi_persite=$pi/$validpos;
					$self->{_Pi}{$chrom}[$locus][$pop]=$pi_persite;
					$self->{_Covered}{$chrom}[$locus][$pop]=$validpos;
					print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\tpi: $pi_persite\n";
					}
				else {warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n"}
				$rt+=$nind;
				}
			}
		}
	return;
	}


sub calculateDxy{	# number of pairwise differences between populations, calculated per site
	my ($self, $genotypes, $ref_missing, $mincov, $minprop, $minind)=@_;

	my $npops=$ref_missing->getNPop();
	if ($npops>2){die("More than two populations in the dataset!")}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		LOCUS: for my $locus (0..$#{$self->{_Loci}{$chrom}}){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $validpos=0;
			my $valid;
			my $dxy=0;
			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
				my @nsamples;
				my @alleles=([0,0,0,0,0], [0,0,0,0,0]);
				my $rt=0;
				POP: for my $pop (0..($npops-1)){
					my $nind=$ref_missing->getNInd($pop);
					my $validind=0;
					IND: for my $index (0..2*$nind-1){
						my $allind=int($index/2)+$rt;
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
						if (defined $coverage && $coverage>=$mincov){
							if ($snp){
								my $call=$genotypes->getValue($chrom, $pos, $pop, $index);
								if ($call eq 0 || $call eq 1 || $call eq 2 || $call eq 3){$alleles[$pop][$call]++; $validind++}
								else {$alleles[$pop][4]++}
								}
							else {$validind++}
							}
						}
					if ($validind<2*$minind->[$pop]){next POSITION}
					push @nsamples, $validind;
					$rt+=$nind;
					}
				$validpos++;
				if ($snp){
					my $total=0;
					my $denom=$nsamples[0]*$nsamples[1];
					for my $call (0..3){
						if ($alleles[0][$call]>0 && $alleles[1][$call]>0){
							$total+=($alleles[0][$call]*$alleles[1][$call])/$denom;
							}
						}
					$dxy+=(1-$total);
#					print STDERR "$chrom: $pos - ", join(":", @{$alleles[0]}, @{$alleles[1]}), "\t$total\n";
					}
				}
			my $propcov=$validpos/$locuslength;
			if ($propcov>=$minprop){
				my $dxy_persite=$dxy/$validpos;
				$self->{_Dxy}{$chrom}[$locus]=$dxy_persite;
				print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\tdxy: $dxy_persite\n";
				}
			else {warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n"}
			}
		}
	return;
	}

sub calculateDxyPhasedOld{	# number of pairwise differences between populations, calculated per haplotype
	my ($self, $ref_fasta, $ref_genotypes, $ref_missing, $ref_groups, $ref_mingroupind, $mincov, $minprop, $excluderefN)=@_;

	my $npops=$ref_missing->getNPop();
	my $ref_popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	foreach (@$ref_popsizes){ push @popsum, $rt; $rt+=$_ }
	my @groupsizes=map { 2 * Misc::sum(@{ $ref_popsizes }[@$_]) } (@$ref_groups);
	foreach (0..$#groupsizes){ die "Zero individuals in group $_!\n" unless $groupsizes[$_] }

	my $poppair=0;
	for my $group1 ( 0..($#{$ref_groups}-1) ){
		for my $group2 ( ($group1+1)..$#{$ref_groups} ){
			CHROM: for my $chrom ( $self->allKeys() ){
				my ($end, $prevend);
				my $firststart;
				LOCUS: for my $locus ( $self->allKeys($chrom) ){
					my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
					$end=$self->{_Loci}{$chrom}[$locus]{'end'};
					$firststart=$start unless defined $firststart;
					if (defined $prevend && $start-$prevend>1){ $firststart=$start }
					my $locuslength=$end-$start+1;
					my $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1);
					if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $chrom:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }
					my @diffs=(0) x ($groupsizes[$group1]*$groupsizes[$group2]);
					my @valid=(0) x ($groupsizes[$group1]*$groupsizes[$group2]);
					my $validpos=0;

					POSITION: for my $pos ($start..$end){
						my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
						next POSITION if ($excluderefN && $refbase eq 'N');
						my $offset=$pos-$firststart;
						my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
						my @alleles=([ ('N') x $groupsizes[$group1] ], [ ('N') x $groupsizes[$group2] ]);

						for my $groupindex (0, 1){
							my $group=($group1, $group2)[$groupindex];
							my $alleleindex=0;
							my $validind=0;
							for my $subpop (@{ $ref_groups->[$group] }){
								IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
									my $allind=$ind+$popsum[$subpop];
									my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
									unless (defined $coverage && $coverage>=$mincov){ next IND }
									if ($snp){
										my $index=2*$ind;
										my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
										if (defined $allele1 && $allele1 ne '.' && $allele1<=3){ $alleles[$groupindex][$alleleindex]=$allele1 }
										my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
										if (defined $allele2 && $allele2 ne '.' && $allele2<=3 ){ $alleles[$groupindex][$alleleindex+1]=$allele2 }
										}
									else {
										$alleles[$groupindex][$alleleindex]="0";
										$alleles[$groupindex][$alleleindex+1]="0";
										}
									++$validind;
									} continue { $alleleindex+=2 }
								}
							if ($validind<$ref_mingroupind->[$group]){ next POSITION }
							}

						my $comp=0;
#						print STDERR "$poppair - group1: ", join(',', @{ $alleles[0] }), ", group2: ", join(',', @{ $alleles[1] }), "\n";
						for my $allele1 (@{ $alleles[0] }){
							for my $allele2 (@{ $alleles[1] }){
								if ($allele1 ne 'N' && $allele2 ne 'N'){
									++$valid[$comp];
									++$diffs[$comp] if ($allele1 ne $allele2);
									}
								++$comp;
								}
							}

						++$validpos;
						}

					my $propcov=$validpos/$locuslength;
					my $minvalidsites=int($minprop*$locuslength) ? int($minprop*$locuslength) : 1;
					$self->{_Loci}{$chrom}[$locus]{'covered'}[$poppair]=$validpos;
					my $validcomps=0;
					if ($validpos>=$minvalidsites){
						print STDERR "minvalidsites: $minvalidsites\ndiffs: ", join(',', @diffs), "\nvalid: ", join(',', @valid), "\n";
						for my $comp (0..$#valid){
							if ($valid[$comp]>=$minvalidsites){ $diffs[$comp] /= $valid[$comp]; ++$validcomps }
							else { $diffs[$comp]=undef }
							}
						$self->{_Loci}{$chrom}[$locus]{'comps'}[$poppair]=$validcomps;
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]=$validcomps ? ( Misc::meansd(\@diffs) )[0] : "NA";
						$self->{_Loci}{$chrom}[$locus]{'Gmin'}[$poppair]=$validcomps ? Misc::min(\@diffs) : "NA";
						print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\n";
						}
					else {
						$self->{_Loci}{$chrom}[$locus]{'comps'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'Gmin'}[$poppair]="NA";
						warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
						}
					} continue { $prevend=$end }
				}
			++$poppair;
			}
		}
	return;
	}

sub calculateDxyPhased{	# number of pairwise differences between populations, calculated per haplotype
	my ($self, $ref_fasta, $ref_genotypes, $ref_missing, $ref_groups, $ref_mingroupind, $mincov, $minprop, $excluderefN, $allvalid)=@_;

	my $npops=$ref_missing->getNPop();
	my $ref_popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	foreach (@$ref_popsizes){ push @popsum, $rt; $rt+=$_ }
	my @groupsizes=map { 2 * Misc::sum(@{ $ref_popsizes }[@$_]) } (@$ref_groups);
	foreach (0..$#groupsizes){ die "Zero individuals in group $_!\n" unless $groupsizes[$_] }
	my $npairs=(scalar(@$ref_groups)*(scalar(@$ref_groups)-1))/2;

	CHROM: for my $chrom ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			$end=$self->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $minvalidsites=int($minprop*$locuslength) ? int($minprop*$locuslength) : 1;
			my $ref_seqstring=defined $ref_fasta ? $ref_fasta->randomAccess($chrom, $start, $end, 1): \("A" x $locuslength);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $chrom:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }

			my (@diffs, @valid, @ndims);
			my @validpos=(0) x $npairs;
			for my $group1 ( 0..$#groupsizes-1){
				for my $group2 ( ($group1+1)..$#groupsizes){
					push @diffs, [ (0) x ($groupsizes[$group1]*$groupsizes[$group2]) ];
					push @valid, [ (0) x ($groupsizes[$group1]*$groupsizes[$group2]) ];
					push @ndims, [ $groupsizes[$group1], $groupsizes[$group2] ];
					}
				}

			POSITION: for my $pos ($start..$end){
				my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				next POSITION if ($excluderefN && $refbase eq 'N');
				my $offset=$pos-$firststart;
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				my @alleles=map { [ ('N') x $_ ] } @groupsizes;
				my @groupvalid=(0) x scalar(@groupsizes);

				GROUP: for my $group (0..$#alleles){
					my $alleleindex=0;
					my $validind=0;
					for my $subpop (@{ $ref_groups->[$group] }){
						IND: for my $ind (0..$ref_popsizes->[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($chrom, $pos, $subpop, $index);
								if (defined $allele1 && $allele1 ne '.' && $allele1<=3){ $alleles[$group][$alleleindex]=$allele1 }
								my $allele2=$ref_genotypes->getValue($chrom, $pos, $subpop, ++$index);
								if (defined $allele2 && $allele2 ne '.' && $allele2<=3 ){ $alleles[$group][$alleleindex+1]=$allele2 }
								}
							else {
								$alleles[$group][$alleleindex]="0";
								$alleles[$group][$alleleindex+1]="0";
								}
							++$validind;
							} continue { $alleleindex+=2 }
						}
					if ($validind<$ref_mingroupind->[$group]){ $allvalid ? next POSITION : next GROUP }
					$groupvalid[$group]=1;
					}

				my $poppair=0;
				GROUP1: for my $group1 ( 0..$#groupsizes-1){
					GROUP2: for my $group2 ( ($group1+1)..$#groupsizes){
						unless ($groupvalid[$group1] && $groupvalid[$group2]){ next GROUP2 }
						my $comp=0;
						for my $allele1 (@{ $alleles[$group1] }){
							for my $allele2 (@{ $alleles[$group2] }){
								if ($allele1 ne 'N' && $allele2 ne 'N'){
									++$valid[$poppair][$comp];
									++$diffs[$poppair][$comp] if ($allele1 ne $allele2);
									}
								++$comp;
								}
							}
						++$validpos[$poppair]
						} continue { ++$poppair }
					}
				}

			$self->{_Loci}{$chrom}[$locus]{'covered'}=\@validpos;
			for my $poppair (0..$#validpos){
				my $propcov=$validpos[$poppair]/$locuslength;
				my $validcomps=0;
				if ($validpos[$poppair]>=$minvalidsites){
#					print STDERR "poppair: $poppair, minvalidsites: $minvalidsites\ndiffs: ", join(',', @{ $diffs[$poppair] }), "\nvalid: ", join(',', @{ $valid[$poppair] }), "\n";
					for my $comp (0..$#{ $valid[$poppair] } ){
						if ($valid[$poppair][$comp]>=$minvalidsites){ $diffs[$poppair][$comp] /= $valid[$poppair][$comp]; ++$validcomps }
						else { $diffs[$poppair][$comp]=undef }
						}
					$self->{_Loci}{$chrom}[$locus]{'comps'}[$poppair]=$validcomps;
					my ($dxy_mean, $dxy_sd)=$validcomps ? Misc::meansd($diffs[$poppair]) : ("NA", "NA");
					my ($min_dxy, $min_index)=$validcomps ? Misc::min($diffs[$poppair]) : ("NA", "NA");
					if ($min_index ne "NA"){
						my $size1=$ndims[$poppair][0];
						my $size2=$ndims[$poppair][1];
						my $index1=int($min_index/$size2);
						my $index2=$min_index % $size2;
#						print STDERR "Size1: $size1, size2: $size2, min_index: $min_index, index1: $index1, index2: $index2\n";
						my @dxy1_comps=( ($index1*$size2)..($index1*$size2+$size2-1) );
						my ($mean_dxy1, $sd_dxy1)=Misc::meansd([ @{ $diffs[$poppair] }[@dxy1_comps] ]);
#						print STDERR "Indices 1: ", join(',', @dxy1_comps), "\nElements 1: ", join(',', @{ $diffs[$poppair] }[@dxy1_comps]), "\n";
						my @dxy2_comps=map { $_*$size2+$index2 } (0..$size1-1);
						my ($mean_dxy2, $sd_dxy2)=Misc::meansd([ @{ $diffs[$poppair] }[@dxy2_comps] ]);
#						print STDERR "Indices 2: ", join(',', @dxy2_comps), "\nElements 2: ", join(',', @{ $diffs[$poppair] }[@dxy2_comps]), "\n";
						$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap1'}[$poppair]=$mean_dxy1;
						$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap2'}[$poppair]=$mean_dxy2;
						$self->{_Loci}{$chrom}[$locus]{'hap1/hap2'}[$poppair]=(defined $mean_dxy1 && $mean_dxy2) ? $mean_dxy1/$mean_dxy2 : "NA";
						}
					else {
						$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'hap1/hap2'}[$poppair]="NA";
						}
					$self->{_Loci}{$chrom}[$locus]{'dxy_mean'}[$poppair]=$dxy_mean;
					$self->{_Loci}{$chrom}[$locus]{'dxy_sd'}[$poppair]=$dxy_sd;
					$self->{_Loci}{$chrom}[$locus]{'dxy_min'}[$poppair]=$min_dxy;
					$self->{_Loci}{$chrom}[$locus]{'min_index'}[$poppair]=$min_index;
					$self->{_Loci}{$chrom}[$locus]{'Gmin'}[$poppair]=($dxy_mean ne "NA" && $dxy_mean>0) ? $min_dxy / $dxy_mean : "NA";
#					print STDERR "$chrom: $start-$end\tPopulation pair $poppair has $validpos[$poppair] of $locuslength ($propcov) sufficiently covered sites.\n";
					}
				else {
					$self->{_Loci}{$chrom}[$locus]{'comps'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'dxy_mean'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'dxy_sd'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'dxy_min'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'min_index'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap1'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'mean_dxy_hap2'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'hap1/hap2'}[$poppair]="NA";
					$self->{_Loci}{$chrom}[$locus]{'Gmin'}[$poppair]="NA";
					warn "$chrom: $start-$end\tPopulation pair $poppair has only $validpos[$poppair] of $locuslength ($propcov) sufficiently covered sites!\n";
					}
				}
			} continue { $prevend=$end }
		}
	return;
	}


sub thetasfromSFS{
	my ($self, $outfile, $sort, $folded)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "";
			GROUP: for my $metapop (0..$#{ $locus->{'SFS'} }){
				my $validsites=$locus->{'covered'}[$metapop];
				print $output "\t$validsites";
				unless ($validsites){
					print $output "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
					next GROUP;
					}
				my @totcounts;
				CAT: for my $cat ( @{ $locus->{'SFS'}[$metapop] } ){
					$totcounts[$_]+=$cat->[$_] foreach (0..$#{$cat});
					}
				my $nalleles=$folded ? (scalar @{ $locus->{'SFS'}[$metapop][0] }-1)*2 : scalar @{ $locus->{'SFS'}[$metapop][0] }-1;	# Wrong for folded SFS!!!
				my ($denom, $segsites, $thetaPi, $thetaH, $thetaL, $thetas1, $thetasrest)=(0,0,0,0,0,0,0);
				for my $i (1..$nalleles-1){
					my $count=$totcounts[$i];	# $locus->{'SFS'}[$metapop][0][$i];
					$denom+=(1/$i);
					if ($i==1){ $thetas1=$count }
					else { $thetasrest+=$count }
					$segsites+=$count;
					$thetaPi+=$i*($nalleles-$i)*$count;
					$thetaH+=($i**2)*$count;
					$thetaL+=$i*$count;
					}
				my $ndenom=($nalleles*($nalleles-1))/2;
				my $thetaW=$segsites/$denom;
				$thetaPi/=$ndenom;
				$thetaH/=$ndenom;
				my $TajD_denom=Misc::TajD_denominator($nalleles, $segsites);
				my $TajD=($TajD_denom>0) ? ($thetaPi-$thetaW)/$TajD_denom : "NA";
				my $FayWuH_denom=Misc::Hprime_denominator($nalleles, $segsites);
				my $FayWuH=($FayWuH_denom>0) ? ($thetaPi-$thetaH)/$FayWuH_denom : "NA";
				$thetas1/=$validsites;
				$thetasrest/=($denom-1)*$validsites;
				$thetaW/=$validsites;
				$thetaPi/=$validsites;
				$thetaH/=$validsites;
				$thetaL/=($nalleles-1)*$validsites;
				my $Tajd=$thetaPi-$thetaW;
				my $FuLiD=$thetasrest-$thetas1;
				my $FayWuh=$thetaPi-$thetaH;
				my $ZengE=$thetaPi-$thetaL;
				print $output "\t$segsites\t$thetaW\t$thetas1\t$thetasrest\t$thetaPi\t$thetaH\t$thetaL\t$Tajd\t$TajD\t$FuLiD\t$FayWuh\t$FayWuH\t$ZengE";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}

sub SSfrom2DSFS{	# folded SFS not yet properly implemented!
	my ($self, $outfile, $sort, $folded)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t$locus->{'info'}";
			PAIR: for my $poppair (0..$#{ $locus->{'SFS'} }){
				my $validsites=$locus->{'covered'}[$poppair];
				print $output "\t$validsites";
				unless ($validsites){
					print $output "\tNA\tNA\tNA\tNA\tNA\tNA";
					next PAIR;
					}
				my $ref_2dsfs=$locus->{'SFS'}[$poppair];	# $ref_2dsfs->[$i][$j];
				my $nalleles1=$folded ? (scalar @$ref_2dsfs-1)*2 : (scalar @$ref_2dsfs-1);
				my $nalleles2=$folded ? (scalar @{ $ref_2dsfs->[0]}-1)*2 : scalar @{ $ref_2dsfs->[0]}-1;
				my $totcount=0;
				my $denom1=($nalleles1*($nalleles1-1))/2;
				my $denom2=($nalleles2*($nalleles2-1))/2;
				my $denomtot=(($nalleles1+$nalleles2)*($nalleles1+$nalleles2-1))/2;

				my ($thetapi1, $thetapi2, $thetapitot, $dxy);
				for my $i (0..$nalleles1){
					for my $j (0..$nalleles2){
						my $count=$ref_2dsfs->[$i][$j];
						$totcount+=$count;
						$thetapi1+=$i*($nalleles1-$i)*$count;
						$thetapi2+=$j*($nalleles2-$j)*$count;
						$thetapitot+=($i+$j)*($nalleles1+$nalleles2-$i-$j)*$count;
						$dxy+=($i*($nalleles2-$j)+$j*($nalleles1-$i))*$count;
						}
					}
				unless ($totcount==$validsites){ die "Wrong number of site counts in SFS at $chrom:$locus->{'start'}-$locus->{'end'}!\n" }
				$thetapi1/=$denom1*$validsites;
				$thetapi2/=$denom2*$validsites;
				$thetapitot/=$denomtot*$validsites;
				$dxy/=$nalleles1*$nalleles2*$validsites;
				my $fst=($thetapitot-($thetapi1+$thetapi2)/2)/$thetapitot;
				my $df=($ref_2dsfs->[$nalleles1][0]+$ref_2dsfs->[0][$nalleles2])/$validsites;
				print $output "\t$thetapi1\t$thetapi2\t$thetapitot\t$dxy\t$fst\t$df";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub calculateSSWithinPops{
	my ($self, $genotypes, $ref_missing, $groups, $mincov, $minprop, $minind, $haplotize, $MAFfilter, $exclnonbiallelic)=@_;

	my %pops;
	for my $ref_group (@$groups){ ++$pops{$_} foreach (@$ref_group) }
	foreach (keys %pops){ warn "Population $_ occurs more than once!\n" if ($pops{$_}>1) }

	my $npops=$ref_missing->getNPop();
	print STDERR "npops: $npops\n";
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		print STDERR "$temp, $rt\n";
		}
	my @thetawdenom;
	my $tempdenom=0;
	for my $i (1..2*$rt){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }

	CHROM: for my $chrom ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			$end=$self->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;

			for my $group ( 0..$#{$groups} ){
				my $npops=@{ $groups->[$group] };
				my $validpos=0;
				my @nsnps=(0) x $npops;
				my @pi=(0) x $npops;
				my @thetaw=(0) x $npops;
				my @meanhets=(0) x $npops;
				POSITION: for my $pos ($start..$end){
					my $offset=$pos-$firststart;
					my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
					for my $subpop (@{$groups->[$group]}){
						my @alleles=(0,0,0,0);
						my ($nsamples, $nhets)=(0,0);
						my $validind=0;
						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.' && $allele<=3){ next IND }
									++$alleles[$allele];
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.' && $allele1<=3){ next IND }
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.' && $allele2<=3 ){ next IND }
									++$alleles[$allele1];
									++$alleles[$allele2];
									++$nhets if ($allele1 ne $allele2);
									}
								}
							++$validind;
							}
						if ($validind<$minind->[$subpop]){ next POSITION }
						if ($haplotize){ $nsamples=$validind }
						else { $nsamples=2*$validind }

						if ($snp){
							if ($MAFfilter){
								my $totn=0; my $minn=0; my $maxn=0;
								for my $n (@alleles){
									$totn+=$n;
									if ($minn>0){
										$minn=$n if $n>0 && $n<$minn;
										}
									else { $minn=$n }
									$maxn=$n if $n>$maxn;
									}
								if ($minn==0){ warn "No data for $chrom:$pos, pop: $subpop, skipping position!\n"; next POSITION }
								elsif ($minn==$totn){ ++$validpos; next POSITION }
								elsif ( $minn/$totn < $MAFfilter ){ next POSITION }
								}
							my $total=0;
							my $nalleles=0;
							my $denom=$nsamples*($nsamples-1);
							for my $n (@alleles){
								$total+=($n*($n-1))/$denom if ($n>1);
								++$nalleles if ($n>0);
								}
							if ($exclnonbiallelic){
								next POSITION if ($nalleles>2);
								}
							$pi[$subpop]+=(1-$total);
							if ($nalleles>1){ $thetaw[$subpop]+=1/$thetawdenom[ $nsamples-2 ]; ++$nsnps[$subpop] }
							$meanhets[$subpop]+=2*$nhets/$nsamples;	# heterozygot calls per valid individual
							}
						}
					++$validpos;
					}
				my $propcov=$validpos/$locuslength;
				if ($validpos && $propcov>=$minprop){
					print STDERR "$chrom: $start-$end, group $group has $validpos of $locuslength ($propcov) sufficiently covered sites.\n";
					for my $subpop (@{ $groups->[$group] }){
						$self->{_Loci}{$chrom}[$locus]{'covered'}[$subpop]=$validpos;
						$self->{_Loci}{$chrom}[$locus]{'nsnps'}[$subpop]=$nsnps[$subpop];
						$self->{_Loci}{$chrom}[$locus]{'pi'}[$subpop]=$pi[$subpop];
						$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$subpop]=$thetaw[$subpop];
						$self->{_Loci}{$chrom}[$locus]{'tajD'}[$subpop]=($pi[$subpop]-$thetaw[$subpop]);
						$self->{_Loci}{$chrom}[$locus]{'hets'}[$subpop]=$meanhets[$subpop];
						}
					}
				else {
					warn "$chrom: $start-$end, group $group has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
					for my $subpop (@{ $groups->[$group] }){
						$self->{_Loci}{$chrom}[$locus]{'covered'}[$subpop]=$validpos;
						$self->{_Loci}{$chrom}[$locus]{'nsnps'}[$subpop]=$nsnps[$subpop];
						$self->{_Loci}{$chrom}[$locus]{'pi'}[$subpop]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$subpop]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajD'}[$subpop]="NA";
						$self->{_Loci}{$chrom}[$locus]{'hets'}[$subpop]="NA";
						}
					}
				}
			} continue { $prevend=$end }
		}
	return;
	}

sub calculateSSWithinGroups{
	my ($self, $genotypes, $ref_missing, $groups, $mincov, $minprop, $minind, $haplotize, $MAFfilter, $exclnonbiallelic)=@_;

	my $npops=$ref_missing->getNPop();
	print STDERR "npops: $npops\n";
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		print STDERR "$temp, $rt\n";
		}
	my @thetawdenom;
	my $tempdenom=0;
	for my $i (1..2*$rt){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }

	CHROM: for my $chrom ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			$end=$self->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;

			for my $group ( 0..$#{$groups} ){
				my ($validpos, $nsnps, $pi, $thetaw, $meanhets)=(0,0,0,0,0);
				POSITION: for my $pos ($start..$end){
					my $offset=$pos-$firststart;
					my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
					my @alleles=(0,0,0,0);
					my ($nsamples, $nhets)=(0,0);
					for my $subpop (@{$groups->[$group]}){
						my $validind=0;
						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ($snp){
								if ($haplotize){
									my $index=2*$ind+int(rand(2));
									my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele && $allele ne '.' && $allele<=3){ next IND }
									++$alleles[$allele];
									}
								else {
									my $index=2*$ind;
									my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
									unless (defined $allele1 && $allele1 ne '.' && $allele1<=3){ next IND }
									my $allele2=$genotypes->getValue($chrom, $pos, $subpop, ++$index);
									unless (defined $allele2 && $allele2 ne '.' && $allele2<=3 ){ next IND }
									++$alleles[$allele1];
									++$alleles[$allele2];
									++$nhets if ($allele1 ne $allele2);
									}
								}
							++$validind;
							}
						if ($validind<$minind->[$subpop]){ next POSITION }
						if ($haplotize){ $nsamples+=$validind }
						else { $nsamples+=(2*$validind) }
						}
					if ($snp){
						if ($MAFfilter){
							my $totn=0; my $minn=0; my $maxn=0;
							for my $n (@alleles){
								$totn+=$n;
								if ($minn>0){
									$minn=$n if $n>0 && $n<$minn;
									}
								else { $minn=$n }
								$maxn=$n if $n>$maxn;
								}
							if ($minn==0){ warn "No data for $chrom:$pos, group: $group, skipping position!\n"; next POSITION }
							elsif ($minn==$totn){ ++$validpos; next POSITION }
							elsif ( $minn/$totn < $MAFfilter ){ next POSITION }
							}
						my $total=0;
						my $nalleles=0;
						my $denom=$nsamples*($nsamples-1);
						for my $n (@alleles){
							$total+=($n*($n-1))/$denom if ($n>1);
							++$nalleles if ($n>0);
							}
						if ($exclnonbiallelic){
							next POSITION if ($nalleles>2);
							}
						$pi+=(1-$total);
						if ($nalleles>1){ $thetaw+=1/$thetawdenom[ $nsamples-2 ]; ++$nsnps }
						$meanhets+=2*$nhets/$nsamples;	# heterozygot calls per valid individual
						}
					++$validpos;
					}
				my $propcov=$validpos/$locuslength;
				$self->{_Loci}{$chrom}[$locus]{'covered'}[$group]=$validpos;
				$self->{_Loci}{$chrom}[$locus]{'nsnps'}[$group]=$nsnps;
				if ($validpos && $propcov>=$minprop){
					$self->{_Loci}{$chrom}[$locus]{'pi'}[$group]=$pi;
					$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$group]=$thetaw;
					$self->{_Loci}{$chrom}[$locus]{'tajD'}[$group]=($pi-$thetaw);
					$self->{_Loci}{$chrom}[$locus]{'hets'}[$group]=$meanhets;
					print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\n";
					}
				else {
					$self->{_Loci}{$chrom}[$locus]{'pi'}[$group]="NA";
					$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$group]="NA";
					$self->{_Loci}{$chrom}[$locus]{'tajD'}[$group]="NA";
					$self->{_Loci}{$chrom}[$locus]{'hets'}[$group]="NA";
					warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
					}
				}
			} continue { $prevend=$end }
		}
	return;
	}


sub calculateSS{
	my ($self, $genotypes, $ref_missing, $groups, $mincov, $minprop, $minind, $haplotize, $MAFfilter, $exclnonbiallelic)=@_;

	my $npops=$ref_missing->getNPop();
	print STDERR "npops: $npops\n";
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		print STDERR "$temp, $rt\n";
		}
	my @thetawdenom;
	my $tempdenom=0;
	for my $i (1..2*$rt){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }

	my $poppair=0;
	for my $pop1 ( 0..($#{$groups}-1) ){
		for my $pop2 ( ($pop1+1)..$#{$groups} ){
			CHROM: for my $chrom ( $self->allKeys() ){
				my ($end, $prevend);
				my $firststart;
				LOCUS: for my $locus ( $self->allKeys($chrom) ){
					my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
					$end=$self->{_Loci}{$chrom}[$locus]{'end'};
					$firststart=$start unless defined $firststart;
					if (defined $prevend && $start-$prevend>1){ $firststart=$start }
					my $locuslength=$end-$start+1;
					my $validpos=0;
					my ($nsnpstot, $nsnps1, $nsnps2)=(0, 0, 0);
					my $pi1=0; my $pi2=0; my $dxy=0; my $da=0; my $pitot=0;
					my $thetaw1=0; my $thetaw2=0; my $thetawtot=0;
					my $sumvara=0; my $sumvart=0;

					POSITION: for my $pos ($start..$end){
						my $offset=$pos-$firststart;
						my $snp=defined $genotypes->getValue($chrom, $pos) ? 1 : 0;
						my @alleles=([0,0,0,0], [0,0,0,0]);
						my @totalleles=(0,0,0,0);
						my @nsamples=(0,0);
						my $popindex=0;
						for my $metapop ($pop1, $pop2){
							for my $subpop (@{$groups->[$metapop]}){
								my $validind=0;
								IND: for my $ind (0..$popsizes[$subpop]-1){
									my $allind=$ind+$popsum[$subpop];
									my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
									unless (defined $coverage && $coverage>=$mincov){next IND}
									if ($snp){
										if ($haplotize){
											my $index=2*$ind+int(rand(2));
											my $allele=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele && $allele ne '.' && $allele<=3){next IND}
											++$alleles[$popindex][$allele]; ++$totalleles[$allele];
											}
										else {
											my $index=2*$ind;
											my $allele1=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele1 && $allele1 ne '.' && $allele1<=3){next IND}
											++$index;
											my $allele2=$genotypes->getValue($chrom, $pos, $subpop, $index);
											unless (defined $allele2 && $allele2 ne '.' && $allele2<=3 ){next IND}
											++$alleles[$popindex][$allele1]; ++$totalleles[$allele1];
											++$alleles[$popindex][$allele2]; ++$totalleles[$allele2];
											}
										}
									++$validind;
									}
#								print STDERR "$metapop, $subpop, $chrom, $locus, $offset, $validind\n";
								if ($validind<$minind->[$subpop]){next POSITION}
								if ($haplotize){ $nsamples[$popindex]+=$validind }
								else { $nsamples[$popindex]+=(2*$validind) }
								}
							++$popindex;
							}
						if ($snp){
							if ($MAFfilter){
								my $totn=0; my $minn=0; my $maxn=0;
								for my $call (0..3){
									$totn+=$totalleles[$call];
									if ($minn>0){
										$minn=$totalleles[$call] if $totalleles[$call]>0 && $totalleles[$call]<$minn;
										}
									else { $minn=$totalleles[$call] }
									$maxn=$totalleles[$call] if $totalleles[$call]>$maxn;
									}
								if ($minn==0){ warn "No data for $chrom:$pos, $pop1 vs. $pop2, skipping position!\n"; next POSITION }
								elsif ($minn==$totn){ ++$validpos; next POSITION }
								elsif ( $minn/$totn < $MAFfilter ){ next POSITION }
								}
							my ($total1, $total2, $totalpair, $totaltot)=(0, 0, 0, 0);
							my ($nalleles1, $nalleles2, $nallelestot)=(0, 0, 0);
							my $denom1=$nsamples[0]*($nsamples[0]-1);
							my $denom2=$nsamples[1]*($nsamples[1]-1);
							my $denompair=$nsamples[0]*$nsamples[1];
							my $validtot=$nsamples[0]+$nsamples[1];
							my $denomtot=$validtot*($validtot-1);
							for my $call (0..3){
								if ($alleles[0][$call]){
									$total1+=($alleles[0][$call]*($alleles[0][$call]-1))/$denom1;
									++$nalleles1;
									}
								if ($alleles[1][$call]){
									$total2+=($alleles[1][$call]*($alleles[1][$call]-1))/$denom2;
									++$nalleles2;
									}
								if ($totalleles[$call]){
									$totaltot+=($totalleles[$call]*($totalleles[$call]-1))/$denomtot;
									++$nallelestot;
									}
								$totalpair+=($alleles[0][$call]*$alleles[1][$call])/$denompair;
								}
							if ($exclnonbiallelic){
								next POSITION if ($nallelestot>2);
								}
							$pi1+=(1-$total1);
							$pi2+=(1-$total2);
							$dxy+=(1-$totalpair);
							$da+=(1-$totalpair)-(((2-$total1-$total2))/2);
							$pitot+=(1-$totaltot);
							if ($nalleles1>1){ $thetaw1+=1/$thetawdenom[ $nsamples[0]-2 ]; ++$nsnps1 }
							if ($nalleles2>1){ $thetaw2+=1/$thetawdenom[ $nsamples[1]-2 ]; ++$nsnps2 }
							if ($nallelestot>1){ $thetawtot+=1/$thetawdenom[ $validtot-2 ]; ++$nsnpstot }
							my $ssdwp=(($nsamples[0]-1)/2)*(1-$total1)+(($nsamples[1]-1)/2)*(1-$total2);
							my $ssdtot=(($validtot-1)/2)*(1-$totaltot);
							my $sumn=$validtot-(($nsamples[0]*$nsamples[0])/$validtot+($nsamples[1]*$nsamples[1])/$validtot);
							$sumvart+=($validtot-2)*$ssdtot-($validtot-1-$sumn)*$ssdwp;
							$sumvara+=($validtot-2)*$ssdtot-($validtot-1)*$ssdwp;
							}
						++$validpos;
						}
					my $propcov=$validpos/$locuslength;
					$self->{_Loci}{$chrom}[$locus]{'covered'}[$poppair]=$validpos;
					$self->{_Loci}{$chrom}[$locus]{'nsnpstot'}[$poppair]=$nsnpstot;
					$self->{_Loci}{$chrom}[$locus]{'nsnps1'}[$poppair]=$nsnps1;
					$self->{_Loci}{$chrom}[$locus]{'nsnps2'}[$poppair]=$nsnps2;
					if ($validpos && $propcov>=$minprop){
						$self->{_Loci}{$chrom}[$locus]{'pi1'}[$poppair]=$pi1/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'pi2'}[$poppair]=$pi2/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'pitot'}[$poppair]=$pitot/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaW1'}[$poppair]=$thetaw1/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaW2'}[$poppair]=$thetaw2/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}[$poppair]=$thetawtot/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajD1'}[$poppair]=($pi1-$thetaw1)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajD2'}[$poppair]=($pi2-$thetaw2)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajDtot'}[$poppair]=($pitot-$thetawtot)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]=$dxy/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'da'}[$poppair]=$da/$validpos;
						my $piS=(0.25*$pi1 + 0.25*$pi2) / 0.5;
						my $piT=0.25*$pi1 + 0.25*$pi2 + 0.5*$dxy;
						$self->{_Loci}{$chrom}[$locus]{'fst'}[$poppair]=$piT ? 1-($piS/$piT) : "NA";
						$self->{_Loci}{$chrom}[$locus]{'phist'}[$poppair]=$sumvart ? $sumvara/$sumvart : "NA";
						print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\n";
						}
					else {
						$self->{_Loci}{$chrom}[$locus]{'pi1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'pi2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'pitot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaW1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaW2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajD1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajD2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajDtot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'da'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'fst'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'phist'}[$poppair]="NA";
						warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
						}
					} continue { $prevend=$end }
				}
			++$poppair;
			}
		}
	return;
	}


sub calculateSSfreq{
	my ($self, $ref_frequencies, $ref_popsizes, $groups, $minprop, $minind, $MAFfilter)=@_;

	my @thetawdenom;
	my $tempdenom=0;
	my $popsum=Misc::sum(@$ref_popsizes);
	for my $i (1..2*$popsum){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }
	my $limit=defined $MAFfilter ? $MAFfilter : 0;

	my $poppair=0;
	for my $pop1 ( 0..($#{$groups}-1) ){
		for my $pop2 ( ($pop1+1)..$#{$groups} ){
			CHROM: for my $chrom ( $self->allKeys() ){
				LOCUS: for my $locus ( $self->allKeys($chrom) ){
					my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
					my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
					my $locuslength=$end-$start+1;
					my $validpos=0;
					my ($nsnpstot, $nsnps1, $nsnps2)=(0, 0, 0);
					my ($pi1, $pi2, $dxy, $da, $pitot, $afd, $derplus1, $derplus2)=(0) x 8;
					my $thetaw1=0; my $thetaw2=0; my $thetawtot=0;
					my $sumvara=0; my $sumvart=0;

					POSITION: for my $pos ($start..$end){
						my $ref_maf=$ref_frequencies->getValue($chrom, $pos);
						my (@maf, @nsamples);
						for my $metapop ($pop1, $pop2){
							my ($minorind, $metaind);
							for my $subpop (@{$groups->[$metapop]}){
								my $popind=$ref_frequencies->getNSamples($chrom, $pos, $subpop);
								if ($popind<$minind->[$subpop]){ next POSITION }
								$metaind+=$popind;
								my $popfreq=defined $ref_maf->[$subpop] ? $ref_maf->[$subpop] : 0;
								$minorind+=$popfreq*$popind;
								}
							push @maf, $minorind/$metaind;
							push @nsamples, $metaind;
							}
						if (defined $ref_maf){
							die "Need a minimum of two individuals per metapopulation!\n" if ($nsamples[0]<2 || $nsamples[1]<2);
							$maf[0]=0 if ($maf[0]<$limit);
							$maf[1]=0 if ($maf[1]<$limit);
							my $snp1=($maf[0]>0 && $maf[0]<1) ? 1 : 0;
							my $snp2=($maf[1]>0 && $maf[1]<1) ? 1 : 0;
							my $snptot=($snp1 || $snp2 || $maf[0]!=$maf[1]) ? 1 : 0;
							my $multi1=$nsamples[0]/($nsamples[0]-1);
							my $multi2=$nsamples[1]/($nsamples[1]-1);
							my $validtot=$nsamples[0]+$nsamples[1];
							my $multitot=$validtot/($validtot-1);
							my $q1=1-$maf[0];
							my $q2=1-$maf[1];
							my $ptot=($maf[0]*$nsamples[0]+$maf[1]*$nsamples[1])/$validtot;
							my $qtot=1-$ptot;
							my $temppi1=$multi1*(1-( $maf[0]**2+$q1**2 ) );
							my $temppi2=$multi2*(1-( $maf[1]**2+$q2**2 ) );
							my $temppitot=$multitot*(1-( $ptot**2+$qtot**2 ) );
							my $tempdxy=$maf[0]*$q2+$maf[1]*$q1;
							$pi1+=$temppi1;
							$pi2+=$temppi2;
							$pitot+=$temppitot;
							$dxy+=$tempdxy;
							$da+=$tempdxy-( ($temppi1+$temppi2)/2 );
							$afd+=abs($maf[0]-$maf[1]) if ($snptot);
							$derplus1+=$maf[0]-$maf[1] if ($maf[0]>$maf[1]);
							$derplus2+=$maf[1]-$maf[0] if ($maf[1]>$maf[0]);
							$thetaw1+=1/$thetawdenom[ $nsamples[0]-2 ] if $snp1;
							$thetaw2+=1/$thetawdenom[ $nsamples[1]-2 ] if $snp2;
							$thetawtot+=1/$thetawdenom[ $validtot-2 ] if ($snp1 || $snp2);
							my $ssdwp=(($nsamples[0]-1)/2)*$temppi1+(($nsamples[1]-1)/2)*$temppi2;
							my $ssdtot=(($validtot-1)/2)*$temppitot;
							my $sumn=$validtot-(($nsamples[0]*$nsamples[0])/$validtot+($nsamples[1]*$nsamples[1])/$validtot);
							$sumvart+=($validtot-2)*$ssdtot-($validtot-1-$sumn)*$ssdwp;
							$sumvara+=($validtot-2)*$ssdtot-($validtot-1)*$ssdwp;
							++$nsnps1 if $snp1;
							++$nsnps2 if $snp2;
							++$nsnpstot if ($snptot);
							}
						++$validpos;
						}
					my $propcov=$validpos/$locuslength;
					$self->{_Loci}{$chrom}[$locus]{'covered'}[$poppair]=$validpos;
					$self->{_Loci}{$chrom}[$locus]{'nsnpstot'}[$poppair]=$nsnpstot;
					$self->{_Loci}{$chrom}[$locus]{'nsnps1'}[$poppair]=$nsnps1;
					$self->{_Loci}{$chrom}[$locus]{'nsnps2'}[$poppair]=$nsnps2;
					if ($validpos && $propcov>=$minprop){
						$self->{_Loci}{$chrom}[$locus]{'pi1'}[$poppair]=$pi1/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'pi2'}[$poppair]=$pi2/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'pitot'}[$poppair]=$pitot/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaW1'}[$poppair]=$thetaw1/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaW2'}[$poppair]=$thetaw2/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}[$poppair]=$thetawtot/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajD1'}[$poppair]=($pi1-$thetaw1)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajD2'}[$poppair]=($pi2-$thetaw2)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'tajDtot'}[$poppair]=($pitot-$thetawtot)/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]=$dxy/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'da'}[$poppair]=$da/$validpos;
						$self->{_Loci}{$chrom}[$locus]{'afd'}[$poppair]=$afd;
						$self->{_Loci}{$chrom}[$locus]{'derplus1'}[$poppair]=$derplus1;
						$self->{_Loci}{$chrom}[$locus]{'derplus2'}[$poppair]=$derplus2;
						my $piS=(0.25*$pi1 + 0.25*$pi2) / 0.5;
						my $piT=0.25*$pi1 + 0.25*$pi2 + 0.5*$dxy;
						eval { $self->{_Loci}{$chrom}[$locus]{'fst'}[$poppair]=1-($piS/$piT) };
						eval { $self->{_Loci}{$chrom}[$locus]{'phist'}[$poppair]=$sumvara/$sumvart };
						print STDERR "$chrom: $start-$end\t$validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites.\n";
						}
					else {
						$self->{_Loci}{$chrom}[$locus]{'pi1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'pi2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'pitot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaW1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaW2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajD1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajD2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'tajDtot'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'dxy'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'da'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'afd'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'derplus1'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'derplus2'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'fst'}[$poppair]="NA";
						$self->{_Loci}{$chrom}[$locus]{'phist'}[$poppair]="NA";
						warn "$chrom: $start-$end has only $validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites!\n";
						}
					}
				}
			++$poppair;
			}
		}
	return;
	}


sub calculatePrivatePi{
	my ($self, $ref_frequencies, $ref_popsizes, $groups, $minprop, $minind, $MAFfilter)=@_;

	my @thetawdenom;
	my $tempdenom=0;
	my $popsum=Misc::sum(@$ref_popsizes);
	for my $i (1..2*$popsum){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }
	my $limit=defined $MAFfilter ? $MAFfilter : 0;

	CHROM: for my $chrom ( $self->allKeys() ){
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $validpos=0;
			for my $group ( 0..$#{ $groups } ){
				$self->{_Loci}{$chrom}[$locus]{'pi'}[$group]=0;
				$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$group]=0;
				$self->{_Loci}{$chrom}[$locus]{'nsnps'}[$group]=0;
				}

			POSITION: for my $pos ($start..$end){
				my $ref_maf=$ref_frequencies->getValue($chrom, $pos);
				next POSITION unless (defined $ref_maf);
				my ($nsamples, $derfreq)=(0, 0);
				my @snps;
				my ($allsnp, $allder, $privatesnp)=(0, 1, 0);
				my $polygroup;
				GROUP: for my $metapop ( 0..$#{ $groups } ){
					my ($derind, $metaind)=(0, 0);
					for my $subpop (@{ $groups->[$metapop] }){
						my $popind=$ref_frequencies->getNSamples($chrom, $pos, $subpop);
						if ($popind<$minind->[$subpop]){ next POSITION }
						$metaind+=$popind;
						my $popfreq=defined $ref_maf->[$subpop] ? $ref_maf->[$subpop] : 0;
						$derind+=$popfreq*$popind;
						}
					my $groupderfreq=$derind/$metaind;
					$snps[$metapop]=($groupderfreq>0 && $groupderfreq<1) ? 1 : 0;
					$allsnp=1 if ($groupderfreq>0);
					$allder=0 if ($groupderfreq<1);
					if ($groupderfreq>$limit){
						if ($derfreq>0){ $privatesnp=0; next GROUP }	# SNP not private, but go on to check if enough samples.
						$privatesnp=1;
						$derfreq=$groupderfreq;
						$nsamples=$metaind;
						$polygroup=$metapop;
						}
					}

				if ($privatesnp){
#					print STDERR "$chrom:$pos - $polygroup: $derfreq\n";
					my $snp=($derfreq>0 && $derfreq<1) ? 1 : 0;
					my $multi=$nsamples/($nsamples-1);
					my $q=1-$derfreq;
					$self->{_Loci}{$chrom}[$locus]{'pi'}[$polygroup]+=$multi*(1-( $derfreq**2+$q**2 ) );
					$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$polygroup]+=1/$thetawdenom[$nsamples-2] if $snp;
					++$self->{_Loci}{$chrom}[$locus]{'nprsnps'}[$polygroup] if $snp;
					}
				for my $group (0..$#snps){
					++$self->{_Loci}{$chrom}[$locus]{'nsnps'}[$group] if $snps[$group];
					}
				++$self->{_Loci}{$chrom}[$locus]{'nsnpstot'} if ($allsnp && !$allder);
				++$validpos;
				}

			my $propcov=$validpos/$locuslength;
			$self->{_Loci}{$chrom}[$locus]{'covered'}=$validpos;

			if ($validpos && $propcov>=$minprop){
				$self->{_Loci}{$chrom}[$locus]{'pi'}[$_]/=$validpos foreach ( 0..$#{ $groups } );
				$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$_]/=$validpos foreach ( 0..$#{ $groups } );
				print STDERR "$chrom: $start-$end\t$validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites.\n";
				}
			else {
				$self->{_Loci}{$chrom}[$locus]{'pi'}[$_]="NA" foreach ( 0..$#{ $groups } );
				$self->{_Loci}{$chrom}[$locus]{'thetaW'}[$_]="NA" foreach ( 0..$#{ $groups } );
				warn "$chrom: $start-$end has only $validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites!\n";
				}
			}
		}
	return;
	}


sub printPrivateSNPs{
	my ($self, $outfolder, $outprefix, $ref_frequencies, $ref_popsizes, $ref_groups, $ref_groupnames, $minprop, $minind, $MAFfilter)=@_;

	my $limit=defined $MAFfilter ? $MAFfilter : 0;
	my @outputs;
	for my $metapop ( 0..$#{ $ref_groups } ){
		open $outputs[$metapop], ">>", "$outfolder/${outprefix}_" . $ref_groupnames->[$metapop] . ".bed" or die "Could not open output file!\n";
		}

	CHROM: for my $chrom ( $self->allKeys() ){
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $validpos=0;

			POSITION: for my $pos ($start..$end){
				my $ref_maf=$ref_frequencies->getValue($chrom, $pos);
				next POSITION unless (defined $ref_maf);
				my @derfreqs=(0) x scalar(@$ref_groups);
				my $derfreq=0;
				my $privatesnp=0;
				my $polygroup;
				GROUP: for my $metapop ( 0..$#{ $ref_groups } ){
					my ($derind, $metaind)=(0, 0);
					for my $subpop (@{ $ref_groups->[$metapop] }){
						my $popind=$ref_frequencies->getNSamples($chrom, $pos, $subpop);
						if ($popind<$minind->[$subpop]){ next POSITION }
						$metaind+=$popind;
						my $popfreq=defined $ref_maf->[$subpop] ? $ref_maf->[$subpop] : 0;
						$derind+=$popfreq*$popind;
						}
					$derfreqs[$metapop]=$derind/$metaind;
					if ($derfreqs[$metapop]>$limit){
						if ($derfreq>0){ $privatesnp=0; next GROUP }	# SNP not private, but go on to check if enough samples.
						$privatesnp=1;
						$derfreq=$derfreqs[$metapop];
						$polygroup=$metapop;
						}
					}

				if ($privatesnp){
#					print STDERR "$chrom:$pos - $polygroup: $derfreq\n";
					my $snp=($derfreq>0 && $derfreq<1) ? 1 : 0;
					if ($snp){
						print { $outputs[$polygroup] } "$chrom\t$pos\t", $pos+1;
						print { $outputs[$polygroup] } "\t", sprintf("%.4f", $_) foreach (@derfreqs);
						print { $outputs[$polygroup] } "\n";
						}
					}
				++$validpos;
				}

			my $propcov=$validpos/$locuslength;
			if ($validpos && $propcov>=$minprop){
				print STDERR "$chrom: $start-$end\t$validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites.\n";
				}
			else {
				warn "$chrom: $start-$end has only $validpos of $locuslength (", sprintf("%.2f", $propcov*100), "%) sufficiently covered sites!\n";
				}
			}
		}
	return;
	}


sub ABBA_BABA_perwindow{
	my ($self, $ref_frequencies, $ref_popsizes, $groups, $minprop, $minind, $MAFfilter)=@_;

	unless (@{$groups}==3){ die "ABBA-BABA test requires to define exactly 3 groups of populations, excluding the outgroup!\n" }
	my $npops=scalar @$ref_popsizes;
	print STDERR "npops: $npops\n";
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_popsizes->[$pop];
		push @popsizes, $temp;
		$rt+=$temp;
		print STDERR "$temp, $rt\n";
		}
	my @thetawdenom;
	my $tempdenom=0;
	for my $i (1..2*$rt){ $tempdenom+=1/$i; push @thetawdenom, $tempdenom }

	CHROM: for my $chrom ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		my $locindex=0;
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			$end=$self->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ ++$locindex; $firststart=$start }
			my $locuslength=$end-$start+1;
			my $validpos=0;
			my ($abba, $baba)=(0, 0);
			my ($nsnps, $thetaw1, $thetaw2, $thetaw3, $thetawtot)=(0, 0, 0, 0, 0);

			POSITION: for my $pos ($start..$end){
				my $ref_maf=$ref_frequencies->getValue($chrom, $pos);
				my (@maf, @nsamples);
				for my $metapop (@$groups){
					my ($minorind, $metaind);
					for my $subpop (@$metapop){
						my $popind=$ref_frequencies->getNSamples($chrom, $pos, $subpop);
						if ($popind<$minind->[$subpop]){ next POSITION }
						$metaind+=$popind;
						my $popfreq=defined $ref_maf->[$subpop] ? $ref_maf->[$subpop] : 0;
						$minorind+=$popfreq*$popind;
						}
					push @maf, $minorind/$metaind;
					push @nsamples, $metaind;
					}
				if (defined $ref_maf){
					my $limit=defined $MAFfilter ? $MAFfilter : 0;
					my $snp1=$maf[0]>$limit ? 1 : 0;
					my $snp2=$maf[1]>$limit ? 1 : 0;
					my $snp3=$maf[2]>$limit ? 1 : 0;
					my $validtot=$nsamples[0]+$nsamples[1]+$nsamples[2];
					$abba+=(1-$maf[0])*$maf[1]*$maf[2]*1;
					$baba+=$maf[0]*(1-$maf[1])*$maf[2]*1;
					if ($snp1){ $thetaw1+=1/$thetawdenom[ $nsamples[0]-2 ] }
					if ($snp2){ $thetaw2+=1/$thetawdenom[ $nsamples[1]-2 ] }
					if ($snp3){ $thetaw3+=1/$thetawdenom[ $nsamples[2]-2 ] }
					if ($snp1 || $snp2 || $snp3){ ++$nsnps; $thetawtot+=1/$thetawdenom[ $validtot-2 ] }
					}
				++$validpos;
				}
			my $propcov=$validpos/$locuslength;
			$self->{_Loci}{$chrom}[$locus]{'covered'}=$validpos;
			$self->{_Loci}{$chrom}[$locus]{'nsnps'}=$nsnps;
			if ($propcov>=$minprop){
				$self->{_Loci}{$chrom}[$locus]{'abba'}=$abba/$validpos;
				$self->{_Loci}{$chrom}[$locus]{'baba'}=$baba/$validpos;
				$self->{_Loci}{$chrom}[$locus]{'thetaW1'}=$thetaw1/$validpos;
				$self->{_Loci}{$chrom}[$locus]{'thetaW2'}=$thetaw2/$validpos;
				$self->{_Loci}{$chrom}[$locus]{'thetaW3'}=$thetaw3/$validpos;
				$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}=$thetawtot/$validpos;
				print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\n";
				}
			else {
				$self->{_Loci}{$chrom}[$locus]{'abba'}="NA";
				$self->{_Loci}{$chrom}[$locus]{'baba'}="NA";
				$self->{_Loci}{$chrom}[$locus]{'thetaW1'}="NA";
				$self->{_Loci}{$chrom}[$locus]{'thetaW2'}="NA";
				$self->{_Loci}{$chrom}[$locus]{'thetaW3'}="NA";
				$self->{_Loci}{$chrom}[$locus]{'thetaWtot'}="NA";
				warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
				}
			} continue { $prevend=$end }
		}
	return;
	}


sub MQperWindow{
	my ($self, $ref_missing, $mincov, $minprop, $minind)=@_;

	my $npops=$ref_missing->getNPop();
	print STDERR "npops: $npops\n";
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		print STDERR "$temp, $rt\n";
		}

	CHROM: for my $chrom ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		my $locindex=0;
		LOCUS: for my $locus ( $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			$end=$self->{_Loci}{$chrom}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ ++$locindex; $firststart=$start }
			my $locuslength=$end-$start+1;
			my $mq_rt=0; my $validpos=0;
			POSITION: for (my $pos=$start; $pos<=$end; ++$pos){
				my $offset=$pos-$firststart;
				if (defined $minind){
					POP: for my $pop ( 0..($npops-1) ){
						next POP unless ($minind->[$pop]);
						my $validind=0;
						IND: for my $ind (0..$popsizes[$pop]-1){
							my $allind=$ind+$popsum[$pop];
							my $coverage=$ref_missing->getValue($chrom, $firststart, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							++$validind;
							}
						if ($validind<$minind->[$pop]){ next POSITION }
						}
					}
				my $mq=$ref_missing->getMQ($chrom, $locindex, $offset);
				unless (defined $mq && $mq){ next POSITION }
				$mq_rt+=$mq; ++$validpos;
				}
			my $propcov=$validpos/$locuslength;
			if ($propcov>=$minprop){
				eval { $self->{_Loci}{$chrom}[$locus]{'mq'}=$mq_rt/$validpos };
				print STDERR "$chrom: $start-$end\t$validpos of $locuslength ($propcov) sufficiently covered sites.\n";
				}
			else {
				$self->{_Loci}{$chrom}[$locus]{'mq'}="NA";
				warn "$chrom: $start-$end has only $validpos of $locuslength ($propcov) sufficiently covered sites!\n";
				}
			} continue { $prevend=$end }
		}
	return;
	}


sub ExonDensityperWindow{
	my ($self, $ref_exons)=@_;

	ID: for my $id ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $ref_locus ( $self->allValues($id) ){
			unless (exists $ref_exons->{_Loci}{$id}){
				print "No exons defined for scaffold $id!\n";
				$ref_locus->{'exondensity'}[0]=0;
				next LOCUS;	
				}
			my $start=$ref_locus->{'start'};
			$end=$ref_locus->{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $sumexons=0;

			for my $exon ( sort { $a->{'start'}<=>$b->{'start'} } $ref_exons->allValues($id) ){
				my $exstart=$exon->{'start'};
				my $exend=$exon->{'end'};
				next if ($exend<$start);
				last if ($exstart>$end);
				$exstart=$start if ($exstart<$start);
				$exend=$end if ($exend>$end);
				$sumexons+=$exend-$exstart+1;
				}
			if ($locuslength){ $ref_locus->{'exondensity'}[0]=$sumexons/$locuslength }
			else { $ref_locus->{'exondensity'}[0]="NA" }
			print STDERR "$id\t$start\t", $end+1, "\t", sprintf("%.3f", $sumexons/$locuslength), "\n";
			} continue { $prevend=$end }
		}
	return;
	}


sub MeanCoverageperWindow{
	my ($self, $ref_fasta, $ref_missing, $ref_groups, $excluderefN)=@_;

	my $npops=$ref_missing->getNPop();
	my $ref_popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	foreach (@$ref_popsizes){ push @popsum, $rt; $rt+=$_ }
	my @groupsizes=map { Misc::sum(@{ $ref_popsizes }[@$_]) } (@$ref_groups);
	foreach (0..$#groupsizes){ die "Zero individuals in group $_!\n" unless $groupsizes[$_] }

	ID: for my $id ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $locus ( $self->allKeys($id) ){
			my $start=$self->{_Loci}{$id}[$locus]{'start'};
			$end=$self->{_Loci}{$id}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$ref_fasta->randomAccess($id, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $id:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my ($validpos, $hardmasked)=(0, 0);
			my @groupcoverage=(0) x @$ref_groups;

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$firststart;
				my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){ ++$hardmasked; next POSITION }	# ignore position if base in reference sequence is hardmasked.
				GROUP: for my $group (0..$#{ $ref_groups }){
					my $groupsum=0;
					POP: for my $pop (@{ $ref_groups->[$group] }){
						my $allind=$popsum[$pop];
						IND: for my $ind (0..$ref_popsizes->[$pop]-1){
							$groupsum+=$ref_missing->getValue($id, $firststart, $allind++, $offset);
							}
						}
					$groupcoverage[$group]+=$groupsum/$groupsizes[$group];
					}
				++$validpos;
				}

			if ($validpos){
				$self->{_Loci}{$id}[$locus]{'coverage'}[$_]=$groupcoverage[$_]/$validpos foreach (0..$#groupcoverage);
				}
			else {
				$self->{_Loci}{$id}[$locus]{'coverage'}[$_]='NA' foreach (0..$#groupcoverage);
				}
			print STDERR "$id:$start-$end, valid positions: $validpos\n";
			print STDERR "Group $_, mean individual coverage: ", sprintf("%.2f", $self->{_Loci}{$id}[$locus]{'coverage'}[$_]), "\n" foreach (0..$#groupcoverage);
			} continue { $prevend=$end }
		}
	return;
	}	

sub SNPswithHetsperWindow{
	my ($self, $ref_genotypes, $ref_popsizes, $ref_minind, $ref_groups, $mincov, $nhetsallowed)=@_;

	ID: for my $id ( $self->allKeys() ){
		LOCUS: for my $locus ( $self->allKeys($id) ){
			my $start=$self->{_Loci}{$id}[$locus]{'start'};
			my $end=$self->{_Loci}{$id}[$locus]{'end'};
			my @grouphets=(0) x @$ref_groups;
			my @nsnps=(0) x @$ref_groups;

			POSITION: for my $pos ($start..$end){
				next POSITION unless ( defined $ref_genotypes->getValue($id, $pos) );
				GROUP: for my $group (0..$#{ $ref_groups }){
					my $groupcounts=0;
					POP: for my $pop (@{ $ref_groups->[$group] }){
						my $validind=0;
						IND: for my $ind (0..$ref_popsizes->[$pop]-1){
							my $coverage=$ref_genotypes->getCoverage($id, $pos, $pop, $ind);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							my $index=2*$ind;
							my $allele1=$ref_genotypes->getValue($id, $pos, $pop, $index);
							unless (defined $allele1 && $allele1 ne "."){ next IND }
							my $allele2=$ref_genotypes->getValue($id, $pos, $pop, ++$index);
							unless (defined $allele2 && $allele2 ne "."){ next IND }
							if ($allele1 ne $allele2){ ++$groupcounts }
							++$validind;
							}
						unless ( $validind>=$ref_minind->[$pop] ){ next GROUP }
						}
					if ($groupcounts>$nhetsallowed){ ++$grouphets[$group] }
					++$nsnps[$group];
					}
				}

			foreach (0..$#grouphets){
				if ($nsnps[$_]){ $self->{_Loci}{$id}[$locus]{'prophets'}[$_]=$grouphets[$_]/$nsnps[$_] }
				else { $self->{_Loci}{$id}[$locus]{'prophets'}[$_]='NA' }
				}
			print STDERR "$id:$start-$end:\n";
			print STDERR "Group $_, valid SNPs: $nsnps[$_], proportion of SNPs with more than $nhetsallowed heterozygotes: ",
				sprintf("%.3f", $self->{_Loci}{$id}[$locus]{'prophets'}[$_]), "\n" foreach (0..$#grouphets);
			}
		}
	return;
	}

sub ObsHetperWindow{
	my ($self, $ref_fasta, $ref_missing, $ref_genotypes, $ref_minind, $ref_groups, $ref_mingroupind, $mincov, $excluderefN)=@_;

	my $npops=$ref_missing->getNPop();
	my $ref_popsizes=$ref_missing->getNInd();
	my @popsum;
	my $rt=0;
	foreach (@$ref_popsizes){ push @popsum, $rt; $rt+=$_ }
	my @groupsizes=map { Misc::sum(@{ $ref_popsizes }[@$_]) } (@$ref_groups);
	foreach (0..$#groupsizes){ die "Zero individuals in group $_!\n" unless $groupsizes[$_] }

	ID: for my $id ( $self->allKeys() ){
		my ($end, $prevend);
		my $firststart;
		LOCUS: for my $locus ( $self->allKeys($id) ){
			my $start=$self->{_Loci}{$id}[$locus]{'start'};
			$end=$self->{_Loci}{$id}[$locus]{'end'};
			$firststart=$start unless defined $firststart;
			if (defined $prevend && $start-$prevend>1){ $firststart=$start }
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$ref_fasta->randomAccess($id, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $id:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }
			my @validpos=(0) x @$ref_groups;
			my @validsnps=(0) x @$ref_groups;
			my @grouphet=(0) x @$ref_groups;
			my $hardmasked=0;

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$firststart;
				my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){ ++$hardmasked; next POSITION }	# ignore position if base in reference sequence is hardmasked.
				GROUP: for my $group (0..$#{ $ref_groups }){
					my $groupcounts=0;
					my $validgroupind=0;
					my %alleles;
					POP: for my $pop (@{ $ref_groups->[$group] }){
						my $allind=$popsum[$pop];
						my $validind=0;
						IND: for my $ind (0..$ref_popsizes->[$pop]-1){
							my $coverage=$ref_missing->getValue($id, $firststart, $allind++, $offset);
							unless (defined $coverage && $coverage>=$mincov){ next IND }
							if ( defined $ref_genotypes->getValue($id, $pos) ){
								my $index=2*$ind;
								my $allele1=$ref_genotypes->getValue($id, $pos, $pop, $index);
								unless (defined $allele1 && $allele1 ne "."){ next IND }
								++$alleles{$allele1};
								my $allele2=$ref_genotypes->getValue($id, $pos, $pop, ++$index);
								unless (defined $allele2 && $allele2 ne "."){ next IND }
								++$alleles{$allele2};
								if ($allele1 ne $allele2){ ++$groupcounts }
								}
							++$validind;
							}
						unless ( $validind>=$ref_minind->[$pop] ){ next GROUP }
						$validgroupind+=$validind;
						}
					unless ( $validgroupind>=$ref_mingroupind->[$group] ){ next GROUP }
					$grouphet[$group]+=$groupcounts/$validgroupind;
					++$validpos[$group];
					++$validsnps[$group] if (keys %alleles>=2);
					}
				}

			foreach (0..$#grouphet){
				if ($validpos[$_]){ $self->{_Loci}{$id}[$locus]{'Ho'}[$_]=$grouphet[$_]/$validpos[$_] }
				else { $self->{_Loci}{$id}[$locus]{'Ho'}[$_]='NA' }
				if ($validsnps[$_]){ $self->{_Loci}{$id}[$locus]{'SHo'}[$_]=$grouphet[$_]/$validsnps[$_] }
				else { $self->{_Loci}{$id}[$locus]{'SHo'}[$_]='NA' }
				}
			print STDERR "$id:$start-$end:\n";
			print STDERR "Group $_, valid positions: $validpos[$_], observed heterozygosity: ",
				sprintf("%.5f", $self->{_Loci}{$id}[$locus]{'Ho'}[$_]), "\n",
				"valid SNPs: $validsnps[$_], observed site heterozygosity: ",
				sprintf("%.5f", $self->{_Loci}{$id}[$locus]{'SHo'}[$_]), "\n" foreach (0..$#grouphet);
			} continue { $prevend=$end }
		}
	return;
	}


sub heterozygosityTest{
	my ($self, $regions, $ref_fasta, $genotypes, $ref_missing, $groups, $mincov, $minind, $maxcov, $excludeRefN)=@_;

	my $npops=$ref_missing->getNPop();
	my (@popsizes, @popsum);
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		my $temp=$ref_missing->getNInd($pop);
		push @popsizes, $temp;
		$rt+=$temp;
		}

	unless ( exists $self->{_HetbyCov} ){	# to prevent overwriting of previously stored counts.
		for my $metapop ( 0..$#{$groups} ){
			$self->{_HetbyCov}[$metapop]=[ (0) x $maxcov ];
			$self->{_ValidbyCov}[$metapop]=[ (0) x $maxcov ];
#			$self->{_Singletons}[$metapop]=[ (0) x $maxcov ];
#			$self->{_Doubletons}[$metapop]=[ (0) x $maxcov ];
			}
		}
	unless ( exists $self->{_HetbyInd} ){	# to prevent overwriting of previously stored counts.
		$self->{_HetbyInd}=[ (0) x $rt ];
		$self->{_ValidbyInd}=[ (0) x $rt ];
		}

	unless ( exists $self->{_HetbyIndCov} ){	# to prevent overwriting of previously stored counts.
		$self->{_HetbyIndCov}[$_]=[ (0) x $maxcov ] foreach (0..$rt-1);
		$self->{_ValidbyIndCov}[$_]=[ (0) x $maxcov ] foreach (0..$rt-1);
		}

	CHROM: for my $chrom ( sort {Misc::expand($a) cmp Misc::expand($b)} $regions->allKeys() ){
		LOCUS: for my $locus ( $regions->allKeys($chrom) ){
			my $start=$regions->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$regions->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				if ($excludeRefN){
					next POSITION if (uc( substr($$ref_seqstring, 0, 1, "") ) eq 'N');	# if reference has a N at the current position, skip position.
					}
				for my $metapop ( 0..$#{$groups} ){
#					my ($polym, $singleton, $doubleton)=(0, 0, 0);
					my $singcov;
					for my $subpop (@{$groups->[$metapop]}){
						my @coverages=( (0) x $popsizes[$subpop] );
						my $validind=0;
						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $allind=$ind+$popsum[$subpop];
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							++$validind if (defined $coverage && $coverage>=$mincov);
							$coverages[$ind]=$coverage if (defined $coverage);
							}
						unless ( $validind>=$minind->[$subpop] ){ next POSITION }
						my $gtstring=$genotypes->getValue($chrom, $pos, $subpop);
#						print STDERR "$pos, $subpop:\t$gtstring\n" if (defined $gtstring);

						IND: for my $ind (0..$popsizes[$subpop]-1){
							my $coverage=$coverages[$ind]>$maxcov ? $maxcov : $coverages[$ind];
							unless ($coverage>0){ next IND }
							my $allind=$ind+$popsum[$subpop];
							if (defined $gtstring){
								my $index=2*$ind;
								my $allele1=substr($gtstring, $index, 1);
								unless (defined $allele1 && $allele1 ne '.'){ next IND }
								my $allele2=substr($gtstring, ++$index, 1);
								unless (defined $allele2 && $allele2 ne '.'){ next IND }
								++$self->{_AltbyCov}[$metapop][$coverage-1] if ($allele1>0);
								++$self->{_AltbyCov}[$metapop][$coverage-1] if ($allele2>0);
#								if ( $polym==0 && ($allele1>0 || $allele2>0) ){
#									$polym=1;
#									$singleton=1 if ($allele1==0 || $allele2==0);
#									$doubleton=1 if ($allele1==$allele2);
#									$singcov=$coverage;
#									}
#								elsif ( $polym==1 && ($allele1>0 || $allele2>0) ){
#									$singleton=0;
#									$doubleton=0;
#									}
								unless ($allele1 eq $allele2){
									++$self->{_HetbyCov}[$metapop][$coverage-1];
									++$self->{_HetbyInd}[$allind] if ($coverage>=$mincov);
									++$self->{_HetbyIndCov}[$allind][$coverage-1];
									}
								}
							++$self->{_ValidbyCov}[$metapop][$coverage-1];
							++$self->{_ValidbyInd}[$allind] if ($coverage>=$mincov);
							++$self->{_ValidbyIndCov}[$allind][$coverage-1];
							}
						}
#					++$self->{_Singletons}[$metapop][$singcov-1] if $singleton;
#					++$self->{_Doubletons}[$metapop][$singcov-1] if $doubleton;
					}
				++$self->{_ValidPos};
				}
			}
		}
	return;
	}


sub generateFasta{
	my ($self, $ref_fasta, $ref_fastaarray, $ref_genotypes, $ref_missing, $ref_minind, $groups, $ref_mingroupind, $mincov, $excluderefN, $haplotize, $ref_subind)=@_;

	my $verbose=0;
	my @individuals;
	my (@popsizes, @popsum, @groupsizes, @groupsum);
	my $npops=$ref_missing->getNPop();
	my $totsize=0;
	for my $group ( 0..$#{$groups} ){
		push @groupsum, $totsize;
		my $groupsize=0;
		for my $pop ( @{$groups->[$group]} ){
			push @popsum, $totsize;
			my $temp=$ref_missing->getNInd($pop);
			push @popsizes, $temp;
			if (defined $ref_subind){ push @individuals, ( map { $totsize+$_ } @{ $ref_subind->[$pop] } ) }
			else { push @individuals, ($totsize..$totsize+$temp-1) }
			$totsize+=$temp;
			$groupsize+=$temp;
			print STDERR "pop $pop: popsize $temp, popsum $popsum[$pop], totsize: $totsize\n" if $verbose;
			}
		push @groupsizes, $groupsize;
		print STDERR "group $group: groupsize $groupsize, groupsum $groupsum[$group], totsize: $totsize\n" if $verbose;
		}
	die "Wrong length of array of fasta objects provided!\n" unless (@$ref_fastaarray==@individuals);
	print STDERR join(',', @individuals), "\n" if $verbose;

	CHROM: for my $chrom ( $self->allKeys() ){
		my $prevend;
		LOCUS: for my $locus ( sort { $self->{_Loci}{$chrom}[$a]{'start'}<=>$self->{_Loci}{$chrom}[$b]{'start'} } $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			if (defined $prevend){
				if ($start-$prevend>1){ $start=$prevend+1; warn "Gap detected between two consecutive loci on the same chromosome. Start changed to $start.\n" }
				elsif (!$start-$prevend){ $start=$prevend+1; warn "Overlap detected between two consecutive loci on the same chromosome. Start changed to $start.\n"}
				}
			my $locuslength=$end-$start+1;
			my @seq;
			my $ref_seqstring;
			if (defined $ref_fasta){ $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1) }
			else { my $seqstring="A" x $locuslength; $ref_seqstring=\$seqstring }
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length at $chrom:$start-$end. ", length($$ref_seqstring), " vs. $locuslength!\n" }
#			for my $allind (0..$totsize-1){ $seq[$allind]=$$ref_seqstring }	# set the sequences of all individuals to reference.
			my ($totalpos, $hardmasked)=(0, 0);
			my $basestring;
			POSITION: for my $pos ($start..$end){
				$basestring="";
				my $refbase=substr($$ref_seqstring, 0, 1, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				my $offset=$pos-$start;
				if ($excluderefN && $refbase eq 'N'){ ++$hardmasked; $basestring="N" x $totsize; next POSITION }	# discard if base in reference sequence is hardmasked.
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				my $groupbases;
				GROUP: for my $group ( 0..$#{$groups} ){
					$groupbases="";
					my $validgroupind=0;
					POP: for my $pop ( @{$groups->[$group]} ){
						my $validind=0;
						IND: for my $ind (0..$popsizes[$pop]-1){
							my $allind=$ind+$popsum[$pop];
							my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
							unless (defined $coverage && $coverage>=$mincov){ $groupbases.="N"; next IND }	# missing data, replace reference base at offset with N.
							if ($snp){	# if position not variable but covered, retain reference base.
								my $index=2*$ind;
								my $allind=$ind+$popsum[$pop];
								my $allele1=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
								unless (defined $allele1 && $allele1 ne "N"){ $groupbases.="N"; next IND }
								my $allele2=$ref_genotypes->getBase($chrom, $pos, $pop, ++$index);
								unless (defined $allele2 && $allele2 ne "N"){ $groupbases.="N"; next IND }
								if ($allele1 eq $allele2){ $groupbases.=$allele1 }
								elsif ($haplotize){ $groupbases.=($allele1, $allele2)[ int(rand(2)) ] } 
								else { $groupbases.=Misc::iupac($allele1, $allele2) }
								}
							else { $groupbases.=$refbase }
							++$validind
							}
						if ($validind<$ref_minind->[$pop]){	# not enough individuals callable in a population, set bases to N for all individuals in the group;
							$groupbases="N" x $groupsizes[$group];
							next GROUP;
							}
						$validgroupind+=$validind;
						}
					if ($validgroupind<$ref_mingroupind->[$group]){	# not enough individuals callable in a group, set bases to N for all individuals in the group;
						$groupbases="N" x $groupsizes[$group];
						next GROUP;
						}
					die "Wrong length of groupbases string at $chrom:$pos, ", length($groupbases), " vs. $groupsizes[$group]!\n$groupbases\n" unless ( length($groupbases)==$groupsizes[$group] );
					} continue { $basestring.=$groupbases }
				die "Wrong length of basestring at $chrom:$pos, ", length($basestring), " vs. $totsize!\n$basestring\n" unless (length($basestring)==$totsize);
				} continue { $seq[$_].=substr($basestring, 0, 1, "") foreach (0..$totsize-1) }
			for my $selindex (0..$#individuals){
				my $selind=$individuals[$selindex];
				unless ( length($seq[$selind])==$locuslength){ die "Wrong length of sequence string $chrom:$start-$end of individual $selind, ",
					length($seq[$selind]), " vs. $locuslength, concatenation aborted!\n" }
				$ref_fastaarray->[$selindex]->concatSeq($chrom, \$seq[$selind]);
				}
			$prevend=$end;
			}
		}
	return;
	}


sub generateGPhocs{
	my ($self, $outfile, $ref_fasta, $ref_genotypes, $ref_missing, $ref_minind, $mincov, $excluderefN, $subsample)=@_;

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
	my $nloci=0;
	$nloci+=$self->getSize($_) foreach ( $self->allKeys() );

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">", $outfile or die "Could not open output file!\n" }
	print $output "$nloci\n\n";

	CHROM: for my $chrom ( $self->allKeys() ){
		my $prevend;
		LOCUS: for my $locus ( sort { $self->{_Loci}{$chrom}[$a]{'start'}<=>$self->{_Loci}{$chrom}[$b]{'start'} } $self->allKeys($chrom) ){
			my $start=$self->{_Loci}{$chrom}[$locus]{'start'};
			my $end=$self->{_Loci}{$chrom}[$locus]{'end'};
			my $locuslength=$end-$start+1;
			my @seq;
			my $ref_seqstring=$ref_fasta->randomAccess($chrom, $start, $end, 1);
			if (length($$ref_seqstring)!=$locuslength){ die "Retrieved reference string has incorrect length. ", length($$ref_seqstring), " vs. $locuslength!\n" }

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

			print $output "locus", ++$locusnumber, " $subtotsize $locuslength\n";			

			POSITION: for my $pos ($start..$end){
				my $offset=$pos-$start;
				my $refbase=uc( substr($$ref_seqstring, 0, 1, "") );	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
				if ($excluderefN && $refbase eq 'N'){	# if reference has a N at the current position, skip position.
					for my $pop (0..$#popsizes){ $seq[$pop][$_].="N" foreach ( 0..$#{ $indlist[$pop] } ) };
					next POSITION;
					}
				my $snp=defined $ref_genotypes->getValue($chrom, $pos) ? 1 : 0;
				POP: for my $pop (0..$#popsizes){
					my $validind=0;
					my @states;
					IND: for my $ind ( @{ $indlist[$pop] } ){
						my $allind=$ind+$popsum[$pop];
						my $coverage=$ref_missing->getValue($chrom, $start, $allind, $offset);
#						my $coverage=$ref_missing->getValue3($chrom, $pos, $allind);
						unless (defined $coverage && $coverage>=$mincov){ push(@states, "N"); next IND }	# missing data, replace reference base at offset with N.
						if ($snp){	# if position is variable and covered, obtain individual genotypes.
							my $index=2*$ind;
							my $allind=$ind+$popsum[$pop];
							my $allele1=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
							unless (defined $allele1){ push(@states, "N"); next IND }
							++$index;
							my $allele2=$ref_genotypes->getBase($chrom, $pos, $pop, $index);
							unless (defined $allele2){ push(@states, "N"); next IND }
							if ($allele1 eq $allele2){ push(@states, $allele1) }
							else { push( @states, Misc::iupac($allele1, $allele2) ) }
							}
						else {	# if position not variable but covered, retain reference base.
							push(@states, $refbase);
							}
						++$validind
						}
					unless( @states==@{ $indlist[$pop] } ){ die "Wrong length of states array! (", scalar(@states), " vs. ", scalar(@{ $indlist[$pop] }), ")\n" }
					if ($validind<$ref_minind->[$pop]){	# not enough individuals callable in a population, temporary genotypes with N for all individuals in the population;
						$seq[$pop][$_].="N" foreach (0..$#states);
						}
					else {
						$seq[$pop][$_].=$states[$_] foreach (0..$#states);
						}
					}
				}

			for my $pop (0..$#popsizes){ 
				for my $ind ( 0..$#{ $seq[$pop] } ){
					unless ( length($seq[$pop][$ind])==$locuslength){ die "Wrong length of sequence string for $chrom:$start-$end in population $pop for individual $ind!\n" }
					print $output "pop${pop}_${ind}\t$seq[$pop][$ind]\n";
					}
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}


sub printLoci{
	my ($self, $outfile, $printinfo)=@_;
	my $output;

	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	SCAFF: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		LOCUS: for my $locus (sort {$a->{'start'}<=>$b->{'start'}} @{$self->{_Loci}{$chrom}}){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t", $locus->{'info'} || '' if $printinfo;
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub printLociSSWithinPops{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'info'}" if defined $locus->{'info'};
			GROUP: for my $group (0..$#{ $locus->{'covered'} }){
				print $output "\t$locus->{'covered'}[$group]\t$locus->{'nsnps'}[$group]\t$locus->{'pi'}[$group]\t$locus->{'thetaW'}[$group]\t$locus->{'tajD'}[$group]\t$locus->{'hets'}[$group]";
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}

sub printLociSS{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

#	print $output "CHROM\tSTART\tEND";
#	for my $poppair (0..$#{ $self->{_Loci}{$chrom}[$locus]{covered} }){print $output "\tCovered\tPi1\tPi2\tDxy\tFst"}
#	print $output "\n";

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'info'}" if defined $locus->{'info'};
			POPPAIR: for my $poppair (0..$#{ $locus->{'covered'} }){
				print $output "\t$locus->{'covered'}[$poppair]\t$locus->{'nsnpstot'}[$poppair]\t$locus->{'nsnps1'}[$poppair]\t$locus->{'nsnps2'}[$poppair]\t$locus->{'pi1'}[$poppair]\t$locus->{'pi2'}[$poppair]\t$locus->{'pitot'}[$poppair]\t$locus->{'thetaW1'}[$poppair]\t$locus->{'thetaW2'}[$poppair]\t$locus->{'thetaWtot'}[$poppair]\t$locus->{'tajD1'}[$poppair]\t$locus->{'tajD2'}[$poppair]\t$locus->{'tajDtot'}[$poppair]\t$locus->{'dxy'}[$poppair]\t$locus->{'da'}[$poppair]\t$locus->{'fst'}[$poppair]\t$locus->{'phist'}[$poppair]";
				}
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}

sub printLociSSfreq{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

#	print $output "CHROM\tSTART\tEND";
#	for my $poppair (0..$#{ $self->{_Loci}{$chrom}[$locus]{covered} }){print $output "\tCovered\tPi1\tPi2\tDxy\tFst"}
#	print $output "\n";

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'info'}" if defined $locus->{'info'};
			POPPAIR: for my $poppair (0..$#{ $locus->{'covered'} }){
				print $output "\t$locus->{'covered'}[$poppair]\t$locus->{'nsnpstot'}[$poppair]\t$locus->{'nsnps1'}[$poppair]\t$locus->{'nsnps2'}[$poppair]\t$locus->{'pi1'}[$poppair]\t$locus->{'pi2'}[$poppair]\t$locus->{'pitot'}[$poppair]\t$locus->{'thetaW1'}[$poppair]\t$locus->{'thetaW2'}[$poppair]\t$locus->{'thetaWtot'}[$poppair]\t$locus->{'tajD1'}[$poppair]\t$locus->{'tajD2'}[$poppair]\t$locus->{'tajDtot'}[$poppair]\t$locus->{'dxy'}[$poppair]\t$locus->{'da'}[$poppair]\t$locus->{'afd'}[$poppair]\t$locus->{'derplus1'}[$poppair]\t$locus->{'derplus2'}[$poppair]\t$locus->{'fst'}[$poppair]\t$locus->{'phist'}[$poppair]";
				}
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub printPrivatePi{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'info'}" if defined $locus->{'info'};
			print $output "\t$locus->{'covered'}";
			print $output "\t$locus->{'nsnpstot'}";
			GROUP: for my $group (0..$#{ $locus->{'nsnps'} }){
				print $output "\t$locus->{'nsnps'}[$group]\t$locus->{'nprsnps'}[$group]\t$locus->{'pi'}[$group]\t$locus->{'thetaW'}[$group]";
				}
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub printLociABBA_BABA{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){ print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "", "\t$locus->{'covered'}\t$locus->{'nsnps'}\t$locus->{'abba'}\t$locus->{'baba'}\t$locus->{'thetaW1'}\t$locus->{'thetaW2'}\t$locus->{'thetaW3'}\t$locus->{'thetaWtot'}\n" }
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printLociTreeStats{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			next LOCUS unless ( $locus->{'covered'} );
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1;
			print $output "\t$locus->{'info'}" if defined $locus->{'info'};
			print $output "\t$locus->{'covered'}";
			if ( defined $locus->{'concordance'} ){
				for my $comp ( @{ $locus->{'concordance'} } ){
					print $output "\t$comp";
					}
				}
			if ( defined $locus->{'gsi'} ){
				for my $comp ( @{ $locus->{'gsi'} } ){
					print $output "\t$comp";
					}
				}
			if ( defined $locus->{'dist'} ){
				for my $comp ( @{ $locus->{'dist'} } ){
					print $output "\t$comp";
					}
				}
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}


sub printLociTrees{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom (sort { Misc::expand($a) cmp Misc::expand($b) } keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			next LOCUS unless ( $locus->{'tree'} );
			print $output "[ $chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "", "]\n";
			print $output "$locus->{'tree'}\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printLociMQ{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){ print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t$locus->{'info'}\t$locus->{'mq'}\n" }
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}

sub printLociCoverage{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT } else { open $output, ">>", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "";
			print $output "\t$_" foreach ( @{ $locus->{'coverage'} } );
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}

sub printLociAttribute{
	my ($self, $outfile, $attribute, $sort, $append)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT }
	elsif ($append){ open $output, ">>", $outfile or die "Could not open output file!\n" }
	else { open $output, ">", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "";
			print $output "\t$_" foreach ( @{ $locus->{$attribute} } );
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}

sub printLociAttributes{
	my ($self, $outfile, $ref_attributes, $sort, $append)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){ $output=*STDOUT }
	elsif ($append){ open $output, ">>", $outfile or die "Could not open output file!\n" }
	else { open $output, ">", $outfile or die "Could not open output file!\n" }

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if ($sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t", $locus->{'info'} || "";
			for my $group (0..$#{ $locus->{$ref_attributes->[0]} }){
				print $output "\t$locus->{$_}[$group]" foreach (@$ref_attributes);
				}
			print $output "\n";
			}
		}
	if ($outfile ne "-" && $outfile ne "STDOUT"){ close $output }
	return;
	}

sub printHeterozygosities{
	my ($self, $outfile)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	for my $allind ( 0..$#{ $self->{_ValidbyInd} } ){
		my $valid=$self->{_ValidbyInd}[$allind];
		my $het=defined $valid && $valid>0 ? $self->{_HetbyInd}[$allind] / $valid : "NA";
		print $output "$allind\t$het\t$self->{_ValidbyInd}[$allind]\n";
		}
	my $validpos=$self->{_ValidPos};
	for my $coverage ( 0..$#{ $self->{_ValidbyCov}[0] } ){
		print $output $coverage+1;
		for my $metapop ( 0..$#{ $self->{_ValidbyCov} } ){
			my $valid=$self->{_ValidbyCov}[$metapop][$coverage];
			my ($het, $alt, $singletons, $doubletons);
			if (defined $valid && $valid>0){
				$het=$self->{_HetbyCov}[$metapop][$coverage] / $valid;
				$alt=defined $self->{_AltbyCov}[$metapop][$coverage] ? $self->{_AltbyCov}[$metapop][$coverage] / $valid : "NA";
				}
			else {
				$het="NA";
				$alt="NA";
				}
#			if (defined $validpos && $validpos>0){
#				$singletons=$self->{_Singletons}[$metapop][$coverage] / $validpos;
#				$doubletons=$self->{_Doubletons}[$metapop][$coverage] / $validpos;
#				}
#			else {
#				$singletons="NA";
#				$doubletons="NA";
#				}
#			print $output "\t$het\t$alt\t$singletons\t$doubletons\t$self->{_ValidbyCov}[$metapop][$coverage]";
			print $output "\t$het\t$alt\t$self->{_ValidbyCov}[$metapop][$coverage]";
			}
		print $output "\n";
		}

	for my $coverage ( 0..$#{ $self->{_ValidbyIndCov}[0] } ){
		print $output $coverage+1;
		for my $allind ( 0..$#{ $self->{_ValidbyIndCov} } ){
			my $valid=$self->{_ValidbyIndCov}[$allind][$coverage];
			my $het=(defined $valid && $valid>0) ? $self->{_HetbyIndCov}[$allind][$coverage] / $valid : "NA";
			print $output "\t$het\t$self->{_ValidbyIndCov}[$allind][$coverage]";
			}
		print $output "\n";
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printLociSFS{
	my ($self, $outfile, $sort)=@_;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">>", $outfile or die "Could not open output file!\n"}

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		my @loci;
		if (defined $sort && $sort){ @loci=sort { $a->{'start'}<=>$b->{'start'} } @{ $self->{_Loci}{$chrom} } }
		else { @loci=@{ $self->{_Loci}{$chrom} } }
		LOCUS: for my $locus (@loci){
			print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\t$locus->{'info'}";
			GROUP: for my $metapop (0..$#{ $locus->{'SFS'} }){
				print $output "\t$locus->{'covered'}[$metapop]";
				my @totcounts;
				CAT: for my $cat ( @{ $locus->{'SFS'}[$metapop] } ){
					print $output "\t", join(',', @{$cat}[ 1..$#{$cat}-1 ] );	# omit 0 and 1 categories.
					$totcounts[$_-1]+=$cat->[$_] foreach (1..$#{$cat}-1);
					}
				print $output "\t", join(',', @totcounts);	# print overall sfs
				}
			print $output "\n";
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}
	return;
	}


sub printSubset{
	my ($self, $outfile, $limit)=@_;
	unless (keys %{$self->{_Loci}}){ return 0 }

	my $output;
	my $locuscount=0;

	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	while ($locuscount <= $limit){
		SCAFF: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
			LOCUS: for my $locus (sort {$a->{'start'}<=>$b->{'start'}} @{$self->{_Loci}{$chrom}}){
				print $output "$chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\n";
				$locuscount++;
				}
			delete $self->{_Loci}{$chrom};
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return 1;
	}


sub printSummary{ # under construction!
	my ($self, $fasta, $chrommap, $outfile)=@_;

	my $scaff_loci=$self->translateChromtoScaff($chrommap);
	$fasta->extractRegions($scaff_loci);
	my %seqs;

	my ($seq, $lower, $upper);
	my $seqname=$seq . "_start_" . $lower . "_end_" . $upper;

	my $output;
	if ($outfile eq "-" || $outfile eq "STDOUT"){$output=*STDOUT} else {open $output, ">", $outfile or die "Could not open output file!\n"}

	my $l=0;
	SCAFF: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} keys %{$self->{_Loci}}){
		LOCUS: for my $locus (sort {$a->{'start'}<=>$b->{'start'}} @{$self->{_Loci}{$chrom}}){
			print $output "Locus $l: $chrom\t$locus->{'start'}\t", $locus->{'end'}+1, "\n";
			print $output "Sequence\n";
			$l++;
			}
		}

	if ($outfile ne "-" && $outfile ne "STDOUT"){close $output}

	return;
	}



1;

