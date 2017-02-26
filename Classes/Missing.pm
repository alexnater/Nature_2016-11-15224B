package Missing;

use strict;
use warnings;
use Classes::Genotypes;
use Classes::Misc;


sub new{
	my ($class, $allsitesfolder, $allsiteslist, $simulated, $ref_samples)=@_;
	my $self={};
	$self->{_Simulated}=$simulated ? 1 : 0;
	$self->{_Pops}=defined $ref_samples ? $ref_samples : [];

	if (defined $allsitesfolder && defined $allsiteslist){
		my ($ref_indnames, $ref_files, $ref_samples)=Misc::getIndnames($allsitesfolder, $allsiteslist);
		$self->{_Files}=$ref_files;
		$self->{_Samples}=$ref_samples;
		$self->{_Pops}=[ map { scalar(@$_) } @$ref_samples ];
		$self->{_Indnames}=$ref_indnames;
		}
	else { warn "No list and folder of allsites vcf files provided!\n" }

	bless $self, $class;
	return $self;
	}


sub setPops{
	my ($self, $ref_subpops, $ref_subsamples)=@_;

	my $files=$self->{_Files};
	my $samples=$self->{_Samples};
	unless (defined $files && defined $samples){ die "List and folder of allsites vcf files have not been defined!\n" }
	$samples=$ref_subsamples if (defined $ref_subsamples);
	my ($subfiles, $subsamples, $popsizes);

	for my $i (@$ref_subpops){
		push @{ $subfiles }, $files->[$i];
		push @{ $subsamples }, $samples->[$i];
		push @{ $popsizes }, scalar(@{ $samples->[$i] });
		print STDERR "$i - $files->[$i]: ", join( ',', @{ $samples->[$i] } ), "\n"; 
		}
	$self->{_Files}=$subfiles;
	$self->{_Samples}=$subsamples;
	$self->{_Pops}=$popsizes;
	return;
	}


sub getValue{
	my ($self, $id, $startpos, $allind, $offset)=@_;

	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==4){ return($self->{_Missing}{$id}{$startpos}[$allind]) }
	elsif (@_==5){
		my $substring=substr($self->{_Missing}{$id}{$startpos}[$allind], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $startpos, $allind, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub getValue2{	# old method.
	my ($self, $id, $locus, $allind, $offset)=@_;

	if ( $self->{_Simulated} ){ return 255 }
	elsif (@_==4){ return($self->{_Missing}{$id}[$locus][$allind]) }
	elsif (@_==5){
		my $substring=substr($self->{_Missing}{$id}[$locus][$allind], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $locus, $allind, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub getValuebyPos{	# to access missing data acquired with readAllsitesBlocks
	my ($self, $id, $pos, $allind, $end)=@_;
	my $offset=$pos % 100;
	my $startpos=$pos-$offset;
	unless ( exists $self->{_Missing}{$id}{$startpos} ){ warn "Invalid access: $id, $pos!\n" }
	if (@_==3){
		if (@_==3 && $offset){ die "Individual not specified, but position is not divisible by 100!\n" }
		else { return( $self->{_Missing}{$id}{$startpos} ) }
		}
	elsif (@_==4){
		my $substring=substr( $self->{_Missing}{$id}{$startpos}[$allind], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $pos, $allind!\n" }
		return( ord($substring) );
		}
	elsif (@_==5){	# returns string of characters from start to end.
		my $string;
		my $overhang=100-$offset;
		my $strlen=$end-$pos+1;
		if ($offset>0){
			my $sublength=$strlen<$overhang ? $strlen : $overhang;
			$string.=substr($self->{_Missing}{$id}{$startpos}[$_], -$overhang, $sublength);
			$startpos+=100;
			$strlen-=$overhang;
			}
		while ($strlen>0){
			my $sublength=$strlen<100 ? $strlen : 100;
			$string.=substr($self->{_Missing}{$id}{$startpos}[$_], 0, $sublength);
			$startpos+=100;
			$strlen-=100;
			}
		unless ( length($string)==($end-$pos+1) ){ warn "Invalid access: $id:$pos-$end, $allind!\n" }
		return($string); 
		}
	else { die "Invalid number of arguments!\n" }
	}


sub setValue{
	my ($self, $id, $startpos, $allind, $offset, $value)=@_;
	unless (@_==6){die "Wrong number of arguments!\n"}

	if ($value>255){$value=255}
	substr($self->{_Missing}{$id}{$startpos}[$allind], $offset, 1, chr($value) );
	return($value);
	}


sub setValue2{	# old method.
	my ($self, $id, $locus, $allind, $offset, $value)=@_;
	unless (@_==6){die "Wrong number of arguments!\n"}

	if ($value>255){$value=255}
	substr($self->{_Missing}{$id}[$locus][$allind], $offset, 1, chr($value) );
	return($value);
	}


sub setValue3{	# to set missing data acquired with readAllsitesBlocks
	my ($self, $id, $pos, $allind, $value)=@_;
	unless (@_==5){ die "Wrong number of arguments!\n" }

	if ($value>255){ $value=255 }
	my $offset=$pos % 100;
	my $startpos=$pos-$offset;
	unless ( exists $self->{_Missing}{$id}{$startpos} ){ warn "Invalid access: $id, $pos!\n" }
	substr($self->{_Missing}{$id}{$startpos}[$allind], $offset, 1, chr($value) );
	return($value);
	}


sub getSize{
	my ($self, $id, $pos, $pop)=@_;

	if (@_==1){ return(scalar(keys %{$self->{_Missing}})) }
	elsif (@_==2){ return(scalar(keys %{$self->{_Missing}{$id}})) }
	elsif (@_==3){ return(scalar @{$self->{_Missing}{$id}{$pos}}) }
	elsif (@_==4){ return(scalar @{$self->{_Missing}{$id}{$pos}[$pop]}) }
	else { die "Invalid number of arguments!\n" }
	}


sub getPosition{
	my ($self, $id)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}

	my $position=each %{$self->{_Missing}{$id}};
	keys %{$self->{_Missing}{$id}};

	return($position);
	}


sub getNPop{
	my $self=shift;
	return(scalar @{$self->{_Pops}});
	}


sub getNInd{
	my ($self, $pop)=@_;
	if (@_==1){ return($self->{_Pops}) }
	elsif (@_==2){ return($self->{_Pops}[$pop]) }
	else { die "Wrong number of arguments!\n" }
	}


sub getIndname{
	my ($self, $allind)=@_;
	unless (@_==2){die "Wrong number of arguments!\n"}
	return($self->{_Indname}[$allind]);
	}


sub getMQ{
	my ($self, $id, $locus, $offset)=@_;

	if (@_==3){ return($self->{_MQ}{$id}[$locus]) }
	elsif (@_==4){
		my $substring=substr($self->{_MQ}{$id}[$locus], $offset, 1);
		unless (defined $substring){ warn "Invalid access: $id, $locus, $offset!\n" }
		return( ord($substring) );
		}
	else { die "Invalid number of arguments!\n" }
	}


sub readPileup{ # Currently not valid!
	my $self=shift @_;
	my $loci_file=shift @_;
	my $bamfolder=shift @_;
	my $bamlist=shift @_;
	my $refgenom=shift @_;
	my $minmq=shift @_;
	my $minvq=shift @_;
	my $npops=shift @_;
	my $popsize=shift @_;
	my $totsize=Misc::sum(@$popsize);

	my @bamfiles;
	open my $fh, "<", $bamlist or die "Could not open list of bam files $bamlist!\n";
	while (<$fh>){
		chomp;
		push @bamfiles, "$bamfolder/$_";
		}
	close $fh;

	my $command="/home/aim/anater/bin/samtools mpileup -C 50 -Q 15 -D -l $loci_file -uf $refgenom @bamfiles | /home/aim/anater/bin/bcftools view -cg -";

	my %missing;
	my @header;
	my @samples;
	my $nsamples;
	my ($mq, @label, $dp);
	my @line;
	my ($scaffold, $pos, $prevscaff, $prevpos);
	my @individuals;
	my $switch=0;

	open my $pileup, "-|", $command or die "Could not open pileup file!\n";

	READ: while (<$pileup>){
		chomp;

		if(/^#CHROM/){
			@header=split("\t", $_);
			@samples=@header[9..$#header];
			$switch=1;
			}

		elsif(!$switch){next READ}

		elsif(/^\w+\s/){
			@line=split("\t", $_);
			$scaffold=$line[0];
			$pos=$line[1]-1; # The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.

			if ($line[7]=~/INDEL/){next READ}; # warn "Warning: $scaffold: $pos contains an INDEL!\n"; 
#			print STDERR "$line[7]\n";

			$nsamples=$#line-8;
			if ($nsamples!=$totsize){die "Number or samples in pileup does not match information on command line!\n"}

			if ($line[5] eq '.' || $line[5]<$minvq){
				push @{$missing{$scaffold}{$pos}}, ( (0) x $totsize );
				next READ;
				}

			($mq)=($line[7]=~/MQ=(\d+)/);
			unless (defined $mq && $mq>=$minmq){
				push @{$missing{$scaffold}{$pos}}, ( (0) x $totsize );
				next READ;
				}

			if ($line[4] eq '.'){$dp=1} else {$dp=2}
#			@label=split(":", $line[8]);
#			LABEL: for my $lab (0..$#label){if ($label[$lab] eq "DP"){$dp=$lab; last LABEL}}

			if (defined $prevpos && $pos-$prevpos!=1 && $scaffold eq $prevscaff && $pos-$prevpos<10000){
#				warn "Gap inbetween positions $prevscaff: $prevpos and $scaffold: $pos!\n";
				for my $gap ($prevpos+1..$pos-1){push @{$missing{$scaffold}{$gap}}, ( (0) x $totsize )};
				}
			$prevscaff=$scaffold;
			$prevpos=$pos;

			@individuals=@line[-$nsamples..-1];
			for my $allind (0..$totsize-1){
				my @data=split(":", (shift @individuals));
				push @{$missing{$scaffold}{$pos}}, $data[$dp];
#				if ($data[$dp]<$mincov){push @{$missing{$scaffold}{$pos}[$pop]}, (0, 0)} else {push @{$missing{$scaffold}{$pos}[$pop]}, (1, 1)}
				}
			}
		}

	close $pileup;

	$self->{_Missing}=\%missing;
	$self->{_Header}=\@header;
	$self->{_Pops}=$popsize;
	return;
	}


sub readRawVcf{	# This method is missing proper tests to make sure that missing data strings have the correct size, especially when loci start or end with ref-Ns!
	my $self=shift @_;
	my $loci_file=shift @_;
	my $ref_loci=shift @_;
	my $bamfolder=shift @_;
	my $bamlist=shift @_;
	my $refgenom=shift @_;
	my $minmq=shift @_;
	my $minvq=shift @_;
	my $npops=shift @_;
	my $popsize=shift @_;
	my $sampleorder=shift @_;
	my $totsize=Misc::sum(@$popsize);

	my %missing; # $missing{$chrom/$scaff}[$locus][$pop][$ind]="";
	my %offsets;
#	my $bams=join(" -I ", @bamfiles);
	my $gatk="java -Xmx8g -jar /sw/apps/bioinfo/GATK/2.7.2/GenomeAnalysisTK.jar";
	my $command="$gatk -T UnifiedGenotyper -nt 2 -nct 4 -R $refgenom -I $bamlist -glm SNP -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES -L $loci_file -log GATK_log.txt";
	# -Djava.io.tmpdir=/proj/b2010010/nobackup/alexn/temp/JAVA -stand_emit_conf $minvq -stand_call_conf $minvq -out_mode EMIT_ALL_CONFIDENT_SITES

	my (@header, @samples, @line, @individuals);
	my ($nsamples, $mq, $vq, @label, $dp, $snp);
	my ($poscol, $altcol, $vqcol, $infocol, $formatcol);
	my ($scaffold, $pos, $prevscaff, $prevpos, $locindex);
	my $switch=0;

	open my $vcf, "-|", $command or die "Could not open raw vcf file!\n$!\n";

	READ: while (<$vcf>){
		chomp;

		if(/^#CHROM/){
			@header=split("\t", $_);
			COL: for my $col (0..$#header){
				if ($header[$col] eq "POS"){$poscol=$col}
				elsif ($header[$col] eq "ALT"){$altcol=$col}
				elsif ($header[$col] eq "QUAL"){$vqcol=$col}
				elsif ($header[$col] eq "INFO"){$infocol=$col}
				elsif ($header[$col] eq "FORMAT"){$formatcol=$col; last COL}
				}
#			print STDERR "$poscol, $altcol, $vqcol, $infocol, $formatcol\n";
			@samples=@header[($formatcol+1)..$#header];
			print STDERR join(',', @samples), "\n";
			$switch=1;
			}

		elsif(!$switch){next READ}
		elsif(/^INFO/){next READ}

		elsif(/^\w+\s/){
			@line=split("\t", $_);
			$scaffold=$line[0];
			$pos=$line[$poscol]-1; # The sam and vcf formats are 1-based. All internal arrays and hashes use a 0-based indexing.
			unless (defined $prevpos && $scaffold eq $prevscaff && $pos-$prevpos<10000){
				$locindex=$ref_loci->findLocus($scaffold, $pos);
				$offsets{$scaffold}[$locindex]=$ref_loci->getValue($scaffold, $locindex, 'start');
				}

			$nsamples=$#line-$formatcol;
			if ($nsamples!=$totsize){die "Number or samples in raw vcf file does not match information on command line!\n"}

			if ($line[$altcol] eq '.'){$snp=0} else {$snp=1}

			$vq=$line[$vqcol];
			unless ($snp || ($vq ne '.' && $vq>=$minvq) ){
				for my $allind (0..$totsize-1){$missing{$scaffold}[$locindex][$allind].=chr(0)}
				next READ;
				}

			($mq)=($line[$infocol]=~/MQ=(\d+)/);
			unless (defined $mq && $mq>=$minmq){
				for my $allind (0..$totsize-1){$missing{$scaffold}[$locindex][$allind].=chr(0)}
				next READ;
				}

			if ($snp){$dp=2} else {$dp=1}
#			@label=split(":", $line[$formatcol]);
#			LABEL: for my $lab (0..$#label){if ($label[$lab] eq "DP"){$dp=$lab; last LABEL}}

			if (defined $prevpos && $scaffold eq $prevscaff && $pos-$prevpos!=1 && $pos-$prevpos<10000){
#				warn "Gap inbetween positions $prevscaff: $prevpos and $scaffold: $pos!\n";
				my $gapsize=$pos-$prevpos-1;
				for my $allind (0..$totsize-1){$missing{$scaffold}[$locindex][$allind].=(chr(0) x $gapsize)}
				}
			$prevscaff=$scaffold;
			$prevpos=$pos;

			@individuals=@line[-$nsamples..-1];
			for my $allind (@$sampleorder){
#			for my $allind (0..$totsize-1){
				my @data=split(":", (shift @individuals));
				if (@data>1 && $data[$dp]<255){$missing{$scaffold}[$locindex][$allind].=chr($data[$dp])}
				elsif (@data>1){$missing{$scaffold}[$locindex][$allind].=chr(255)}
				else {$missing{$scaffold}[$locindex][$allind].=chr(0)}
				}
			}
		}

	close $vcf;
	warn "\$?=$?, \$!=$!\n";

	$self->{_Missing}=\%missing;
	$self->{_Header}=\@header;
	$self->{_Pops}=$popsize;
	$self->{_Offsets}=\%offsets;
	return;
	}


sub readAllsites{	
	my ($self, $ref_loci, $minmq, $minvq, $ref_popsizes, $excludeRefN, $verbose)=@_;

	unless (defined $self->{_Files} && defined $self->{_Samples} ){ die "List and folder of allsites vcf files have not been defined!\n" }
	my @popsizes=@{ $self->{_Pops} };
#	my $totsize=Misc::sum(@popsizes);

	if (defined $ref_popsizes){	# allows for population size check.
		if (@$ref_popsizes!=@popsizes){ warn "Wrong number of populations specified in allsiteslist (", scalar(@$ref_popsizes), " vs. ", scalar(@popsizes), ")!\n" }
		for my $pop (0..$#{ $ref_popsizes } ){
			if ( $ref_popsizes->[$pop]!=$popsizes[$pop] ){ warn "Wrong number of samples specified in allsiteslist for population $pop ($ref_popsizes->[$pop] vs. ", scalar($popsizes[$pop]), ")!\n" }
			}
		}

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my $reference=$ref_loci->getValue($chrom, $locus, 'refseq');
			unless (defined $reference){ warn "Reference sequence for $chrom:$start-$end is not defined!\n" if $verbose; $reference="" }
			my @missing;
			my $indsum=0;
			for my $pop (0..$#{ $self->{_Files} }){
				open my $vcffile, "-|", "tabix $self->{_Files}[$pop] $chrom:$start-$end" or die "Could not open vcf file $self->{_Files}[$pop]!\n$!\n";
				print STDERR "Acquiring coverage data for $chrom:$start-$end, population $pop.\n" if $verbose;
				my $refseq=$reference;	# make hard copy of reference sequence.
				my $prevpos=$start-1;
				my ($scaffold, $pos, $mq, $vq, $dp);
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1];

						if ($pos-$prevpos!=1){
							my $gapsize=$pos-$prevpos-1;
							for my $ind (0..$popsizes[$pop]-1){ $missing[$indsum+$ind].=(chr(0) x $gapsize) }
							substr($refseq, 0, $gapsize, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
							}
						$prevpos=$pos;

						my $refbase=uc( substr($refseq, 0, 1, "") );
						if ($excludeRefN){
							if ($refbase eq 'N' ){	# discard if base in reference sequence is hardmasked.
								for my $ind (0..$popsizes[$pop]-1){ $missing[$indsum+$ind].=chr(0) }
								next LINE;
								}
							}
						$vq=$line[5];
						unless ($vq ne '.' && $vq>=$minvq){
							for my $ind (0..$popsizes[$pop]-1){ $missing[$indsum+$ind].=chr(0) }
							next LINE;
							}
						($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
						unless (defined $mq && $mq>=$minmq){
							for my $ind (0..$popsizes[$pop]-1){ $missing[$indsum+$ind].=chr(0) }
							next LINE
							}

						$dp=$line[4] eq '.' ? 1 : 2;

						my $allind=$indsum;
						for my $ind ( @{ $self->{_Samples}[$pop] } ){
							my @data=split(":", $line[$ind+8]);
							if (@data>1 && $data[$dp]<255){ $missing[$allind].=chr($data[$dp]) }
							elsif (@data>1){ $missing[$allind].=chr(255) }
							else { $missing[$allind].=chr(0) }
							++$allind;
							}
						}
					}
				close $vcffile;
				my $gapsize=defined $pos ? $end-$pos: $locuslength;
				for my $ind (0..$popsizes[$pop]-1){
					if ($gapsize){ $missing[$indsum+$ind].=(chr(0) x $gapsize) }
					my $stringlength=length( $missing[$indsum+$ind] );
					if ($stringlength != $locuslength){
						warn "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				$indsum+=$popsizes[$pop];
				}
			$self->{_Missing}{$chrom}{$start-1}=\@missing;			
			}
		}
	return;
	}


sub readAllsitesBlocks{
	my ($self, $ref_loci, $allsitesfolder, $allsiteslist, $minmq, $minvq, $npops, $popsizes, $excludeRefN)=@_;

#	my $totsize=Misc::sum(@$popsizes);
	my ($files, $samples);
	if ( defined $self->{_Files} && defined $self->{_Samples} ){
		$files=$self->{_Files};
		$samples=$self->{_Samples};
		}
	elsif (defined $allsitesfolder && defined $allsiteslist){
		($files, $samples)=Misc::openvcflist($allsitesfolder, $allsiteslist);
		}
	else { die "List of allsites vcf files have not been provided or set before calling readAllistes!\n" }

	unless (@{$files}==$npops){ die "Wrong number of allsites files specified!\n" }
	for my $i (0..$npops-1){
		if ( $popsizes->[$i]!=@{$samples->[$i]} ){ die "Wrong number of samples specified in allsiteslist for population $i ($popsizes->[$i] vs. ", scalar(@{$samples->[$i]}), ")!\n" }
		}

	my @indnames;
	for my $pop (0..$#{$files}){
		open my $vcfheader, "-|", "tabix $files->[$pop] -H" or die "Could not open vcf header $files->[$pop]!\n$!\n";
		while (<$vcfheader>){
			if(/^#CHROM/){
				chomp;
				my @header=split("\t", $_);
				for my $sample ( @{$samples->[$pop]} ){ push @indnames, $header[8+$sample] }
				}
			}
		close $vcfheader;
		}

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my @missstring;
			my $indsum=0;
			for my $pop (0..$#{$files}){
				open my $vcffile, "-|", "tabix $files->[$pop] $chrom:$start-$end" or die "Could not open vcf file $files->[$pop]!\n$!\n";
				print STDERR "Acquiring coverage data for $chrom:$start-$end, population $pop.\n";
				my $refseq=$ref_loci->getValue($chrom, $locus, 'refseq');
				unless (defined $refseq){ warn "Reference sequence for $chrom:$start-$end is not defined!\n"; $refseq="" }
				my $prevpos=$start-1;
				my ($scaffold, $pos, $mq, $vq, $dp);
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1];

						if ($pos-$prevpos!=1){
							my $gapsize=$pos-$prevpos-1;
							for my $ind (0..$popsizes->[$pop]-1){ $missstring[$indsum+$ind].=(chr(0) x $gapsize) }
							substr($refseq, 0, $gapsize, "");	# chew up string to keep current position at the beginning of the string to speed up substr extraction.
							}
						$prevpos=$pos;

						my $refbase=uc( substr($refseq, 0, 1, "") );
						if ($excludeRefN){
							if ($refbase eq 'N' ){	# discard if base in reference sequence is hardmasked.
								for my $ind (0..$popsizes->[$pop]-1){ $missstring[$indsum+$ind].=chr(0) }
								next LINE;
								}
							}
						$vq=$line[5];
						unless ($vq ne '.' && $vq>=$minvq){
							for my $ind (0..$popsizes->[$pop]-1){ $missstring[$indsum+$ind].=chr(0) }
							next LINE;
							}
						($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
						unless (defined $mq && $mq>=$minmq){
							for my $ind (0..$popsizes->[$pop]-1){ $missstring[$indsum+$ind].=chr(0) }
							next LINE
							}

						$dp=$line[4] eq '.' ? 1 : 2;

						my $allind=$indsum;
						for my $ind (@{$samples->[$pop]}){
							my @data=split(":", $line[$ind+8]);
							if (@data>1 && $data[$dp]<255){ $missstring[$allind].=chr($data[$dp])}
							elsif (@data>1){ $missstring[$allind].=chr(255)}
							else { $missstring[$allind].=chr(0)}
							++$allind;
							}
						}
					}
				close $vcffile;
				my $gapsize=defined $pos ? $end-$pos: $locuslength;
				for my $ind (0..$popsizes->[$pop]-1){
					if ($gapsize){ $missstring[$indsum+$ind].=(chr(0) x $gapsize) }
					my $stringlength=length( $missstring[$indsum+$ind] );
					if ($stringlength != $locuslength){
						warn "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				$indsum+=$popsizes->[$pop];
				}

			my $offset=($start-1) % 100;
			my $startpos=($start-1)-$offset;
			my $overhang=($offset>0) ? 100-$offset : 0;
#			print STDERR "$chrom: $start-$end, $offset, $startpos, $overhang\n";
			unless ( $self->{_Missing}{$chrom}{$startpos} ){	# there is already data in the block.
				@{ $self->{_Missing}{$chrom}{$startpos} }=( (chr(0) x 100) x scalar(@missstring) );
				}
			if ($locuslength>$overhang){
				if ($overhang>0){
					substr( $self->{_Missing}{$chrom}{$startpos}[$_], -$overhang, $overhang)=substr($missstring[$_], 0, $overhang, "") foreach (0..$#missstring);
					$startpos+=100;
					}
				for my $block ( 1..($locuslength-$overhang)/100 ){
					$self->{_Missing}{$chrom}{$startpos}[$_]=substr($missstring[$_], 0, 100, "") foreach (0..$#missstring);
					$startpos+=100;
					}
				my $reslength=length( $missstring[0] );
#				print STDERR "residual string length: $reslength, starting position: $startpos\n";
				if ($reslength>0){
					@{ $self->{_Missing}{$chrom}{$startpos} }=( (chr(0) x 100) x scalar(@missstring) ) unless ( $self->{_Missing}{$chrom}{$startpos} );
					substr( $self->{_Missing}{$chrom}{$startpos}[$_], 0, $reslength)=substr($missstring[$_], 0, $reslength, "") foreach (0..$#missstring);
					}
				}
			else { substr( $self->{_Missing}{$chrom}{$startpos}[$_], -$overhang, $locuslength)=substr($missstring[$_], 0, $locuslength, "") foreach (0..$#missstring) }
			for my $indstring (@missstring){ die "Wrong residual string length!\n" unless length($indstring)==0 }
			}
		}

	$self->{_Pops}=$popsizes;
	$self->{_Indnames}=\@indnames;
	return;
	}


sub readAllsitesGeno{
	my ($self, $ref_loci, $allsitesfolder, $allsiteslist, $npops, $popsize, $minmq, $minvq, $mingtq, $suballeles, $recordbases, $recordmq)=@_;

	my ($files, $samples)=Misc::openvcflist($allsitesfolder, $allsiteslist);
	unless (@{$files}==$npops){die "Wrong number of allsites files specified!\n"}
	for my $i (0..$npops-1){
		if ( $popsize->[$i]!=@{$samples->[$i]} ){die "Wrong number of samples specified in allsiteslist for population $i!\n"}
		}

	my $chromosome=$ref_loci->nextIterator();
	my @indnames;
	for my $pop (0..$#{$files}){
		open my $vcfheader, "-|", "tabix $files->[$pop] -h $chromosome:1-1" or die "Could not open vcf header $files->[$pop]!\n$!\n";
		while (<$vcfheader>){
			if(/^#CHROM/){
				chomp;
				my @header=split("\t", $_);
				for my $sample ( @{$samples->[$pop]} ){ push @indnames, $header[8+$sample] }
				}
			}
		close $vcfheader;
		}

	my %missing;	# $missing{chrom/scaff}[locus][allind]="char string";
	my %mq;	# $mq{chrom/scaff}[locus]="char string";
	my $genotypes=Genotypes->new();

	CHROM: for my $chrom ( $ref_loci->allKeys() ){
		LOCUS: for my $locus ( $ref_loci->allKeys($chrom) ){
			my $start=$ref_loci->getValue($chrom, $locus, 'start')+1;	# vcf files are 1-based, internal storage is 0-based.
			my $end=$ref_loci->getValue($chrom, $locus, 'end')+1;
			my $locuslength=$end-$start+1;
			my $indsum=0;
			my $iterator;
			if ($suballeles){ $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus) }	# but data is not yet subsampled!
			for my $pop (0..$#{$files}){
				open my $vcffile, "-|", "tabix $files->[$pop] $chrom:$start-$end" or die "Could not open vcf file $files->[$pop]!\n$!\n";
				print STDERR "Aquiring coverage and genotype data for $chrom:$start-$end, population $pop.\n";
				my $prevpos=$start-2;
				my ($pos, $scaffold, $snp, $dp);
				LINE: while (<$vcffile>){
					if(/^\w+\s/){
						my @line=split("\t", $_);
						$scaffold=$line[0];
						$pos=$line[1]-1;	# vcf files are 1-based, internal storage is 0-based.
						my $vq=$line[5];

						my $gapsize=$pos-$prevpos-1;
						if ($gapsize>0){	# Attention: negative values evalute to true!
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=(chr(0) x $gapsize) }
							$mq{$scaffold}[$locus].=(chr(0) x $gapsize) if ($recordmq);
							}
						$prevpos=$pos;

						my ($mq)=($line[7]=~/MQ=(\d+\.?\d*)/);
						if ($recordmq){
							if (defined $mq && $mq<255){ $mq{$scaffold}[$locus].=chr( int($mq) ) }
							elsif (defined $mq){ $mq{$scaffold}[$locus].=chr(255) }
							else { $mq{$scaffold}[$locus].=chr(0) }
							}
						unless ($vq ne '.' && $vq>=$minvq){
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=chr(0) }
							next LINE;
							}
						unless (defined $mq && $mq>=$minmq){
							for my $ind (0..$popsize->[$pop]-1){ $missing{$scaffold}[$locus][$indsum+$ind].=chr(0) }
							next LINE
							}

						if ($line[4] eq '.'){
							$snp=0; $dp=1;
							}
						else {
							$snp=1; $dp=2;
							$genotypes->{_Alleles}{$scaffold}{$pos}=0 unless defined $genotypes->{_Alleles}{$scaffold}{$pos};
							unless ( $pop && exists $genotypes->{_Genotypes}{$scaffold}{$pos} ){
								for my $i (0..$npops-1){ $genotypes->{_Genotypes}{$scaffold}{$pos}[$i]='00' x $popsize->[$i] }
								}
							$genotypes->{_Genotypes}{$scaffold}{$pos}[$pop]='';
							if ($recordbases){ $genotypes->{_Bases}{$scaffold}{$pos}="$line[3]" . join('', split(',' , $line[4]) ) }
							}

						my $allind=$indsum;
						my $index=0;
						for my $ind (@{$samples->[$pop]}){
							my @data=split(":", $line[$ind+8]);
							if (@data>1 && $data[$dp]<255){ $missing{$scaffold}[$locus][$allind].=chr($data[$dp]) }
							elsif (@data>1){ $missing{$scaffold}[$locus][$allind].=chr(255) }
							else { $missing{$scaffold}[$locus][$allind].=chr(0) }
							if ($snp){
								if ($data[0]=~/([0-3])\/([0-3])/){
									if ( $data[3]>=$mingtq){ $genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=$1 . $2 }
									else { $genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=".." }
									if ($suballeles && defined $iterator->[$pop][$index]){	# only records allele states for subsampled individuals.
										if (2*$ind==$iterator->[$pop][$index]){$genotypes->setAlleles($scaffold, $pos, $1); ++$index}
										if (2*$ind+1==$iterator->[$pop][$index]){$genotypes->setAlleles($scaffold, $pos, $2); ++$index}
										}
									else { $genotypes->setAlleles($scaffold, $pos, $1, $2) }
									}
								else {$genotypes->{_Genotypes}{$scaffold}{$pos}[$pop].=".."}
								}
							++$allind;
							}
						}
					}
				close $vcffile;
				my $gapsize=defined $pos ? $end-$pos-1 : $locuslength;
				if ($recordmq && $gapsize){ $mq{$scaffold}[$locus].=(chr(0) x $gapsize) }
				for my $ind (0..$popsize->[$pop]-1){
					if ($gapsize){ $missing{$chrom}[$locus][$indsum+$ind].=(chr(0) x $gapsize) }
					my $stringlength=length($missing{$chrom}[$locus][$indsum+$ind]);
					if ($stringlength != $locuslength){
						warn "Non-matching string length for $chrom:$start-$end! $stringlength vs. $locuslength\n";
						}
					}
				$indsum+=$popsize->[$pop];
				}
			}
		}
	$self->{_Missing}=\%missing;
	$self->{_Pops}=$popsize;
	$self->{_Indnames}=\@indnames;
	$self->{_MQ}=\%mq;
	return($genotypes);
	}


sub missingReport_byposition{ # currently not valid!
	my ($self, $listfile, $folder, $file_prefix, $ref_loci, $mincov)=@_;
	my $l=0;

	open my $filelist, ">", $listfile or die "Could not open filelist $listfile!\n";
	mkdir $folder;

	CHROM: for my $chrom (sort $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if $ref_loci->getExcluded($chrom, $locus)==1;
			my $outfile=$folder . "/" . $file_prefix . "_" . $l . ".txt";
			print $filelist "$outfile\n";
			open my $output, ">", $outfile or die "Could not open report file $outfile!\n";

			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $npops=$self->getNPop();
			my $n=1;
			my @missing;
			my @rangestart;
			my @rangeend;

			for my $pos ($start..$end){
				my $allind=0;
				for my $pop (0..$npops-1){
					my $iterator=$ref_loci->getreftoIndividuals('position', $chrom, $pos);
					if (!defined $iterator){die "Array reference to chromosome $chrom, position $pos is not defined!\n"}
					for my $ind (@{$iterator->[$pop]}){
						unless (defined $self->{_Missing}{$chrom}{$pos}[$pop][$ind] && $self->{_Missing}{$chrom}{$pos}[$pop][$ind]>=$mincov){
							if (!defined $rangeend[$allind]){$rangestart[$allind]=$n; $rangeend[$allind]=$n}
							elsif (($n-1)==$rangeend[$allind]){$rangeend[$allind]=$n}
							else {
								push @{$missing[$allind][0]}, ($rangestart[$allind]-1)/$locuslength;
								push @{$missing[$allind][1]}, $rangeend[$allind]/$locuslength;
								$rangestart[$allind]=$n; $rangeend[$allind]=$n;
								}
							}
						if ($n==$locuslength && defined $rangeend[$allind]){
							push @{$missing[$allind][0]}, ($rangestart[$allind]-1)/$locuslength;
							push @{$missing[$allind][1]}, $rangeend[$allind]/$locuslength;
							}
						$allind++;
						}
					}
				$n++;
				}
			for my $allind (0..$#missing){
				for my $range (0..$#{$missing[$allind][0]}){
					print $output "$allind\t$missing[$allind][0][$range]\t$missing[$allind][1][$range]\n";
					}
				}
			print STDERR "$chrom\t$start\t$end\n";
			close $output;
			$l++;
			}
		}
	close $filelist;
	return;
	}


sub missingReport_bylocus{
	my ($self, $listfile, $folder, $file_prefix, $ref_loci, $mincov)=@_;
	my $npops=$self->getNPop();
	my $l=0;

	mkdir $folder;
	open my $filelist, ">", $listfile or die "Could not open filelist $listfile!\n";
	open my $locuslist, ">", "$folder/../locuslist.txt" or die "Could not open locuslist $listfile!\n";

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if $ref_loci->getExcluded($chrom, $locus)==1;
			my $outfile=$folder . "/" . $file_prefix . "_" . $l . ".txt";
			print $filelist "$outfile\n";
			open my $output, ">", $outfile or die "Could not open report file $outfile!\n";

			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
			if (!defined $iterator){die "Array reference to chromosome $chrom, locus $locus: $start-$end is not defined!\n"}
			my (@missing, @nvalid);
			my $rt=0;
			my $lineage=0;

			for my $pop (0..$npops-1){
				for my $index (@{$iterator->[$pop]}){
					my $allind=int($index/2)+$rt;
					my $mstring=$self->getValue($chrom, $start, $allind);
					unless (length($mstring)==$locuslength){ die "Wrong length of missing string!\n" }
					my $rangestart;
					my $rangeend;
					for my $n (1..$locuslength){
						my $coverage=ord( substr($mstring, 0, 1, "") );
						if (defined $coverage && $coverage>=$mincov){ ++$nvalid[$n-1] }
						else {
							if (!defined $rangeend){$rangestart=$n; $rangeend=$n}
							elsif (($n-1)==$rangeend){$rangeend=$n}
							else {
#								push @{$missing[$lineage]}, ($rangestart-1)/$locuslength . "\t" . $rangeend/$locuslength;
								print $output "$lineage\t", ($rangestart-1)/$locuslength, "\t", $rangeend/$locuslength, "\n";
								$rangestart=$n; $rangeend=$n;
								}
							}
						}
					if (defined $rangeend){
#						push @{$missing[$lineage]}, ($rangestart-1)/$locuslength . "\t" . $rangeend/$locuslength;
						print $output "$lineage\t", ($rangestart-1)/$locuslength, "\t", $rangeend/$locuslength, "\n";
						}
					++$lineage;
					}
				$rt+=$self->{_Pops}[$pop];
				}
#			for my $allind (0..$#missing){
#				for my $range (0..$#{$missing[$allind]}){
#					print $output "$allind\t$missing[$allind][$range]\n";
#					}
#				}
			print $locuslist "Locus $l\t$chrom - $locus: $start-$end\n";
			close $output;
			++$l;
			}
		}
	close $filelist;
	close $locuslist;
	return;
	}

sub missingReport_bylocusBlocks{
	my ($self, $listfile, $folder, $file_prefix, $ref_loci, $mincov)=@_;
	my $npops=$self->getNPop();
	my $l=0;

	mkdir $folder;
	open my $filelist, ">", $listfile or die "Could not open filelist $listfile!\n";
	open my $locuslist, ">", "$folder/../locuslist.txt" or die "Could not open locuslist $listfile!\n";

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if $ref_loci->getExcluded($chrom, $locus)==1;
			my $outfile=$folder . "/" . $file_prefix . "_" . $l . ".txt";
			print $filelist "$outfile\n";
			open my $output, ">", $outfile or die "Could not open report file $outfile!\n";

			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
			if (!defined $iterator){die "Array reference to chromosome $chrom, locus $locus: $start-$end is not defined!\n"}
			my (@missing, @nvalid);
			my $rt=0;
			my $lineage=0;

			for my $pop (0..$npops-1){
				for my $index (@{$iterator->[$pop]}){
					my $allind=int($index/2)+$rt;
					my $rangestart;
					my $rangeend;
					for my $n (1..$locuslength){
						my $coverage=$self->getValuebyPos($chrom, $start+$n-1, $allind);
						if (defined $coverage && $coverage>=$mincov){ ++$nvalid[$n-1] }
						else {
							if (!defined $rangeend){$rangestart=$n; $rangeend=$n}
							elsif (($n-1)==$rangeend){$rangeend=$n}
							else {
								print $output "$lineage\t", ($rangestart-1)/$locuslength, "\t", $rangeend/$locuslength, "\n";
								$rangestart=$n; $rangeend=$n;
								}
							}
						}
					if (defined $rangeend){
						print $output "$lineage\t", ($rangestart-1)/$locuslength, "\t", $rangeend/$locuslength, "\n";
						}
					++$lineage;
					}
				$rt+=$self->{_Pops}[$pop];
				}
			print $locuslist "Locus $l\t$chrom - $locus: $start-$end\n";
			close $output;
			++$l;
			}
		}
	close $filelist;
	close $locuslist;
	return;
	}


sub validSites{
	my ($self, $outfile, $ref_loci, $ref_groups, $ref_minind, $mincov)=@_;

	my $npops=$self->getNPop();
	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		$rt+=$self->getNInd($pop);
		}

	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";
	my @nvalid;	# $nvalid[group][totlocus]
	my $totlocus=0;

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if $ref_loci->getExcluded($chrom, $locus)==1;
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
			if (!defined $iterator){die "Array reference to chromosome $chrom, locus $locus: $start-$end is not defined!\n"}
			GROUP: for my $group (0..$#{$ref_groups}){
				my $validpos=0;
				POSITION: for my $offset (0..$locuslength-1){
					for my $pop (@{$ref_groups->[$group]}){
						my $validind=0;
						for my $index (@{$iterator->[$pop]}){
							my $allind=int($index/2)+$popsum[$pop];
							my $coverage=$self->getValue($chrom, $start, $allind, $offset);
							++$validind if (defined $coverage && $coverage>=$mincov);
							}
						next POSITION if $validind<$ref_minind->[$pop];
						}
					++$validpos;
					}
				$nvalid[$group][$totlocus]=$validpos;
				}
			++$totlocus;
			}
		}
	for my $group (0..$#nvalid){
		print $output "$group\t", join(',', @{$ref_groups->[$group]}), "\t", join(',', @{$nvalid[$group]}), "\n";
		}
	close $output;
	return;
	}


sub validSitesBlocks{
	my ($self, $outfile, $ref_loci, $ref_groups, $ref_minind, $mincov)=@_;

	my $npops=$self->getNPop();
	my @popsum;
	my $rt=0;
	for my $pop ( 0..($npops-1) ){
		push @popsum, $rt;
		$rt+=$self->getNInd($pop);
		}

	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";
	my @nvalid;	# $nvalid[group][totlocus]
	my $totlocus=0;

	CHROM: for my $chrom (sort {Misc::expand($a) cmp Misc::expand($b)} $ref_loci->allKeys()){
		LOCUS: for my $locus ($ref_loci->allKeys($chrom)){
			next LOCUS if $ref_loci->getExcluded($chrom, $locus)==1;
			my $start=$ref_loci->getValue($chrom, $locus, 'start');
			my $end=$ref_loci->getValue($chrom, $locus, 'end');
			my $locuslength=$end-$start+1;
			my $iterator=$ref_loci->getreftoIndividuals('locus', $chrom, $locus);
			if (!defined $iterator){die "Array reference to chromosome $chrom, locus $locus: $start-$end is not defined!\n"}
			GROUP: for my $group (0..$#{$ref_groups}){
				my $validpos=0;
				POSITION: for my $pos ($start..$end){
					for my $pop (@{$ref_groups->[$group]}){
						my $validind=0;
						for my $index (@{$iterator->[$pop]}){
							my $allind=int($index/2)+$popsum[$pop];
							my $coverage=$self->getValuebyPos($chrom, $pos, $allind);
							++$validind if (defined $coverage && $coverage>=$mincov);
							}
						next POSITION if $validind<$ref_minind->[$pop];
						}
					++$validpos;
					}
				$nvalid[$group][$totlocus]=$validpos;
				}
			++$totlocus;
			}
		}
	for my $group (0..$#nvalid){
		print $output "$group\t", join(',', @{$ref_groups->[$group]}), "\t", join(',', @{$nvalid[$group]}), "\n";
		}
	close $output;
	return;
	}


1;

