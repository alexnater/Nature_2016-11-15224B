#!/usr/bin/env perl

use strict;
use warnings;
use Dumpvalue;

# written by Alexander Nater, November 2013

my $task=shift @ARGV;

if ($task eq "Filter" && scalar (@ARGV) != 6){die "Arguments missing!\nusage: [Task] [Name of vcf file] [Name of output file] [MQRSFilter threshold] [RPRSFilter threshold] [Coverage filter threshold] [distance between clusters]"}
elsif ($task eq "Filter"){
	my ($inname, $outname, $MQRSthres, $RPRSthres, $dpthres, $clusterdist)=@ARGV;
	my $depthlimit=DPLimit($inname, $dpthres);
	FilterVcf($inname, $outname, $MQRSthres, $RPRSthres, $depthlimit, $clusterdist);
	}

elsif ($task eq "CRISP" && scalar (@ARGV) != 6){die "Arguments missing!\nusage: [Task] [Name of vcf file] [Name of output file] [minimal variant quality] [minimal coverage] [maximal coverage] [min Quality-by-depth filter]"}
elsif ($task eq "CRISP"){
	my ($inname, $outname, $minqual, $mindp, $maxdp, $minqd)=@ARGV;
#	my $maxdp=DPLimit($inname, $dpthres);
	FilterCRISP($inname, $outname, $minqual, $mindp, $maxdp, $minqd);
	}

elsif ($task eq "Count" && scalar (@ARGV) != 3){die "Arguments missing!\nusage: [Task] [Name of vcf file] [Number of samples] [minimal sample coverage]"}
elsif ($task eq "Count"){
	my ($inname, $totsize, $minimalcov)=@ARGV;
	countgenotypes($inname, $totsize, $minimalcov);
	}

elsif ($task eq "DPStats" && scalar (@ARGV) != 2){die "Arguments missing!\nusage: [Task] [Name of doc file] [thinning interval]"}
elsif ($task eq "DPStats"){
	my ($inname, $thinint)=@ARGV;
#	my $dplimit = DPLimitDoC($inname, $thinint, 5);
	my ($mean, $sd, $n) = meansd_stream($inname, $thinint);
	print "Number of positions considered: $n\nMean: $mean\nStandard deviation: $sd\n";
	}

elsif ($task eq "VQSR" && scalar (@ARGV) != 3){die "Arguments missing!\nusage: [Task] [Name of vcf file] [proportion of top variants for variant quality threshold]"}
elsif ($task eq "VQSR"){
	my ($inname, $outname, $qualthres)=@ARGV;
	my $quallimit=QualLimit($inname, $qualthres);
	TopVariants($inname, $outname, $quallimit);
	}

elsif ($task eq "HardFilter" && scalar (@ARGV) != 8){die "Arguments missing!\nusage: [Task] [Name of vcf file] [Name of output file] [Number of samples] [minimal mapping quality] [minimal distance between SNPs] [minimal sample coverage] [minimal number of alternative reads] [minimal genotype quality]"}
elsif ($task eq "HardFilter"){
	my ($inname, $outname, $totsize, $minmq, $mindist, $mindp, $minaltreads, $minaltgtqual)=@ARGV;
	HardfilterVcf($inname, $outname, $totsize, $minmq, $mindist, $mindp, $minaltreads, $minaltgtqual);
	}

else {die "Task $task not known!"}


sub FilterCRISP{
	my ($filename, $outfile, $minqual, $mindp, $maxdp, $minqd)=@_;
	my $totalcount=0;
	my $interval=100000;
	my @count=(0, 0, 0, 0, 0, 0);

	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";
	open my $input, "-|", "tabix -H $filename" or die "Could not open vcf file $filename!\n";	
	while (<$input>){ print $output $_ if(/^#/) }
	close $input;
	open $input, "-|", "tabix $filename ." or die "Could not open vcf file $filename!\n";	

	LINE: while (<$input>){
		my @call=split(/\s/, $_);
		unless ($call[7]=~/VT=SNV;/){ ++$count[0]; next }
		my %filter = map { $_=>1 } split(/;/, $call[6]);
		if (exists $filter{'StrandBias'} || exists $filter{'LowMQ20'} || exists $filter{'LowMQ10'}){ ++$count[1]; next }
		if ($call[5]<$minqual){
			++$count[2];
			next;
			}
		if ($call[7]=~/DP=(\d+),(\d+),(\d+)/){
			my $dp = $1+$2+$3;
			if ($dp<$mindp){
				++$count[3];
				next;
				}
			if ($dp>$maxdp){
				++$count[4];
				next;
				}
			if ($dp && $call[5]/$dp<$minqd){
				++$count[5];
				next;
				}
			}
		else { print "Could not find coverage information for variant $call[0]:$call[1]!\n" }
		print $output $_;
		} continue { ++$totalcount; unless ($totalcount % $interval){ print "$totalcount variants processed.\n" } }

	close $input;
	close $output;

	print "Total number of variants: $totalcount.\n";
	print "$count[0] (", sprintf("%.2f", $count[0]/$totalcount*100), "%) variants removed due to indel.\n";
	print "$count[1] (", sprintf("%.2f", $count[1]/$totalcount*100), "%) removed due to filter tag.\n";
	print "Min QUAL filter applied to $count[2] (", sprintf("%.2f", $count[2]/$totalcount*100), "%) variants.\n";
	print "Min DP filter applied to $count[3] (", sprintf("%.2f", $count[3]/$totalcount*100), "%) variants.\n";
	print "Max DP filter applied to $count[4] (", sprintf("%.2f", $count[4]/$totalcount*100), "%) variants.\n";
	print "Min QD filter applied to $count[5] (", sprintf("%.2f", $count[5]/$totalcount*100), "%) variants.\n";
	return;
	}


sub FilterVcf{
	my ($filename, $outfile, $MQRS, $RPRS, $dplimit, $mindist)=@_;
	my @header;
	my @call;
	my $iscluster=0;
	my $nclusters=0;
	my $varincluster=0;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	my @depth;
	my @count=(0, 0, 0);
	my @prevcall;
	my $modline;
	my $filterline=0;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";
	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";

	LINE: while (<$vcffile>){
		if(/^##FILTER/){
			$filterline=1;
			print $output $_;
			}
		elsif(/^#CHROM/){
			@header=split(/\s/, $_);
			print $output $_;
			}
		elsif(!/^#/ && /^\w+\s/){
			@call=split(/\s/, $_);
			if ($call[6] eq '.'){$call[6]="PASS"}
			if (defined $prevcall[0]){
				my $clustertest=0;
				$clustertest=1 if ($call[0] eq $prevcall[0] && ($call[1]-$prevcall[1])<$mindist);
				if ($clustertest){
					unless ($iscluster){
						$nclusters++;
						$varincluster++;
						}
					$varincluster++;
					$iscluster=1;
					}
				if ($iscluster){
					if ($prevcall[6] eq "PASS"){$prevcall[6]="ClusterFilter"}
					else {$prevcall[6]=join(';', $prevcall[6], "ClusterFilter")}
					}
				unless ($clustertest){$iscluster=0}
				$modline=join("\t", @prevcall);
				print $output "$modline\n";
				}

			if ($call[7]=~/MQRankSum=(-?\d+\.?\d*)/){
				if ($1<$MQRS){
					if ($call[6] eq "PASS"){$call[6]="MQRSFilter"}
					else {$call[6]=join(';', $call[6], "MQRSFilter")}
					$count[0]++;
					}
				}
#			else {print "Could not find MQRankSum information for variant $call[0]:$call[1]!\n"}

			if ($call[7]=~/ReadPosRankSum=(-?\d+\.?\d*)/){
				if ($1<$RPRS){
					if ($call[6] eq "PASS"){$call[6]="RPRSFilter"}
					else {$call[6]=join(';', $call[6], "RPRSFilter")}
					$count[1]++;
					}
				}
#			else {print "Could not find ReadPosRankSum information for variant $call[0]:$call[1]!\n"}

			if ($call[7]=~/DP=(\d+)/){
				if ($1>$dplimit){
					if ($call[6] eq "PASS"){$call[6]="DPFilter"}
					else {$call[6]=join(';', $call[6], "DPFilter")}
					$count[2]++;
					}
				}
			else {print "Could not find coverage information for variant $call[0]:$call[1]!\n"}

			++$totalcount;
			--$counter;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			@prevcall=@call;
			}
		else {
			if ($filterline){
				print $output '##FILTER=<ID=MQRSFilter,Description="MQRankSum < ', $MQRS,'">', "\n";
				print $output '##FILTER=<ID=RPRSFilter,Description="ReadPosRankSum < ', $RPRS,'">', "\n";
				print $output '##FILTER=<ID=DPFilter,Description="DP > ', $dplimit,'">', "\n";
				print $output '##FILTER=<ID=ClusterFilter,Description="Located within ', $mindist, ' bp of other SNP">', "\n";
				$filterline=0;
				}
			print $output $_;
			}
		}

	$modline=join("\t", @prevcall);
	print $output "$modline\n";
	close $vcffile;
	close $output;

	print "Total number of variants: $totalcount.\n";
	print "MQRSFilter applied to $count[0] (", $count[0]/$totalcount*100, "%) variants.\n";
	print "RPRSFilter applied to $count[1] (", $count[1]/$totalcount*100, "%) variants.\n";
	print "DPFilter applied to $count[2] (", $count[2]/$totalcount*100, "%) variants.\n";
	print "Number of variants in clusters: $varincluster, in a total of $nclusters clusters.\n";
	return;
	}


sub FilterVcf_old{
	my ($filename, $outfile, $MQRS, $RPRS, $dplimit)=@_;
	my @header;
	my @call;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	my @depth;
	my @count=(0, 0, 0);
	my $modline;
	my $filterline=0;
	
        open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";
        open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";

	LINE: while (<$vcffile>){
		if(/^##FILTER/){
			$filterline=1;
			print $output $_;
			}
		elsif(/^#CHROM/){
			@header=split(/\s/, $_);
			print $output $_;
			}
		elsif(!/^#/ && /^\w+\s/){
			chomp;
			@call=split(/\s/, $_);
			if ($call[7]=~/MQRankSum=(-?\d+\.?\d*)/){
				if ($1<$MQRS){
					if ($call[6] eq "PASS"){$call[6]="MQRSFilter"}
					else {$call[6]=join(';', $call[6], "MQRSFilter")}
					$count[0]++;
					}
				}
#			else {print "Could not find MQRankSum information for variant $call[0]:$call[1]!\n"}

			if ($call[7]=~/ReadPosRankSum=(-?\d+\.?\d*)/){
				if ($1<$RPRS){
					if ($call[6] eq "PASS"){$call[6]="RPRSFilter"}
					else {$call[6]=join(';', $call[6], "RPRSFilter")}
					$count[1]++;
					}
				}
#			else {print "Could not find ReadPosRankSum information for variant $call[0]:$call[1]!\n"}

			if ($call[7]=~/DP=(\d+)/){
				if ($1>$dplimit){
					if ($call[6] eq "PASS"){$call[6]="DPFilter"}
					else {$call[6]=join(';', $call[6], "DPFilter")}
					$count[2]++;
					}
				}
			else {print "Could not find coverage information for variant $call[0]:$call[1]!\n"}

			$totalcount++;
			$modline=join("\t", @call);
			print $output "$modline\n";
			$counter--;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			}
		else {
			if ($filterline){
				print $output '##FILTER=<ID=MQRSFilter,Description="MQRankSum < ', $MQRS,'">', "\n";
				print $output '##FILTER=<ID=RPRSFilter,Description="ReadPosRankSum < ', $RPRS,'">', "\n";
				print $output '##FILTER=<ID=DPFilter,Description="DP > ', $dplimit,'">', "\n";
				$filterline=0;
				}
			print $output $_;
			}
		}

	close $vcffile;
	close $output;

	print "Total number of variants: $totalcount.\n";
	print "MQRSFilter applied to $count[0] (", $count[0]/$totalcount*100, "%) variants.\n";
	print "RPRSFilter applied to $count[1] (", $count[1]/$totalcount*100, "%) variants.\n";
	print "DPFilter applied to $count[2] (", $count[2]/$totalcount*100, "%) variants.\n";
	return;
	}


sub DPLimit{
	my ($filename, $threshold)=@_;
	my @header;
	my @call;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	my @depth;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";

	LINE: while (<$vcffile>){
		if(/^#CHROM/){
			@header=split("\t", $_);
			}
		elsif(!/^#/ && /^\w+\s/){
			@call=split(/\s/, $_);
			if ($call[7]=~/DP=(\d+)/){push @depth, $1}
			else {print "Could not find coverage information for variant $call[0]:$call[1]!\n"}
			$totalcount++;
			$counter--;
			unless ($counter){ print "$totalcount variants processed.\n"; $counter=$interval }
			}
		}

	close $vcffile;
	my ($mean, $sd)=meansd(\@depth);
	my $dplimit=$mean+$threshold*$sd;

	print "Total number of variants: $totalcount\n";
	print "Mean depth: $mean, standard deviation of depth: $sd, threshold value for depth: $dplimit\n";
	return $dplimit;
	}

sub DPLimitDoC{
	my ($docfile, $thininterval, $threshold)=@_;
	my $totalcount = 0;
	my $outinterval = 100000;
	my $counter = $outinterval;
	my $interval = $thininterval;
	my @depth;

	open my $input, "-|", "tabix $docfile ." or die "Could not open doc file $docfile!\n$!\n";	
	my @header = split("\t", <$input>);

	while (<$input>){
		unless (/^chr\d+/){ next }
		if ($interval--){ next }
		$interval = $thininterval;
		my @line = split("\t", $_);
		push @depth, $line[2];
		} continue { unless ($counter--){ print STDERR "Processed ", $totalcount+=$outinterval, " lines ...\n"; $counter = $outinterval } }

	close $input;
	my ($mean, $sd)=meansd(\@depth);
	my $dplimit=$mean+$threshold*$sd;

	print "Total number of variants: $totalcount\n";
	print "Mean depth: $mean, standard deviation of depth: $sd, threshold value for depth: $dplimit\n";
	return $dplimit;
	}

sub DPLimit_CRISP{
	my ($filename, $threshold)=@_;
	my @header;
	my @call;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	my @depth;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";

	LINE: while (<$vcffile>){
		if(/^#CHROM/){
			@header=split("\t", $_);
			}
		elsif(!/^#/ && /^\w+\s/){
			@call=split(/\s/, $_);
			if ($call[7]=~/MQS=(\d+),(\d+),(\d+),(\d+)/){ push @depth, $1+$2+$3+$4 }
			else {print "Could not find coverage information for variant $call[0]:$call[1]!\n"}
			++$totalcount;
			--$counter;
			unless ($counter){ print "$totalcount variants processed.\n"; $counter=$interval }
			}
		}

	close $vcffile;
	my ($mean, $sd)=meansd(\@depth);
	my $dplimit=$mean+$threshold*$sd;

	print "Total number of variants: $totalcount\n";
	print "Mean depth: $mean, standard deviation of depth: $sd, threshold value for depth: $dplimit\n";
	return $dplimit;
	}


sub QualLimit{
	my ($filename, $threshold)=@_;
	my @header;
	my @call;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	my @qual;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";

	LINE: while (<$vcffile>){
		if(/^#CHROM/){
			@header=split("\t", $_);
			}
		elsif(!/^#/ && /^\w+\s/){
			@call=split(/\s/, $_);
			if ($call[5]=~/(\d+\.?\d*)/){push @qual, $1}
			else {print "Could not find quality information for variant $call[0]:$call[1]!\n"}
			$totalcount++;
			$counter--;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			}
		}

	close $vcffile;
	my ($mean, $sd)=meansd(\@qual);
	my $cutoff=int($threshold*$totalcount);
	my @qualsorted=sort {$a <=> $b} @qual;
	my $limit=$qualsorted[-$cutoff];

	print "Total number of variants: $totalcount\n";
#	print "Mean variant quality: $mean, standard deviation of variant quality: $sd\n";
	print "Number of top variants selected: $cutoff\n";
	print "Threshold value for variant quality: $limit\n";
	return $limit;
	}


sub countgenotypes{
	my $filename=shift @_;
	my $totsize=shift @_;
	my $mindp=shift @_;

	my @header;
	my $startcol;
	my $endcol;
	my $homaltcount=0;
	my $tworeadscount=0;
	my $altdpcount=0;
	my $totalcount=0;
	my $interval=100000;
	my $counter=$interval;
	
        open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";

	while (<$vcffile>){
		chomp;
		if(/^#CHROM/){
			@header=split(/\s/, $_);
			$startcol=$#header-($totsize-1);
			$endcol=$#header;
			}

		elsif(!/^#/ && /^\w+\s/){
			my @call=split(/\s/, $_);
			my $homalt=0;
			my $tworeads=0;
			my $altdp=0;
			for my $i ($startcol..$endcol){
				my @gt=split(':', $call[$i]);
				if ($gt[0] eq "1/1"){$homalt=1}
				elsif ($gt[0] eq "0/1" || $gt[0] eq "1/0"){
					if ( (split(',', $gt[1]))[1]>=2 ){$tworeads=1}
					if ( $gt[2]>=$mindp ){$altdp=1}
					}
				}
			$homaltcount++ if $homalt;
			$tworeadscount++ if $tworeads;
			$altdpcount++ if $altdp;
			$totalcount++;
			$counter--;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			}
		}

	close $vcffile;
	print "Total number of variants: $totalcount\n";
	print "Number of variants with homozygote alternative genotype: $homaltcount\n";
	print "Number of variants with genotype with more than two alternative reads: $tworeadscount\n";
	print "Number of variants with alternative genotype and more than $mindp coverage: $altdpcount\n";
	return;
	}


sub TopVariants{
	my ($filename, $outfile, $minqual)=@_;
	my @header;
	my $filterline=0;
	my $totalcount=0;
	my $accepted=0;
	my $modline;
	my $interval=100000;
	my $counter=$interval;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";
	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";

	while (<$vcffile>){
		if(/^##FILTER/){
			$filterline=1;
			print $output $_;
			}
		elsif(/^#CHROM/){
			@header=split(/\s/, $_);
			print $output $_;
			}
		elsif(!/^#/ && /^\w+\s/){
			my @call=split(/\s/, $_);
			if ($call[6] eq '.'){$call[6]="PASS"}
			if ($call[5]<$minqual){
				if ($call[6] eq "PASS"){$call[6]="NotTopQual"}
				else {$call[6]=join(';', $call[6], "NotTopQual")}
				}
			else {$accepted++}
			$totalcount++;
			$modline=join("\t", @call);
			print $output "$modline\n";
			$counter--;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			}
		else {
			if ($filterline){
				print $output '##FILTER=<ID=NotTopQual,Description="Variant quality < ', $minqual,'">', "\n";
				$filterline=0;
				}
			print $output $_;
			}
		}
	close $vcffile;
	close $output;
	print "Total number of variants: $totalcount\n";
	print "Total number of accepted variants: $accepted\n";
	return;
	}


sub HardfilterVcf{
	my ($filename, $outfile, $totsize, $minmq, $mindist, $mindp, $minaltreads, $minaltgtqual)=@_;
	my @header;
	my $startcol;
	my $endcol;
	my $filterline=0;
	my $iscluster=0;
	my $nclusters=0;
	my $varincluster=0;
	my $homaltcount=0;
	my $minreadscount=0;
	my $mingtcount=0;
	my $singletoncount=0;
	my $totalcount=0;
	my @prevcall;
	my $modline;
	my $interval=100000;
	my $counter=$interval;
	
	open my $vcffile, "<", $filename or die "Could not open vcf file $filename!\n";
	open my $output, ">", $outfile or die "Could not open outfile $outfile!\n";

	while (<$vcffile>){
		if(/^##FILTER/){
			$filterline=1;
			print $output $_;
			}
		elsif(/^#CHROM/){
			@header=split(/\s/, $_);
			$startcol=$#header-($totsize-1);
			$endcol=$#header;
			print $output $_;
			}
		elsif(!/^#/ && /^\w+\s/){
			my @call=split(/\s/, $_);
			if ($call[7]=~/MQ=(\d+\.?\d*)/){
				if ($1<$minmq){
					if ($call[6] eq "PASS" || $call[6] eq '.'){$call[6]="MQFilter"}
					else {$call[6]=join(';', $call[6], "MQFilter")}
					}
				elsif ($call[6] eq '.'){$call[6]="PASS"}
				}
			if (defined $prevcall[0]){
				my $clustertest=0;
				$clustertest=1 if ($call[0] eq $prevcall[0] && ($call[1]-$prevcall[1])<$mindist);
				if ($clustertest){
					unless ($iscluster){
						$nclusters++;
						$varincluster++;
						}
					$varincluster++;
					$iscluster=1;
					}
				if ($iscluster){
					if ($prevcall[6] eq "PASS"){$prevcall[6]="ClusterFilter"}
					else {$prevcall[6]=join(';', $prevcall[6], "ClusterFilter")}
					}
				unless ($clustertest){$iscluster=0}
				$modline=join("\t", @prevcall);
				print $output "$modline\n";
				}

			my $homalt=0;
			my $minreads=0;
			my $mingtqual=0;
			my $altdp=0;
			my $first=1;
			my $issingleton=0;
			COLL: for my $i ($startcol..$endcol){
				my @gt=split(':', $call[$i]);
				if ($gt[0] eq "0/0" || $gt[0] eq "./."){next COLL}
				if ($gt[2]<$mindp){next COLL}
				my @alcov=split(',', $gt[1]);
				if (scalar(@alcov)>2){next COLL}
				if ($gt[0] eq "1/1" && $gt[3]>=$minaltgtqual){$homalt=1; $first=0}
				elsif ( ($gt[0] eq "0/1" || $gt[0] eq "1/0") && $gt[3]>=$minaltgtqual){
					if ($first){$issingleton=1; $first=0}
					else {$issingleton=0}
					if ($alcov[1]>=$minaltreads){$minreads=1}
					$mingtqual=1;					
					}
				}

			unless ($homalt){
				if ($call[6] eq "PASS"){$call[6]="HardFilter"}
				else {$call[6]=join(';', $call[6], "HardFilter")}
				}
			$homaltcount++ if $homalt;
			$minreadscount++ if $minreads;
			$mingtcount++ if $mingtqual;
			$singletoncount++ if $issingleton;
			$totalcount++;
			$counter--;
			unless ($counter){print "$totalcount variants processed.\n"; $counter=$interval}
			@prevcall=@call;
			}
		else {
			if ($filterline){
				print $output '##FILTER=<ID=MQFilter,Description="MQ < ', $minmq,'">', "\n";
				print $output '##FILTER=<ID=HardFilter,Description="No homozygote alternative genotype">', "\n";
				print $output '##FILTER=<ID=ClusterFilter,Description="Located within ', $mindist, ' bp of other SNP">', "\n";
				$filterline=0;
				}
			print $output $_;
			}
		}
	$modline=join("\t", @prevcall);
	print $output "$modline\n";
	close $vcffile;
	close $output;
	print "Total number of variants: $totalcount\n";
	print "Number of variants with homozygote alternative genotype: $homaltcount\n";
	print "Number of variants with genotype with more than $minaltreads alternative reads: $minreadscount\n";
	print "Number of variants with high quality alternative heterozygote genotype: $mingtcount\n";
	print "Number of variants in clusters: $varincluster, in a total of $nclusters clusters.\n";
	print "Percentage of singleton sites: ", $singletoncount/$totalcount, "\n";
	return;
	}


sub show{
	my $dumper=Dumpvalue->new(tick => q("), compactDump => 1, veryCompact => 1);
	$dumper->dumpValue(@_);
	}


sub meansd{
	my $ref=shift @_;
	my $rt=0;
	my $sum2=0;
	my $am;
	my $sd;

	for my $i (@$ref){
		$rt+=$i;
		}

	$am=$rt / scalar(@$ref);

	for my $j (@$ref){
		$sum2+=($j-$am)**2;
		}

	$sd=sqrt($sum2 / scalar(@$ref));
	return($am, $sd);
	}

sub meansd_stream{
	my $docfile = shift;
	my $thinint = shift;
	open my $input, "-|", "tabix $docfile ." or die "Could not open doc file $docfile!\n$!\n";

	my $totvariants = 0;
	my $interval = $thinint;
	my @header = split("\t", <$input>);
	my @line = split("\t", <$input>);
	my $xk = $line[2];
	my $k = 1;
	my $Mk = $xk;
	my $Qk = 0;

	while (<$input>){
		unless (/^chr\d+/){ next }
		if ($interval--){ next }
#      print STDERR "Processed ", $totvariants+=$thinint, " lines ...\n";
		$interval = $thinint;
		@line = split("\t", $_);
		$xk = $line[2];
		++$k;
		my $d = $xk-$Mk;
		$Qk += ($k-1)*$d*$d/$k;
		$Mk += $d/$k;
		}

	return($Mk, sqrt($Qk/$k), $k);
	}

