package Chrommap;

use strict;
use warnings;
#use functions;


sub new{
	my ($class, $chrom)=@_;
	my $self={};
	$self->{_Chrom}=defined $chrom ? $chrom : {};
	bless $self, $class;
	return $self;
	}


sub getValue{
	my $self=shift;
	my ($chrom, $scaff, $key)=@_;

	if (scalar(@_)==1){return($self->{_Chrom}{$chrom})}
	elsif (scalar(@_)==2){return($self->{_Chrom}{$chrom}{$scaff})}
	elsif (scalar(@_)==3){return($self->{_Chrom}{$chrom}{$scaff}{$key})}
	else {die "Invalid number of arguments!\n"}
	}


sub setValue{
	my $self=shift;
	unless (@_==4){die "Wrong number of arguments!\n"}
	my ($chrom, $scaff, $key, $value)=@_;

	$self->{_Chrom}{$chrom}{$scaff}{$key}=$value;

	return($self->{_Chrom}{$chrom}{$scaff}{$key});
	}


sub allKeys{
	my $self=shift;
	my ($chrom, $scaff)=@_;

	if (defined $scaff){return(keys %{$self->{_Chrom}{$chrom}{$scaff}})}
	elsif (defined $chrom){return(keys %{$self->{_Chrom}{$chrom}})}
	else {return(keys %{$self->{_Chrom}})}
	}


sub allValues{
	my $self=shift;
	my ($chrom, $scaff)=@_;

	if (defined $scaff){return(values %{$self->{_Chrom}{$chrom}{$scaff}})}
	elsif (defined $chrom){return(values %{$self->{_Chrom}{$chrom}})}
	else {return(values %{$self->{_Chrom}})}
	}


sub getScaffolds{
	my ($self, @chrom)=@_;
	my @scaffolds;
	for my $chrom (@chrom){
		unless ( exists $self->{_Order}{$chrom} ){ warn "Chromosome $chrom not found in chromosome map! Skipping chromosome.\n"; next } 
		push @scaffolds, @{ $self->{_Order}{$chrom} };
		}
	return(@scaffolds);
	}

sub getChrom{
	my ($self, $scaff)=@_;
	if ( defined $self->{_Scaffolds}{$scaff} ){
		return( $self->{_Scaffolds}{$scaff}{'chrom'} );
		}
	else {
		warn "Warning: scaffold $scaff has not been located on a chromosome!\n";
		return;
		}
	}

sub findPosonChrom{
	my ($self, $scaff, $pos)=@_;

	my ($chrom, $trpos);
	if ( defined $self->{_Scaffolds}{$scaff} ){
		$chrom=$self->{_Scaffolds}{$scaff}{'chrom'};
		my $dir=$self->{_Scaffolds}{$scaff}{'dir'};
		if ($dir ne '+' && $dir ne '-') { print "Warning: no information about direction for scaffold $chrom: $scaff, assuming forward direction!\n" }
		if ($dir eq '-'){
			$trpos=$self->{_Scaffolds}{$scaff}{'end'}-$pos;
			}
		else {
			$trpos=$pos+$self->{_Scaffolds}{$scaff}{'start'};
			}
		}
	else {
		warn "Warning: scaffold $scaff has not been located on a chromosome!\n";
		}

	return($chrom, $trpos);
	}


sub readChromMap{
	my ($self, $chrommap, $gapsize, $keepinfo)=@_; # bed file format with chromosome locations
	my %mapping;
	my %scaffolds;

	open my $bed, "<", $chrommap or die "Could not open bed file $chrommap!\n";
	my $chrompos=0;
	my $prevchrom="";
	my ($start, $end);

	while (<$bed>){
		chomp;
		if (/^([^\n\s]+)\s+([^\n\s]+)\s+(\d+)\s+([\+\-\?])\s*(.*)/){
			push @{ $self->{_Order}{$1} }, $2;
			$self->{_Lengths}{$2}=$3;
			unless ($1 eq $prevchrom){$chrompos=0; $prevchrom=$1}
			$mapping{$1}{$2}={'start'=>$chrompos, 'end'=>$chrompos+$3-1, 'dir'=>$4};
			$scaffolds{$2}={'chrom'=>$1, 'start'=>$chrompos, 'end'=>$chrompos+$3-1, 'dir'=>$4};
			if ($keepinfo==1){$mapping{$1}{$2}{'info'}=$5}
#			print STDERR "Scaffold $2 on is on chromosome $1, ", $chrompos,  "-", $chrompos+$3, ".\n";
			$chrompos+=$gapsize+$3;
			}
		else {print "Warning: wrong format on line $. of chromosome map $chrommap!\n"}
		}
	close $bed;
	$self->{_Chrom}=\%mapping;
	$self->{_Scaffolds}=\%scaffolds;
	return;
	}


sub locifromChromMap{
	my ($self, $maxsize, @chromosomes)=@_;
	my $maxset=defined $maxsize ? 1 : 0;
	if (@chromosomes){
		@chromosomes=keys %{ $self->{_Chrom} } if ($chromosomes[0] eq 'all');
		}
	else { @chromosomes=keys %{ $self->{_Chrom} } }
	my %loci;

	CHROM: for my $chrom (@chromosomes){
		unless ( exists $self->{_Order}{$chrom} ){ warn "Chromosome $chrom not found in chromosome map! Skipping chromosome.\n"; next CHROM } 
		SCAFF: for my $scaffold ( @{$self->{_Order}{$chrom}} ){
			if ( defined $loci{$scaffold} ){ warn "Scaffold $scaffold on chromosome $chrom appreared twice in chromosome map! Skipping scaffold.\n"; next SCAFF }
			my $seqlength=$self->{_Lengths}{$scaffold};
			my $dir=$self->{_Chrom}{$chrom}{$scaffold}{'dir'};
			unless ($maxset){$maxsize=$seqlength}

			if ($dir eq '-'){
				my $start;
				my $end=$seqlength+1;
				while ($end>1){
					$start=$end-1;
					$end-=$maxsize;
					if ($end<1){$end=1}
					push @{$loci{$scaffold}}, {'start'=>$end-1, 'end'=>$start-1};
					}
				}
			else {
				my $start;
				my $end=0;
				while ($end<$seqlength){
					$start=$end+1;
					$end+=$maxsize;
					if ($end>$seqlength){$end=$seqlength}
					push @{$loci{$scaffold}}, {'start'=>$start-1, 'end'=>$end-1};
					}
				}
			unless ($dir eq '+' || $dir eq '-') {print "Warning: no information about direction for scaffold $chrom: $scaffold, assuming forward direction!\n"}
			}
		}
	my $regions=Regions->new(\%loci);
	return($regions);
	}


sub readChromMap_old{
	my ($self, $chrommap, $format, $keepinfo)=@_; # bed file format with chromosome locations
	my %mapping;

        open my $bed, "<", $chrommap or die "Could not open bed file $chrommap!\n";

	while (<$bed>){
		chomp;

		if (/^([^\n\t]+)\t(\d+)\t(\d+)\t([^\s\n]+)\t?(.*)/){
			if ($format eq 'bed'){$mapping{$1}{$4}={'start'=>$2, 'end'=>($3-1)}}
			elsif ($format eq 'onebased'){$mapping{$1}{$4}={'start'=>($2-1), 'end'=>($3-1)}}
			else {$mapping{$1}{$4}={'start'=>$2, 'end'=>$3}}
			if ($keepinfo==1){$mapping{$1}{$4}{'info'}=$5}
			}
		}

	close $bed;
	$self->{_Chrom}=\%mapping;
	return;
	}



1;

