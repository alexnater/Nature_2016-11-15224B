package Param;

use strict;
use warnings;
use Math::Random qw(random_uniform random_uniform_integer random_normal random_gamma);

# written by Alexander Nater, May 2013

sub new{
	my $class=shift;
	my $self={
		_Name      => shift,
		_Int       => shift,
		_Log       => shift,
		_Dist      => shift,
		};

	bless $self, $class;
	return $self;
	}


sub setMinMax{
	my ($self, $min, $max)=@_;
	$self->{_Min}=$min;
	$self->{_Max}=$max;
#	print "Min of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Min}\n";
#	print "Max of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Max}\n";
	return $self->{_Min}, $self->{_Max};
	}


sub setMeanStd{
	my ($self, $mean, $std)=@_;
	$self->{_Mean}=$mean;
	$self->{_Std}=$std;
#	print "Mean of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Mean}\n";
#	print "Standard deviation of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Std}\n";
	return $self->{_Mean}, $self->{_Std};
	}


sub setMeanShape{
	my ($self, $mean, $shape)=@_;
	$self->{_Mean}=$mean;
	$self->{_Shape}=$shape;
#	print "Mean of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Mean}\n";
#	print "Shape parameter of $self->{_Dist} distribution of parameter $self->{_Name} is $self->{_Shape}\n";
	return $self->{_Mean}, $self->{_Shape};
	}


sub getRandomNumbers{
	my ($self, $n)=@_;

	if ($self->{_Dist} eq "uniform" && !$self->{_Log}){
		my @a;
		if ($self->{_Int}){
			@a=random_uniform_integer($n, $self->{_Min}, $self->{_Max});
			}
		else {
			@a=random_uniform($n, $self->{_Min}, $self->{_Max});
			}
		return \@a;		
		}

	elsif ($self->{_Dist} eq "uniform" && $self->{_Log}){
		my @a;
		my $logmin=logn(10, $self->{_Min});
		my $logmax=logn(10, $self->{_Max});
		if ($self->{_Int}){
			@a=random_uniform_integer($n, $logmin, $logmax);
			}
		else {
			@a=random_uniform($n, $logmin, $logmax);
			}
		my $out_ref=powerarray(10, \@a);
		return $out_ref;
		}


	elsif ($self->{_Dist} eq "normal" && !$self->{_Log}){
		my @a;
		if ($self->{_Int}){
			@a=random_normal($n, $self->{_Mean}, $self->{_Std});
			}
		else {
			@a=random_normal($n, $self->{_Mean}, $self->{_Std});
			}
		return \@a;
		}

	elsif ($self->{_Dist} eq "normal" && $self->{_Log}){
		my @a;
#		my $logmean=logn(10, $self->{_Mean});
#		my $logstd=logn(10, $self->{_Std});
		if ($self->{_Int}){
			@a=random_normal($n, $self->{_Mean}, $self->{_Std});
			}
		else {
			@a=random_normal($n, $self->{_Mean}, $self->{_Std});
			}
		my $out_ref=powerarray(10, \@a);
		return $out_ref;
		}

	elsif ($self->{_Dist} eq "gamma" && !$self->{_Log}){
		my @a;
		if ($self->{_Int}){
			@a=random_gamma($n, ($self->{_Shape}/$self->{_Mean}), $self->{_Shape});
			}
		else {
			@a=random_gamma($n, ($self->{_Shape}/$self->{_Mean}), $self->{_Shape});
			}
		return \@a;
		}

	elsif ($self->{_Dist} eq "gamma" && $self->{_Log}){
		my @a;
#		my $logmean=logn(10, $self->{_Mean});
#		my $logshape=logn(10, $self->{_Shape});
		if ($self->{_Int}){
			@a=random_gamma($n, ($self->{_Shape}/$self->{_Mean}), $self->{_Shape});
			}
		else {
			@a=random_gamma($n, ($self->{_Shape}/$self->{_Mean}), $self->{_Shape});
			}
		my $out_ref=powerarray(10, \@a);
		return $out_ref;
		}

	else {die "unknown distribution!"}

#	return \@a;
	return 0;
	}



# helper functions

sub power{
	my $base=shift @_;
	my $exp=shift @_;
	return $base**$exp;
	}


sub powerarray{
	my $base=shift @_;
	my $ref_exp=shift @_;

	my @result=();

	for my $exp (@{$ref_exp}){
		push @result, ($base**$exp);
		}

	return \@result;
	}


sub logn{
	my $base=shift @_;
	my $number=shift @_;
	return log($number)/log($base);
	}


sub lognarray{
	my $base=shift @_;
	my $ref_number=shift @_;

	my @result=();

	for my $number (@{$ref_number}){
		push @result, (log($number)/log($base));
		}

	return \@result;
	}


sub meansd{
	my $self=shift @_;
	my $ref=shift @_; # reference to array
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

1;

