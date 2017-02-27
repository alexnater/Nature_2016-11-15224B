#!/aplic/perl-5.16.2/bin/perl -w

use strict;
use Getopt::Long;
use Switch;

my %opt=();

my ($program, $pversion, $pdescription, $pusage);

$program = $0;
$program =~ s/^.*\///;

$pversion='1';
$pdescription = "$program (ver:$pversion) 
		Description of the program: from a given vcf file it generates a matrix,
		each line corresponds to a SNP, each column to a sample,
		it gives 0 if it is like the reference, 1 if it has a different allele\n";
$pusage="$0 ";

if (! defined $ARGV[0]){
  print "
    DESCRIPTION $pdescription 
    USAGE $pusage 
      -vcfF	vcf file with the SNPs, it should have the genotypes of each sample and not contain INDELS
      -mtAl	if it is passed, multiAllelic SNPs are ignored
		optional, by default they are not ignored
      -minDP	[int] minimum RD for each sample to take into account its genotypes
		optional, by default no restriction
      -maxDP	[int] maximum RD for each sample to take into account its genotypes
		optional, by default no restriction
      -excNA	if it is passed, snps with a sample with genotype as ./. or with less/more RD than -minDP/-maxDP in case they were included 
		optional, by default they are not ignored
      -excDiv	if it is passed, snps fixed in all samples are removed, after applying the other filters (those are the ones that are fixed but different from the reference
		optional, by default they are not ignored
      -outF\n";
  exit;
}

if (! &GetOptions(\%opt, "vcfF:s", "mtAl", "minDP:i", "maxDP:i","excNA","excDiv",
                         "outF:s")){
  die "  Command line unsucessfully parsed!\n";
}


#CHECKING THE COMMAND and INPUTS
$opt{'vcfF'} ||= die "  Insert a vcf file with the genotypes of the sample with the command -vcfF\n";
$opt{'outF'} ||= die "  Insert the output file with the command -outF\n";
$opt{'minDP'} ||= 0;
die "  -minDP should be an integer\n" if ($opt{'minDP'} !~ /^-?\d+\z/); 
die "  -maxDP should be an integer\n" if (defined $opt{'maxDP'} && $opt{'maxDP'} !~ /^-?\d+\z/); 

open (IN, $opt{'vcfF'}) or die "could not open $opt{'vcfF'}\n";
open (OUT, ">".$opt{'outF'}) or die "could not open $opt{'outF'}\n";
print OUT "CHROM\tPOS";

my $numColFormat;
my $numColAlt;
my $numGTfirst;
my $numDPfirst;
my $contMA=0;
my $contExc=0;
my $contExcDiv=0;
my $cont=0;
while (my $line = <IN>){
	chomp ($line);
	my @columns = split (/\t/,$line);
  
	next if $line =~ /^##/;
	if ($line =~ /^#CHROM/ && !defined $numColFormat){
		foreach my $i (0..$#columns){
                        switch ($columns[$i]){
                                case "ALT" {$numColAlt=$i;}
                                case "FORMAT" {$numColFormat=$i;}
                        }
                }
                die "  Not found all columns in $opt{'vcfF'}\n" if (!defined $numColFormat || !defined $numColAlt);

		foreach my $i (($numColFormat+1)..$#columns){
			print OUT "\t".$columns[$i];
		}
		print OUT "\n";

                next;
	} 

	my @form = split (/:/, $columns[$numColFormat]);
	my $numGT;
	my $numDP;
	foreach my $i (0..$#form){
        	if ($form[$i] =~ /GT/){
			$numGT = $i;
			$numGTfirst = $i if !defined $numGTfirst;
		}
        	if ($form[$i] =~ /DP/){
			$numDP = $i;
			$numDPfirst = $i if !defined $numDPfirst;
		}
	}
        die "  Not found field \"GT\" in $opt{'vcfF'} in line $line\n" if (!defined $numGT);
	die "  Found different position for field \"GT\" in line $line than the rest\n" if (defined $numGTfirst && $numGTfirst != $numGT);
        die "  Not found field \"DP\" in $opt{'vcfF'} in line $line\n" if ($opt{'minDP'} >0 && !defined $numDP);
	die "  Found different position for field \"DP\" in line $line than the rest\n" if ($opt{'minDP'} >0 && defined $numDPfirst && $numDPfirst != $numDP);

	my @alts = split (/,/, $columns[$numColAlt]);
	$contMA ++ if scalar @alts >1;
	$cont++;

	###Removed the multiallelic filter:
	#next if (defined $opt{'mtAl'} && scalar @alts >1);

	my $snp = "$columns[0]\t$columns[1]";
	
	my $numSamp=0;
	my $numCh=0;
	foreach my $i (($numColFormat+1)..$#columns){
		my @samp = split (/:/, $columns[$i]);
		$samp[$numGT] =~ s/[|\/]//g;
			
		if ($samp[$numGT] =~/^\.+$/){
			$snp.="\tNA";
		}
		elsif ($opt{'minDP'} >0 && $samp[$numDP] < $opt{'minDP'}){
			$snp.="\tNA";
		}
		elsif (defined $opt{'maxDP'} && $samp[$numDP] > $opt{'maxDP'}){
			$snp.="\tNA";
		}
		elsif ($samp[$numGT] =~/^0+$/){
			$snp.="\t0";
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^01$/){
			$snp.="\t1";
			$numCh++;
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^11$/){
			$snp.="\t2";
			$numCh++;
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^02$/){
			$snp.="\t3";
			$numCh++;
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^22$/){
			$snp.="\t4";
			$numCh++;
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^03$/){
			$snp.="\t5";
			$numCh++;
			$numSamp++;
		}
		elsif ($samp[$numGT] =~/^33$/){
			$snp.="\t6";
			$numCh++;
			$numSamp++;
		}
		if (!($samp[$numGT] eq "00" || $samp[$numGT] eq "01" || $samp[$numGT] eq "10" || $samp[$numGT] eq "11")){
			#print $line;
		}
	}
	$snp.= "\n";

	if (defined $opt{'excNA'} && $snp =~ /NA/){
		$contExc ++;
	}
	elsif (defined $opt{'excDiv'} && ($numCh==0 || $numCh == $numSamp)){
		$contExcDiv ++;
	}
	else{
		print OUT $snp;
	}

}
close IN;
close OUT;

print "Number of SNPs = $cont; number of multiallelic = $contMA; number of excluded SNPs due to a bad genotype = $contExc; number of SNPs that are not informative = $contExcDiv\n";

