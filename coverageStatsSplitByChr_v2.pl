#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use POSIX;

my %options;
getopts("i:", \%options);
my $infile = $options{i} or &usage;
sub usage {die "USAGE: " . basename($0) . " [-i BAM infile (with index)]\n";} 

# Ensure bam file has an index.
die "ERROR: Can't find an index file for '$infile'.\nPlease ensure '$infile.bai' is in same directory as '$infile'.\n" 
	if !-e "$infile.bai";

# Do not proceed if samtools lacks idxstats command.
die "ERROR: You are using an older version of samtools. Please upgrade.\n" 
	if !&isRightSamtools;  

# Get sizes of chromsomes within reference.
my $refSizes = getRefSizes($infile);

# Extract raw coverage data.
open (IN, " /pub9/laura/samtools-0.1.16/samtools pileup $infile |") or die "ERROR: could not open $infile.\n";

# Header for outfile.
print "Ref id\tRef size\tMapped\t% Mapped\tMean\tMedian\tSD\tQ1\tQ3\t2.5%\t97.5%\tMin\tMax\n";


my $currentChr = "";
my @data;
while (<IN>) {
	/^(\S+)\t(\d+)\t(\w)\t(\d+)/ or die "ERROR: regex failed with line '$_'.\n";
	my $chr 	= $1;
	my $position 	= $2;
	my $coverage	= $4;

	if ($chr ne $currentChr) {
		process(\@data, $currentChr, $refSizes) if @data;
		undef @data;
		$currentChr = $chr;
		}

	push(@data, $coverage);
	

	}

process(\@data, $currentChr, $refSizes);

close (IN) or die "ERROR: could not close $infile.\n";

################################################################################

# Checks that correct version of samtools is installed (insofar as idxstats
# is required). 

sub isRightSamtools {
	
	open (IN, "samtools 2>&1 | ") or die "ERROR: could not pipe samtools.\n";

	my $isOk = 0;
	while (<IN>) {
		chomp;
		if (/^\s+idxstats\s+/) {
			$isOk = 1;
			last;
			}
		}

	return $isOk;

	} # End of method.

################################################################################

# Determine reference size from accompanying bam index file.
sub getRefSizes {
	
	my $infile = shift;

	my %sizes;

	open (PIPE, " samtools idxstats $infile |") or die "ERROR: could not open pipe.\n";
	while (<PIPE>) {
		chomp;
		/^(\S+)\s+(\d+)/ or die "ERROR: Regex failed with line '$_'.\n";
		$sizes{$1} = $2;
		}
	close (PIPE) or die "ERROR: could not close pipe.\n";

	return \%sizes;

	} # End of method.

################################################################################

sub printData {	
	
	my $data = shift;
	my $chr = shift;

	my $outfile = "$chr\_coverage.dat";

	open (OUT, ">$outfile") or die "ERROR: could not create $outfile.\n";
	foreach (@$data) {
		print OUT "$_\n";	
		}
	close (OUT) or die "ERROR: could not close $outfile.\n";

	} # End of method.

################################################################################

sub process {

	my $data = shift;
	my $chr = shift;
	my $refSizes = shift;

	if (!exists $refSizes->{$chr}) {
		die "ERROR: '$chr' could not be found within reference file.\n";
		}

	@$data = sort {$a<=>$b} @$data;
	
	#printData($data, $chr);

	my $n 			= scalar(@$data);
	my $sum 		= getSum($data);
	my $n50 		= getN50($data, $sum);
	my $mean		= $sum/$n;
	my $deviates 		= getDeviates($data, $mean);
	my $squaredDeviates 	= getSquares($deviates);
	my $sumOfSquares 	= getSum($squaredDeviates);
	my $popVar		= $sumOfSquares/$n; 
	my $popSD 		= sqrt($popVar);
	my $sampleVar		= $sumOfSquares/($n-1);
	my $sampleSD		= sqrt($sampleVar);	
	#my $mode 		= getMode($data);
	my $median		= getMedian($data);
	# Alternative: calculated with R type 6 (Minitab) quantile method.
	#my $median		= getQuantile($data, 0.5);
	my $q1			= getQuantile($data, 0.25);
	my $q3			= getQuantile($data, 0.75);
	my $lowerPercentile	= getQuantile($data, 0.025);
	my $upperPercentile	= getQuantile($data, 0.975);
	my $min			= $data[0];
	my $max 		= $data[$n-1];
	my $range		= $max-$min;

	printAsLine(
		$chr, 
		$refSizes->{$chr}, 
		$n, 
		$mean, 
		$median, 
		$popSD, 
		$q1, 
		$q3,
		$lowerPercentile,
		$upperPercentile, 
		$min, 
		$max
		);

	undef @$data;

	} # End of method.

################################################################################

sub printAsLine {
	
	my $chr 		= shift;
	my $chrSize 		= shift;
	my $n 			= shift;
	my $mean 		= shift;
	my $median 		= shift;
	my $popSD 		= shift;
	my $q1 			= shift;
	my $q3 			= shift;
	my $lowerPercentile	= shift;
	my $upperPercentile 	= shift;
	my $min 		= shift;
	my $max 		= shift;

	printf(
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        	$chr,		                 	# chr id
       		&format($chrSize),                      # size of chr
        	&format($n),                            # mapped region
        	&roundup(&format($n/$chrSize*100)),	# % mapped region
        	&roundup(&format($mean)),               # mean coverage
        	&roundup(&format($median)),             # median coverage
        	&roundup(&format($popSD)),              # sd coverage
        	&roundup(&format($q1)),                 # coverage lower quartile
        	&roundup(&format($q3)),                 # coverage upper quartile
		&roundup(&format($lowerPercentile)),	# coverage at 2.5% percentile
		&roundup(&format($upperPercentile)),	# coverage at 97.5% percentile
        	&format($min),                          # minimum coverage
       	 	&format($max)                           # maximum coverage
        	);


	} # End of method.

################################################################################

sub roundup {
        return sprintf("%.2f", unformat(shift));
        } # End of method.

################################################################################

sub unformat {
        my $x = shift;
        $x =~ s/,//g;
        return $x;
        } # End of method.

################################################################################

sub format {
        my $x = reverse(shift);
        $x =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
        return scalar reverse $x;
        } # End of method.

################################################################################

sub getMode {
	
	my $data = shift;

	my %map;
	foreach (@$data) {
		$map{$_}++;
		}

	my $mode =  (sort {$map{$b} <=> $map{$a}} keys %map)[0];

	if ($map{$mode} == 1) {
		return "na";
		} else {
		return $mode;
		}

	} # End of method.

################################################################################

sub getN50 {
	
	my $sortedData = shift;
	my $sum = shift;

	my $halfSum = $sum / 2;

	for (my $i = @$sortedData-1; $i >= 0; $i--) {
	
		$sum -= $sortedData->[$i];
		return $sortedData->[$i] if $halfSum > $sum;
		}

	die "ERROR: failed to find N50.\n";

	} # End of method.

################################################################################

sub getSum {
	
	my $data = shift;
	
	my $x = 0;
	foreach (@$data) {
		$x += $_;
		}

	return $x;

	} # End of method.

################################################################################

sub getDeviates {
	
	my $data = shift;
	my $mean = shift;

	my @deviates;

	foreach (@$data) {
		push(@deviates, $_ - $mean);
		}
	
	return \@deviates;

	} # End of method.

################################################################################

sub getSquares {
	
	my $data = shift;
	
	my @squares;

	foreach (@$data) {
		push(@squares, $_**2);
		}	

	return \@squares;

	} # End of method.

################################################################################

sub getMedian {
	
	my $sortedData = shift;

	if (@$sortedData % 2) {
		return $sortedData->[@$sortedData/2];
		}
	else {
		return ($sortedData->[@$sortedData/2-1] + $sortedData->[@$sortedData/2]) / 2;
		} 

	} # End of method.

################################################################################

sub getQuantile {
	
	my $sortedData = shift;
	my $quantile = shift;
	
	my $n = @$sortedData;

	# Multiply quantile by number of observations plus 1.
	my $number = $quantile * ($n + 1);
	
	# Next step depends on whether number is a whole number or not.
	my $fraction = getFraction($number);

	if ($fraction == 0) {
		# number is a whole number. Therefore our required value already
		# exists within the observation collection.
		return $sortedData->[$number - 1];
		}
	elsif ($number >= $n) {
		# Number represents largest observation in collection. 
		return $sortedData->[$number - 1];
		}
	elsif ($number <= 0) {
		# Number represents smallest number in collection.
		return $sortedData->[0];
		}
	else {
		# number is a fraction number. Therefore our required value is the
		# weighted average of two values within the observation collection.
		my $x1 = $sortedData->[floor($number) - 1];
		my $x2 = $sortedData->[floor($number)];

		return (1-$fraction) * $x1 + ($fraction) * $x2;
		}

	} # End of method.

################################################################################

sub getFraction {
	
	my $number = shift;
	
	if ($number =~ /^\d+$/) {
		return 0;
		} 
	else {
		$number =~ /\.(\d+)$/;
		return "0.$1";
		}

	} # End of method.

################################################################################
1;
