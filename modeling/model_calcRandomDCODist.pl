#!/usr/bin/perl

use strict;

my $trials = 100000;

my $dsbmax = 26;
my $dsbmin = 18;

# Expected based on 2016 data
my $COPerArm = 0.55;  # based on 541 COs / 980 arms
my $NCOPerArm = 0.56; # based on 549 NCO-gcs / 980 arms

#my $stocks = 52;
#my $armsPerStock = 2;

#my $totalArms = $stocks * $armsPerStock;
#my $totalCOs = $totalArms * $COPerArm; # cos expected
#my $totalNCOs = $totalArms * $NCOPerArm; #ncos expected

my %snpDensity;
my @chrList = qw/ chrX chr2L chr2R chr3L chr3R chr4 /;
#my @chrList = qw/ chr2L chr2R /;

## outputs
my(%coCounts,%ncoCounts);

my %chrs;
$chrs{1} = "chrX";
$chrs{2} = "chr2L";
$chrs{3} = "chr2R";
$chrs{4} = "chr3L";
$chrs{5} = "chr3R";
$chrs{6} = "chr4";

my %chromosomeSizes;
$chromosomeSizes{"chrX"} = 23542271;
$chromosomeSizes{"chr3R"} = 32079331;
$chromosomeSizes{"chr3L"} = 28110227;
$chromosomeSizes{"chr2L"} = 23513712;
$chromosomeSizes{"chr2R"} = 25286936;
$chromosomeSizes{"chr4"} = 1348131;

my %chromosomeFraction;
my $genomeSize = 23542271 + 32079331 + 28110227 + 23513712 + 25286936 + 1348131;
my $low = 0;
my $newHigh = 0;
foreach my $chr (keys %chromosomeSizes) {
	my $fraction = sprintf("%0.3f", $chromosomeSizes{$chr} / $genomeSize);
	$fraction *= 1000;
	$newHigh += $fraction;
	foreach ($low..$newHigh) {
		$chromosomeFraction{$_} = $chr;
	}
	#print "$chr\t$low - $newHigh\n";
	$low = $newHigh + 1;
}

my(%coCounts,%dcoTotal);
foreach my $trial (1..$trials) {
	my(%modeledCOs);
	#print "$trial\n" if $trial =~ /000/;

	my $dsbs = int(rand($dsbmax - $dsbmin)) + $dsbmin;

	# Distribute DSBs
	foreach my $dsbNum (1..$dsbs) {
		my $chrNum = int(rand(1000)) + 1;
		my $chrName = $chromosomeFraction{$chrNum};
		my $maxbp = $chromosomeSizes{$chrName};
		my $pos = int(rand($maxbp)) + 1;

		# is this a CO? assume 6 total COs per genome (0.5 too high)
		my $coChance = int(rand($dsbs)) + 1;
		next if $coChance >= 7;

		# Do we recover the chromatid the DSB is on? 50% chance for DCOs
		my $chromatid = int(rand(4)) + 1;
		next if $chromatid >= 3;

		#print "$dsbNum\t$chrName\t$pos\t$chromatid\n";
	
		$modeledCOs{$chrName} .= "$pos,";
	}

	foreach my $chr (@chrList) {
		if ($modeledCOs{$chr} =~ /\d+\,\d+\,\d+\,\d+\,/) {
			$coCounts{$chr}{"QCO"} += 3; # 2 doubles 
			my($co1,$co2,$co3,$co4) = $modeledCOs{$chr} =~ /(\d+)\,(\d+)\,(\d+)\,(\d+)/;

			my @cos = qw/ $co1 $co2 $co3 $co4 /;
			my @sortedcos = sort { $a <=> $b } @cos;
			my($co1,$co2,$co3,$co4) = ($sortedcos[0],$sortedcos[1],$sortedcos[2],$sortedcos[3]);

			my $span = abs($co1-$co2);
			$dcoTotal{$chr}{"QCO"} += $span;
			my $span = abs($co2-$co3);
			$dcoTotal{$chr}{"QCO"} += $span;
			my $span = abs($co3-$co4);
			$dcoTotal{$chr}{"QCO"} += $span;

		} elsif ($modeledCOs{$chr} =~ /\d+\,\d+\,\d+\,/) {
			$coCounts{$chr}{"TCO"} += 2; # 2 doubles 
			my($co1,$co2,$co3) = $modeledCOs{$chr} =~ /(\d+)\,(\d+)\,(\d+)/;

			my @cos = qw/ $co1 $co2 $co3 /;
			my @sortedcos = sort { $a <=> $b } @cos;
			my($co1,$co2,$co3) = ($sortedcos[0],$sortedcos[1],$sortedcos[2]);

			my $span = abs($co1-$co2);
			$dcoTotal{$chr}{"TCO"} += $span;
			my $span = abs($co2-$co3);
			$dcoTotal{$chr}{"TCO"} += $span;

		} elsif ($modeledCOs{$chr} =~ /\d+\,\d+\,/) {
			$coCounts{$chr}{"DCO"}++; 

			my($low,$high) = $modeledCOs{$chr} =~ /(\d+)\,(\d+)/;
			my $span = abs($low-$high);

			$dcoTotal{$chr}{"DCO"} += $span;

		} elsif ($modeledCOs{$chr} =~ /\d+\,/) {
			#$coCounts{"SCO"}++; 
		} else {
			#$coCounts{"NCO"}++; 
		}
	}
}

print "Counts:\n";
foreach my $chr ("chrX","chr2L","chr2R") {
	my $dcoCount = $coCounts{$chr}{"DCO"};
	my $dcoBases = $dcoTotal{$chr}{"DCO"};
	my $ave = sprintf("%0.1f", $dcoBases / $dcoCount);
	print "$chr\t$dcoCount\t$ave\n";

	$dcoCount = $coCounts{$chr}{"DCO"} + $coCounts{$chr}{"TCO"} + $coCounts{$chr}{"QCO"};
	$dcoBases = $dcoTotal{$chr}{"DCO"} + $dcoTotal{$chr}{"TCO"} + $dcoTotal{$chr}{"QCO"};
	my $ave = sprintf("%0.1f", $dcoBases / $dcoCount);
	print "$chr\t$dcoCount\t$ave <- includes TCOs and QCOs\n";
}
