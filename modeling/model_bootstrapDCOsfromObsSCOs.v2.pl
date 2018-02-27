#!/usr/bin/perl

use strict;
use Getopt::Std;

# bootstrap DCOs from SCO events, use data from 2016
# Is the SCO data from the interchromosomal effect differnt from wt?

#my %opts;
#getopts('c:h', \%opts); 

## Usage Output
#if ($opts{'h'} || !$opts{'c'}) {
	#print "
	#Required:
		#-c Chromosome (chrX, or chr2)
	#
	#Optional:
		#-h This helpful help
	#\n";
	#exit 0;
#}


my $trials = 100000;

my $bootstrapOutput = "chr\tpos\n";
foreach my $chr ("chrX","chr2") {
	my(%scoList,%scoCount);

	open INF,"../data.2016_co_nco_data.tsv" or die "can't open file data.2016_co_nco_data.tsv: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split/\t/, $_;
		next unless $F[0] =~ /$chr/;
		next unless $F[3] eq "sco";
	
		$scoCount{$F[0]}++;
		$scoList{$F[0]}{$scoCount{$F[0]}} = $F[2];
	}
	close INF;

	# Bootstrap DCOs by selecting two SCOs from WT and determine distance between
	my($totalGap,%allDCOData);
	my $tmpTrial = $trials;
	$tmpTrial *= 2 if $chr eq "chr2";
	foreach my $trial (1..$tmpTrial) {
		my $tmpChr = $chr;
		if ($chr =~ /chr2/) {
			my $chr2 = int(rand(2)) + 1; # value of either 1 or 2
			$tmpChr = "chr2L";
	   	   	$tmpChr = "chr2R" if $chr2 == 2;
		}

		my $position1 = int(rand($scoCount{$tmpChr})) + 1;
		my $position2 = int(rand($scoCount{$tmpChr})) + 1;
		next if $position1 == $position2;

		my $gap = abs($scoList{$tmpChr}{$position1} - $scoList{$tmpChr}{$position2});

		$bootstrapOutput .= "$tmpChr\t$gap\n";
		$allDCOData{$gap}++;
		$totalGap += $gap;

		#print "$scoList{$position1} - $scoList{$position2}\t$gap\n";
	}

	my %obsDCOGap;
   	$obsDCOGap{"chrX"}{"IC"} = 8478101; # chrX DCO ave in IC paper - counts TCOs and QCOs as DCOs
   	$obsDCOGap{"chrX"}{"wt"} = 8344960; # chrX DCO in 2016 paper
   	$obsDCOGap{"chr2L"}{"IC"} = 9811734; # chr2L DCO in IC paper
   	$obsDCOGap{"chr2L"}{"wt"} = 11792279; # chr2L DCO in 2016 paper
   	$obsDCOGap{"chr2R"}{"IC"} = 8448263; # chr2R DCO in IC paper - counts TCO as 2 DCOs
   	$obsDCOGap{"chr2R"}{"wt"} = 10850377; # chr2R DCO in 2016 paper

	my $low_2_5 = $trials * 0.025; #100000 * 2.5%
	my $high_2_5 =  $trials * 0.975; #100000 * 97.5%
	my($bottom_2_5,$top_2_5,$gapCount);

	foreach my $gap (sort {$a<=>$b} keys %allDCOData) {
		$gapCount += $allDCOData{$gap};

		$bottom_2_5 = $gap if $gapCount >= $low_2_5 && !$bottom_2_5;
		#print "$gap\t$gapCount\n" if $gapCount >= $high_2_5 && !$top_2_5;
		$top_2_5    = $gap if $gapCount >= $high_2_5 && !$top_2_5;
	}

	my $averageDist = sprintf("%0.0f",$totalGap / $trials);
	my $ninteyFiveCI = "$bottom_2_5 - $top_2_5";


	print "$chr Average: $averageDist\t95% CI: $ninteyFiveCI\n";
	print "Count of events greater than observed DCO gap for:\n";

	my @tmpChrs = qw/ chrX /;
	   @tmpChrs = qw/ chr2L chr2R / if $chr =~ /chr2/;
	foreach my $tmpChr (@tmpChrs) {
		foreach my $class ("wt","IC") {
			my $countGreaterThanDCO;
			my $count;
			foreach my $gap (sort {$a<=>$b} keys %allDCOData) {
				$count+= $allDCOData{$gap};
				$countGreaterThanDCO += $allDCOData{$gap} if $gap > $obsDCOGap{$tmpChr}{$class};
			}

			my $pval = sprintf("%0.3f",$countGreaterThanDCO / $trials);
			print "\t$tmpChr\t$class\t$countGreaterThanDCO\tp = $pval\t$count\n";
		}
	}
	print "\n";

	next;

	# Sample DCOs with replacement
	my(%dcoList,$dcoCount);
	open INF,"data.ICE_dco.tsv";
	while (<INF>) {
		chomp($_);
		next if $_ =~ /^#/;

		my(@F) = split/\t/, $_;
		next unless $F[1] =~ /$chr/;
		my $gap = abs($F[2] - $F[3]);

		$dcoCount++;
		$dcoList{$dcoCount} = $gap; # store the size 
	}
	print "DCOs that will be bootstrapped: $dcoCount\n";

	my(%allDCOs,$allTrialTotal);
	foreach my $trial (1..$trials) {
		my $cycleCount;
		foreach my $event (1..$dcoCount) {
			my $dcoChosen = int(rand($dcoCount)) + 1;

			#print "$event gap: $dcoList{$dcoChosen}\n";
			$cycleCount += $dcoList{$dcoChosen};
		}

		my $trialAve = sprintf("%0.0f",$cycleCount / $dcoCount);
		$allDCOs{$trialAve}++;
		$allTrialTotal += $trialAve;
	}

	my($bottom_2_5,$top_2_5,$gapCount);
	foreach my $gap (sort {$a<=>$b} keys %allDCOs) {
		$gapCount += $allDCOs{$gap};

		$bottom_2_5 = $gap if $gapCount >= $low_2_5 && !$bottom_2_5;
		#print "$gap\t$gapCount\n" if $gapCount >= $high_2_5 && !$top_2_5;
		$top_2_5    = $gap if $gapCount >= $high_2_5 && !$top_2_5;
	}

	my $averageDist = sprintf("%0.0f",$allTrialTotal / $trials);
	my $ninteyFiveCI = "$bottom_2_5 - $top_2_5";

	print "$chr Average: $averageDist\t95% CI: $ninteyFiveCI\n";
}

open OUTF,">out_bootstrapDCOs.tsv";
print OUTF $bootstrapOutput;
close OUTF;
