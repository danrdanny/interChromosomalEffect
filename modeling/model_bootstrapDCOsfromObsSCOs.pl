#!/usr/bin/perl

use strict;
use Getopt::Std;

# bootstrap DCOs from SCO events, use data from 2016
# Is the SCO data from the interchromosomal effect differnt from wt?
# chrX only

my %opts;
getopts('c:h', \%opts); 

## Usage Output
if ($opts{'h'} || !$opts{'c'}) {
	print "
	Required:
		-c Chromosome (chrX, or chr2)
	
	Optional:
		-h This helpful help
	\n";
	exit 0;
}


my $trials = 100000;

my $chr = $opts{'c'};

my(%scoListchrX,%scoListchr2,$scoCount);
open INF,"data.2016_co_nco_data.tsv";
while (<INF>) {
	chomp($_);
	my(@F) = split/\t/, $_;
	next unless $F[0] =~ /$chr/;
	next unless $F[3] eq "sco";

	$scoCount++;
	$scoListchrX{$scoCount} = $F[2] if $chr =~ /chrX/; # store the position
	$scoListchr2{$F[0]}{$scoCount} = $F[2] if $chr =~ /chr2/; 
}
close INF;

# Bootstrap DCOs by selecting two SCOs from WT and determine distance between
my($outData,$totalGap,%allDCOData);
foreach my $trial (1..$trials) {
	my $position1 = int(rand($scoCount)) + 1;
	my $position2 = int(rand($scoCount)) + 1;
	next if $position1 == $position2;

	my $chr2 = int(rand(2)) + 1; # value of either 1 or 2
	my $tmpChr = "chr2L";
	   $tmpChr = "chr2R" if $chr2 == 2;

	my $gap;
	$gap = abs($scoListchrX{$position1} - $scoListchrX{$position2}) if $chr =~ /chrX/;
	$gap = abs($scoListchr2{$tmpChr}{$position1} - $scoListchr2{$tmpChr}{$position2}) if $chr =~ /chr2/;

	$allDCOData{$gap}++;
	$totalGap += $gap;
	$outData .= "$trial\t$gap\n";

	#print "$scoList{$position1} - $scoList{$position2}\t$gap\n";
}

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

