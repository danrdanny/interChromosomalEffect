#!/usr/bin/perl

use strict;

my $sample = "recombining";
   $sample = "balancers";
my $output = "type\tsample\tChr\tPosition\tCOEF\n";
my $outputFile = "out_cmMbintervals.NCOGC.$sample.tsv";

my %chr;
   $chr{'chrX'} = 23542271;
   $chr{'chr2L'} = 23513712;
   $chr{'chr2R'} = 25286936;
   $chr{'chr3L'} = 28110227;
   $chr{'chr3R'} = 32079331;
   #$chr{'chr4'} = 1348131;

my(%individuals,@chrs);

if ($sample =~ /balancers/) {
   @chrs = qw/ chr2L chr2R chr3L chr3R /;

   $individuals{'chr2L'}{'IC'} = 50;
   $individuals{'chr2R'}{'IC'} = 50;
   $individuals{'chr3L'}{'IC'} = 50;
   $individuals{'chr3R'}{'IC'} = 50;
   $individuals{'chr2L'}{'wt'} = 196;
   $individuals{'chr2R'}{'wt'} = 196;
   $individuals{'chr3L'}{'wt'} = 196;
   $individuals{'chr3R'}{'wt'} = 196;

} elsif ($sample =~ /recombining/) {
   @chrs = qw/ chrX chr2L chr2R /;

   $individuals{'chrX'}{'IC'} = 50;
   $individuals{'chr2L'}{'IC'} = 52;
   $individuals{'chr2R'}{'IC'} = 52;
   $individuals{'chrX'}{'wt'} = 196;
   $individuals{'chr2L'}{'wt'} = 196;
   $individuals{'chr2R'}{'wt'} = 196;
}

##----------------------##
## get WT and IC NCOGCs ##
##----------------------##

foreach my $type ("wt","IC") {
	my(%intervalCounts,%NCOGCPos,%TotalNCOGC);
	my $file = "../data.2016_co_nco_data.tsv";
	   $file = "../data.IC_recombining_nco_data.tsv" if $sample eq "recombining" && $type eq "IC";
	   $file = "../data.IC_balancer_nco_data.tsv"    if $sample eq "balancers"   && $type eq "IC";

	open INF,"$file" or die "can't open $file: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;

		next unless $F[2] > 0;
		next unless $F[3] eq "nco";

		my $chr = $F[0];
		my $id = $F[2];
		$id =~ s/\.5//;

		$NCOGCPos{$type}{$chr}{$id} = 1;
		$TotalNCOGC{$type}{$chr}++;
	}
	close INF;

	##-----------------##
	## Calculate cM/Mb ##
	##-----------------##
	foreach my $chr (@chrs) {
		my $min = 0;
		my $max = $chr{$chr};
		my $step = 500000;
		my $currMax = $step * 2;

		while ($currMax < $max) {
			foreach my $id (keys %{$NCOGCPos{$type}{$chr}}) {
				$intervalCounts{$type}{$chr}{$currMax}++ if $id > $min && $id <= $currMax;
				print "$type\t$chr\t$id\tbetween $min - $currMax $intervalCounts{$type}{$chr}{$currMax}\n" if $id > $min && $id <= $currMax;
			}
			$min += $step;
			$currMax += $step;
		}
	}

	print "Chr\tPosition\tCOEF\n";
	foreach my $chr (@chrs) { 
		my $min = 0;
		my $max = $chr{$chr};
		my $step = 500000;
		my $currMax = $step * 2;
		while ($currMax < $max) {
			my $count = $intervalCounts{$type}{$chr}{$currMax};
			$count = 0 if !$count;

			my $coef = 0;
		   	$coef = sprintf("%0.2f", $count / $individuals{$chr}{$type}) * 100; # if $totalCO{$chr} > 0;

			#if (($chr eq "chrX") || 
			    #($chr eq "chr2L") ||
			    #($chr eq "chr2R" && $currMax > 4500000)  ||
			    #($chr eq "chr3L" && $currMax < 24000000) ||
			    #($chr eq "chr3R" && $currMax > 4500000)) {

				$output .= "$type\t$sample\t$chr\t$currMax\t$coef\n";
				print "$type\t$sample\t$chr\t$currMax\t$coef\t$intervalCounts{$type}{$chr}{$currMax}/$individuals{$chr}{$type}\n";
			#}
			$min += $step;
			$currMax += $step;
		}
	}
}

open OUTF,">$outputFile";
print OUTF $output;
close OUTF;
