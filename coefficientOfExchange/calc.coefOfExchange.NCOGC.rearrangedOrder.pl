#!/usr/bin/perl

use strict;

my $sample = "recombining";
   $sample = "balancers";
my $output = "type\tsample\tChr\tPosition\tCOEF\n";
my $outputFile = "out_cmMbintervals.NCOGC.balancerOrder.tsv";

my(@chrs,%individuals);

my %chr; # WT numbers
   $chr{'chrX'} = 23542271;
   $chr{'chr2L'} = 23513712;
   $chr{'chr2R'} = 25286936;
   $chr{'chr3L'} = 28110227;
   $chr{'chr3R'} = 32079331;
   #$chr{'chr4'} = 1348131;

my %chr; # when it's a balancer
   $chr{'chr2L'} = 19711702; # CyO
   $chr{'chr2R'} = 29806844; # CyO
   $chr{'chr3L'} = 14938242; # TM6B
   $chr{'chr3R'} = 47728532; # TM6B
   #$chr{'chr4'} = 1348131;


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

foreach my $chr (@chrs) {
	my @intervals;
	my(%NCOGCPos,%TotalNCOGC);

	if ($chr eq "chr2L") {
		push(@intervals,"chr2L:1..2151403");
		push(@intervals,"chr2L:12731221..9805575");
		push(@intervals,"chr2R:14067771..22689962");
		push(@intervals,"chr2R:6012459..1");
		#CEN	
	} elsif ($chr eq "chr2R") {
		push(@intervals,"chr2L:23513712..12731221");
		push(@intervals,"chr2L:2151403..9805575");
		push(@intervals,"chr2R:14067771..6012459");
		push(@intervals,"chr2R:21972072..25286936");
			
	} elsif ($chr eq "chr3L") {
		push(@intervals,"chr3L:1..233563");
		push(@intervals,"chr3R:12227473..10742047");
		push(@intervals,"chr3R:7048580..10742047");
		push(@intervals,"chr3R:8287181..7048580");
		push(@intervals,"chr3R:8287181..1");
		#CEN	
	} elsif ($chr eq "chr3R") {
		push(@intervals,"chr3L:28110227..18693657");
		push(@intervals,"chr3R:22393827..32000000");
		push(@intervals,"chr3R:20180000..12227473");
		push(@intervals,"chr3L:233562..3173046");
		push(@intervals,"chr3L:16308841..3173046");
		push(@intervals,"chr3L:16308841..18693657");
		push(@intervals,"chr3R:22393827..20180000");
		push(@intervals,"chr3R:32000000..32079331");
	}

	my $count = 1;
	my %positionHash;
	foreach my $segment (@intervals) {
		my($origChr,$start,$stop) = $segment =~ /(\w+)\:(\d+)\.\.(\d+)/;
		print "$chr\t$origChr\t$start-$stop\n";

		if ($start < $stop) {
			foreach my $pos ($start..$stop) {
				$positionHash{$origChr}{$pos} = $count;
				++$count;
			}
		} else {
			for (my $pos=$start; $pos >= $stop; $pos--) {
				$positionHash{$origChr}{$pos} = $count;
				++$count;
			}
		}
	}

	print "$chr is $count bp\n";

	##----------------------##
	## get WT and IC NCOGCs ##
	##----------------------##

	foreach my $type ("wt","IC") {
		next unless $type eq "IC";
		#my $file = "../data.2016_co_nco_data.tsv";
	   	#$file = "../data.IC_recombining_nco_data.tsv" if $sample eq "recombining" && $type eq "IC";
		my $file = "../data.IC_balancer_nco_data.tsv"  if $sample eq "balancers"   && $type eq "IC";

		open INF,"$file" or die "can't open $file: $!";
		while (<INF>) {
			chomp($_);
			my(@F) = split /\t/, $_;

			next unless $F[2] > 0;
			next unless $F[3] eq "nco";

			my $currChr = $F[0];
			my $id = $F[2];
			$id =~ s/\.5//;

			# Get the new coordinates
			my $newPosition = $positionHash{$currChr}{$id};
			next if !$newPosition;

			print "$currChr $id is now $chr $newPosition\n";

			$NCOGCPos{$type}{$chr}{$newPosition} = 1;
			$TotalNCOGC{$type}{$chr}++;
		}
		close INF;

		##-----------------##
		## Calculate cM/Mb ##
		##-----------------##
		my %intervalCounts;
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

#sleep 10;

		print "Chr\tPosition\tCOEF\n";
		my $min = 0;
		my $max = $chr{$chr};
		my $step = 500000;
		my $currMax = $step * 2;
		while ($currMax < $max) {
			my $count = $intervalCounts{$type}{$chr}{$currMax};
			$count = 0 if !$count;

			my $coef = 0;
	   		$coef = sprintf("%0.2f", $count / $individuals{$chr}{$type}) * 100; # if $totalCO{$chr} > 0;

			$output .= "$type\t$sample\t$chr\t$currMax\t$coef\n";
			print "$type\t$sample\t$chr\t$currMax\t$coef\t$intervalCounts{$type}{$chr}{$currMax}/$individuals{$chr}{$type}\n";
			$min += $step;
			$currMax += $step;
		}
	}
}

open OUTF,">$outputFile";
print OUTF $output;
close OUTF;
