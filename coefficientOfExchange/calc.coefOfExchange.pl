#!/usr/bin/perl

use strict;

my %chr;
   $chr{'chrX'} = 23542271;
   $chr{'chr2L'} = 23513712;
   $chr{'chr2R'} = 25286936;
   #$chr{'chr3L'} = 27; #28110227;
   #$chr{'chr3R'} = 33; #32079331;
   #$chr{'chr4'} = 1348131;

my %individuals;
   $individuals{'chrX'} = 50;
   $individuals{'chr2L'} = 52;
   $individuals{'chr2R'} = 52;
   $individuals{'chrX'} = 196;
   $individuals{'chr2L'} = 196;
   $individuals{'chr2R'} = 196;

my $output .= "Chr\tPosition\tCOEF\n";

my(%coPos,%totalCO,%coEvents);
#foreach my $file ("out_coData.tsv","out_dcoData.tsv") {
foreach my $file ("cogcDetail.r6.wt.tsv") {
	#print "$file\n";
	open INF,"../$file" or die "can't open $file: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;

		next unless $F[2] > 0;
		next if $F[3] =~ /nco/i;

		my $chr = $F[0]; ## for wt data
		#my $chr = $F[1]; ## for ice data
		my $id = $F[2];
		$id =~ s/\.5//;

		$coEvents{$chr}{$id} = 1;

		if ($F[3] > 0) {
			$F[3] =~ s/\.5//;
			$coEvents{$chr}{$F[3]} = 1;
		}

		$totalCO{$chr}++;
		#$totalCO{$chr} += 1 if $_ =~ /sco/;
		#$totalCO{$chr} += 2 if $_ =~ /dco/;
		#$totalCO{$chr} += 3 if $_ =~ /tco/;
		#$totalCO{$chr} += 4 if $_ =~ /qco/;
	}
}
close INF;

foreach my $chr (keys %chr) {
	my $min = 0;
	my $max = $chr{$chr};
	my $step = 500000;
	my $currMax = $step * 2;

	while ($currMax < $max) {
		foreach my $id (keys %{$coEvents{$chr}}) {
#print "$chr\t$id\n";
			$coPos{$chr}{$currMax}++ if $id > $min && $id <= $currMax;
			#print "$chr\t$id\tbetween $min - $currMax\n" if $id > $min && $id <= $currMax;
		}

		$min += $step;
		$currMax += $step;
	}
}

my $output .= "Chr\tPosition\tCOEF\n";
print "Chr\tPosition\tCOEF\n";
foreach my $chr (keys %chr) { 
	my $min = 0;
	my $max = $chr{$chr};
	my $step = 500000;
	my $currMax = $step * 2;
	while ($currMax < $max) {
		my $count = $coPos{$chr}{$currMax};
		$count = 0 if !$count;

		my $coef = 0;
		   #$coef = sprintf("%0.2f", $count / $totalCO{$chr}) if $totalCO{$chr} > 0;
		   $coef = sprintf("%0.2f", $count / $individuals{$chr}) * 100; # if $totalCO{$chr} > 0;
	
		$output .= "$chr\t$currMax\t$coef\n";
		print "$chr\t$currMax\t$coef\n";
		$min += $step;
		$currMax += $step;
	}
}

#open OUTF,">coefficentOfExchange.coOnly.wt.dat";
#print OUTF $output;
#close OUTF;
