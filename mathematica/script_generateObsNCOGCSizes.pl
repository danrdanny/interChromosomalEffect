#!/usr/bin/perl

use strict;

my %chromosomeSizes;
$chromosomeSizes{"chr3R"} = 32079331;
$chromosomeSizes{"chr3L"} = 28110227;
$chromosomeSizes{"chr2L"} = 23513712;
$chromosomeSizes{"chrX"} = 23542271;
$chromosomeSizes{"chr2R"} = 25286936;
$chromosomeSizes{"chr4"} = 1348131;

my(%finalOutput,$count);
foreach my $balancer ("CyO","TM6B") {
	my($mergeFile,@chrs,$chrPrefix,$vcfCount);
	if ($balancer =~ /CyO/) {
		@chrs = qw/ chr2L chr2R /;
		$chrPrefix = "chr2";
		$mergeFile = "vcf-merge.CyO.vcf.gz";
		$vcfCount = 50; 
	} elsif ($balancer =~ /TM6B/) {
		@chrs = qw/ chr3L chr3R /;
		$chrPrefix = "chr3";
		$mergeFile = "vcf-merge.TM6B.vcf.gz";
		$vcfCount = 58;
	}

	my(%snpCount,$stockCount);
	print "opening $mergeFile\n";
	open INF,"gunzip -c ../$mergeFile |" or die "Can't open $mergeFile: $!";
	foreach my $line (<INF>) {
		my(@F) = split /\t/, $line;
		next unless $F[0] =~ /$chrPrefix/;

		#next unless $line =~ /0\/1/;
		next unless $F[5] > 50;
		next if $line =~ /INDEL/;
		next unless $line =~ /AC\=$vcfCount\;/; #we merge 50 or 58 VCFs for each

		$snpCount{$F[0]}{$F[1]}++;
	}
	close INF;

	print "{";
	open INF,"../data.NCOGCs_with_balancers.tsv" or die "Can't open file: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[3] =~ /$chrPrefix/;
		
		my($lastSNP,$nextSNP);
		foreach my $id (sort {$a<=>$b} keys %{$snpCount{$F[3]}}) {
			#print "$id\t$snpCount{$chr}{$id}\n";
	
			$lastSNP = $id if $id < $F[4];
			$nextSNP = $id if !$nextSNP && $id > $F[5];
		}

		my $maxGap = $nextSNP - $lastSNP;
		my $leftGap = $F[4] - $lastSNP;
		my $rightGap = $nextSNP - $F[5];
		#next if $leftGap > 1000;
		#next if $rightGap > 1000;
		print "{$maxGap, $leftGap, $rightGap}, " if $leftGap < $rightGap;
		print "{$maxGap, $rightGap, $leftGap}, " if $rightGap < $leftGap;
		++$count;
	}
	close INF;
	print "\n";
}

print "total events: $count\n";
