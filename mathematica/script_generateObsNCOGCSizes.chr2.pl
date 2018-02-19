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
foreach my $currChr ("chr2L","chr2R") {
	my(%snpCount,%lastSNP);
	my $vcfCount = 52; # we looked at 52 2nd chromosomes in FM7; +/+; TM6B females
	open INF,"../out_uniqueParentalVariants.chr2.tsv" or die "Can't open out_uniqueParentalVariants.chr2.tsv: $!";
	while (<INF>) {
        	my(@F) = split /\t/, $_;
        	next unless $F[1] =~ /[0-9]/;

        	my($chr,$id,$ref,$acns,$avcf,$bcns,$bvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

        	next if $avcf ne $acns && $avcf !~ /\./;
        	next if $bvcf ne $bcns && $bvcf !~ /\./;
        	next if $acns !~ /[A|G|C|T]/ || $bcns !~ /[A|G|C|T]/;
        	next if $acns eq $bcns;

		next unless $chr eq $currChr;

		$snpCount{$id}++
	}

	print "{";
	open INF,"data.chr2NCOGCs.tsv" or die "Can't open file: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[2] eq $currChr;
		my $lowParent = $F[3];
		my $highParent = $F[6];

		my($firstNCOSNP,$lastNCOSNP);
		foreach my $id (sort {$a<=>$b} keys %snpCount) {
			$firstNCOSNP = $id if !$firstNCOSNP && $id > $lowParent;
			$lastNCOSNP = $id if $id < $highParent;
		}
print "$F[2]\t$lowParent\t$highParent\t$firstNCOSNP\t$lastNCOSNP\n";
		my $maxGap = $highParent - $lowParent;
		my $leftGap = $firstNCOSNP - $lowParent;
		my $rightGap = $highParent - $lastNCOSNP;
		#next if $leftGap > 1000;
		#next if $rightGap > 1000;
		print "{$maxGap, $leftGap, $rightGap}, " if $leftGap < $rightGap;
		print "{$maxGap, $rightGap, $leftGap}, " if $rightGap < $leftGap;
		++$count;
	}

}

print "Count: $count\n";
