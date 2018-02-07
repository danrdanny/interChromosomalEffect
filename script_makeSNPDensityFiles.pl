#!/usr/bin/perl

use strict;

my %snpDensity;

my %chrToCheck;
   $chrToCheck{"chrX"} = "out_";
   $chrToCheck{"chr2L"} = "out_";
   $chrToCheck{"chr2R"} = "out_";

#my @chrList = qw/ chr2L chr2R /;
my @chrList = qw/ chrX /;

my %chromosomeSizes;
$chromosomeSizes{"chrX"} = 23542271;
$chromosomeSizes{"chr3R"} = 32079331;
$chromosomeSizes{"chr3L"} = 28110227;
$chromosomeSizes{"chr2L"} = 23513712;
$chromosomeSizes{"chr2R"} = 25286936;
$chromosomeSizes{"chr4"} = 1348131;

#open INF,"out_uniqueOregonRVariants.tsv";
open INF,"out_uniquew1118Variants.tsv";
#open INF,"out_uniqueParentalVariants.tsv";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;

	my($chr,$id,$ref,$wcns,$wvcf,$cscns,$csvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

	next if $wvcf ne $wcns && $wvcf !~ /\./;
	next if $csvcf ne $cscns && $csvcf !~ /\./;
	next if $wcns !~ /[A|G|C|T]/ || $cscns !~ /[A|G|C|T]/;
	next if $wcns eq $cscns;

	$snpDensity{$chr}{$id} = 1;
}
close INF;

my $output = "Chr\tPos\tCount\n";
foreach my $chr (@chrList) {
	my $min = 0;
	my $step = 20000;
	my $curr = $step;

	while ($curr < $chromosomeSizes{$chr}) {
		my $SNPCount = 0;
		foreach my $currPos ($min..$curr) {
			$SNPCount++ if $snpDensity{$chr}{$currPos} == 1;
		}
		$SNPCount = sprintf("%0.5f", $SNPCount / $step);
		$output .= "$chr\t$curr\t$SNPCount\n";
		print "$chr\t$curr\t$SNPCount\n";

		$min = $curr + 1;
		$curr += $step;
	}
}

#open OUTF,">out_parentalSNPDensity.tsv";
#open OUTF,">out_parentalSNPDensity.OregonR.tsv";
open OUTF,">out_parentalSNPDensity.w1118tsv";
print OUTF $output;
close OUTF;
