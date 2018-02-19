#!/usr/bin/perl

use strict;

my %chromosomeSizes;
$chromosomeSizes{"chr3R"} = 32079331;
$chromosomeSizes{"chr3L"} = 28110227;
$chromosomeSizes{"chr2L"} = 23513712;
$chromosomeSizes{"chrX"} = 23542271;
$chromosomeSizes{"chr2R"} = 25286936;
$chromosomeSizes{"chr4"} = 1348131;

my %finalOutput;
foreach my $currChr ("chrX") {

	my $vcfCount = 50; # we looked at 50 X chromosomes in +/+; CyO/+; TM6B females
	my $lastSNP;
	open INF,"../out_uniqueParentalVariants.chrX.tsv" or die "Can't open out_uniqueParentalVariants.chrX.tsv: $!";
	while (<INF>) {
        	my(@F) = split /\t/, $_;
        	next unless $F[1] =~ /[0-9]/;

        	my($chr,$id,$ref,$acns,$avcf,$bcns,$bvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

        	next if $avcf ne $acns && $avcf !~ /\./;
        	next if $bvcf ne $bcns && $bvcf !~ /\./;
        	next if $acns !~ /[A|G|C|T]/ || $bcns !~ /[A|G|C|T]/;
        	next if $acns eq $bcns;

		next unless $chr eq $currChr;
	
		if ($lastSNP) {
			my $gap = $id - $lastSNP;
			$finalOutput{$gap} += $vcfCount if $gap <= 10000;
		}
		$lastSNP = $id;		
	}
}

print "{";
foreach my $gap (1..10000) {
	my $value = $finalOutput{$gap};
	   $value = 0 if !$value;
	print "{$gap, $value}, ";
}

print "}\n";
