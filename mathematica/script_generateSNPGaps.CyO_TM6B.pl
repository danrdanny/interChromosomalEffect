#!/usr/bin/perl

use strict;
use Getopt::Std;

my %opts;
getopts('c:h', \%opts); 

if ($opts{'h'} || !$opts{'c'}) {
	print "
	Required:
		-c Chromosome you want to get data for, a for all

	Optional:
		-h this helpful help

	\n";
	exit 0;
}

my %chromosomeSizes;
$chromosomeSizes{"chr3R"} = 32079331;
$chromosomeSizes{"chr3L"} = 28110227;
$chromosomeSizes{"chr2L"} = 23513712;
$chromosomeSizes{"chrX"} = 23542271;
$chromosomeSizes{"chr2R"} = 25286936;
$chromosomeSizes{"chr4"} = 1348131;

my %finalOutput;
foreach my $balancer ("CyO","TM6B") {
	next unless $balancer eq "TM6B";
	my($mergeFile,@chrs,$chrPrefix,$vcfCount);
	if ($balancer =~ /CyO/) {
		@chrs = qw/ chr2L chr2R /;
		@chrs = qw/ chr2R /;
		$chrPrefix = "chr2R";
		$mergeFile = "vcf-merge.CyO.vcf.gz";
		$vcfCount = 50; 
	} elsif ($balancer =~ /TM6B/) {
		@chrs = qw/ chr3L chr3R /;
		@chrs = qw/ chr3R /;
		$chrPrefix = "chr3";
		$chrPrefix = "chr3R";
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

	my $lastSNP;
	foreach my $chr (@chrs) {
		foreach my $id (sort {$a<=>$b} keys %{$snpCount{$chr}}) {
			#print "$id\t$snpCount{$chr}{$id}\n";
	
			if ($lastSNP) {
				my $gap = $id - $lastSNP;
				$finalOutput{$gap} += $vcfCount if $gap <= 10000;
			}
			$lastSNP = $id;		
		}
	}
}

print "{";
foreach my $gap (1..10000) {
	my $value = $finalOutput{$gap};
	   $value = 0 if !$value;
	print "{$gap, $value}, ";
}

print "}\n";
