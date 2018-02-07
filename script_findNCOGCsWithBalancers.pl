#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options
my %opts;
getopts('t:h', \%opts); # options as above. Values in %opts

## Usage Output
if ($opts{'h'} || !$opts{'t'}) {
        print "

	Required:
		-t What you're making a heatmap of, options are:

			CyO
			3rd_TM6B
			3rd_notTM6B

	Optional:

		-h This helpful help.

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

my(@stocks,%stocks,$outputFile,$description,@chrs,$chrPrefix);
my $output = "stock,chr,pos,percent\n";

if ($opts{'t'} =~ /CyO/) {
	$outputFile = "out_heatmap_2nd_CyO.tsv";
	@chrs = qw/ chr2L chr2R /;

	@stocks = qw/ chrX-01-M chrX-02-M chrX-03-M chrX-04-M chrX-05-M chrX-06-M chrX-07-M chrX-08-M chrX-09-M chrX-10-M chrX-11-M chrX-12-M chrX-13-M chrX-14-M chrX-15-M chrX-16-M chrX-17-M chrX-18-M chrX-19-M chrX-20-M chrX-21-M chrX-22-M chrX-23-M chrX-24-M chrX-25-M chrX-26-M chrX-27-M chrX-28-M chrX-29-M chrX-30-M chrX-31-M chrX-32-M chrX-33-M chrX-34-M chrX-35-M chrX-36-M chrX-37-M chrX-38-M chrX-39-M chrX-40-M chrX-41-M chrX-42-M chrX-43-M chrX-44-M chrX-45-M chrX-46-M chrX-47-M chrX-48-M chrX-49-M chrX-50-M /;

} elsif ($opts{'t'} =~ /3rd_TM6B/) {
	$outputFile = "out_heatmap_3rd_TM6B.tsv";
	@chrs = qw/ chr3L chr3R /;
	$chrPrefix = "chr3";
	$description = "Unique SNPs found on third chr balancer TM6B";

	@stocks = qw/ chr2-03-M-_wt_Y___-recomb_ISO1-__wt_TM6B chr2-13-F-_wt_ISO1-recomb_ISO1-__wt_ISO1 chr2-15-F-_wt_ISO1-recomb_ISO1-TM6B_ISO1 chr2-16-M-_wt_Y___-recomb_ISO1-__wt_ISO1 chr2-19-F-_wt_ISO1-recomb_ISO1-TM6B_ISO1 chr2-21-F-FM7_ISO1-recomb_ISO1-TM6B_ISO1 chr2-25-F-FM7_ISO1-recomb_ISO1-TM6B_ISO1 chr2-26-F-_wt_ISO1-recomb_ISO1-TM6B_ISO1 chrX-01-M chrX-02-M chrX-03-M chrX-04-M chrX-05-M chrX-06-M chrX-07-M chrX-08-M chrX-09-M chrX-10-M chrX-11-M chrX-12-M chrX-13-M chrX-14-M chrX-15-M chrX-16-M chrX-17-M chrX-18-M chrX-19-M chrX-20-M chrX-21-M chrX-22-M chrX-23-M chrX-24-M chrX-25-M chrX-26-M chrX-27-M chrX-28-M chrX-29-M chrX-30-M chrX-31-M chrX-32-M chrX-33-M chrX-34-M chrX-35-M chrX-36-M chrX-37-M chrX-38-M chrX-39-M chrX-40-M chrX-41-M chrX-42-M chrX-43-M chrX-44-M chrX-45-M chrX-46-M chrX-47-M chrX-48-M chrX-49-M chrX-50-M /;

} elsif ($opts{'t'} =~ /3rd_notTM6B/) {
	$outputFile = "out_heatmap_3rd_notTM6B.tsv";
	@chrs = qw/ chr3L chr3R /;
	$chrPrefix = "chr3";
	$description = "Unique SNPs found on non-TM6B 3rd chromosomes";
	@stocks = qw/ chr2-01-M-_wt_Y___-recomb_ISO1-__wt_ISO1 chr2-04-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-05-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-06-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-07-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-08-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-10-F-_wt_ISO1-recomb_ISO1-__wt_ISO1 chr2-12-F-FM7_ISO1-recomb_ISO1-TM6B_ISO1 chr2-18-M-FM7_Y___-recomb_ISO1-__wt_ISO1 chr2-24-M-_wt_Y___-recomb_ISO1-__wt_ISO1 chr2-27-M-_wt_Y___-recomb_ISO1-__wt_ISO1 chr2-29-F-FM7_ISO1-recomb_ISO1-__wt_ISO1 chr2-35-M chr2-36-M chr2-37-M chr2-38-M chr2-39-M chr2-40-M chr2-41-M chr2-42-M chr2-43-M chr2-44-M chr2-45-M chr2-46-M chr2-47-M chr2-50-M chr2-51-M chr2-52-M chr2-55-M chr2-56-M chr2-57-M chr2-58-M chr2-59-M chr2-60-M chr2-61-M chr2-62-M chr2-63-M chr2-64-M chr2-65-M chr2-66-M chr2-67-M chr2-68-M chr2-69-M /;

}

foreach (@stocks) {
	chomp($_);
	$_ =~ s/\///;
	$stocks{$_} = 1;
}


if ($opts{'t'} =~ /CyO/ || $opts{'t'} =~ /3rd_notTM6B/ || $opts{'t'} =~ /3rd_TM6B/) {
	foreach my $chr (@chrs) {
		my(%vcfData,%allSNPs,%vcfScores,%max,%snpCount,$stockCount,%cnsData);

		foreach my $stock (sort keys %stocks) {
			print "opening $stock.vcf.gz\n";
			$stockCount++;
			open INF,"gunzip -c ~/projects/interChrom/$stock/$stock.vcf.gz |" or die "Can't open $stock.vcf.gz: $!";
			foreach my $line (<INF>) {
				my(@F) = split /\t/, $line;
				next unless $F[0] eq $chr;

				#next unless $F[1] > 6858000 && $F[1] < 6858600;

				#next unless $line =~ /0\/1/;
				next unless $F[5] > 50;
				next if $line =~ /INDEL/;
				my($prefix) = $F[1] =~ /^(\d\d)/;
				$prefix = "a" if $F[1] < 1000000;

				$vcfData{$F[0]}{$stock}{$prefix}{$F[1]} = $F[4];
				$vcfScores{$F[0]}{$stock}{$prefix}{$F[1]} = $F[5];
				$allSNPs{$F[0]}{$prefix}{$F[1]}{$F[4]}++;
				$snpCount{$F[0]}{$F[1]}++;
				$max{$F[0]} = $F[1] if $F[1] > $max{$F[0]};
			}

			print "opening $stock.cns.$chr.gz\n";
			open INF,"gunzip -c ~/projects/interChrom/$stock/$stock.cns.$chr.gz |" or die "Can't open $stock.cns.$chr.gz: $!";
			foreach my $line (<INF>) {
				my(@F) = split /\t/, $line;
				next unless $F[0] > 0;
				$F[1] = "N" unless $F[2] > 80;

				$cnsData{$stock}{$chr}{$F[0]} = $F[1];
			}
		}
	

		# find instances where one SNP is missing among all chromosomes
		foreach my $pos (1..$max{$chr}) {
			#next; # comment this out when really running!!!
			next unless $snpCount{$chr}{$pos} > 1;
			#print "$chr\t$pos\t$snpCount{$chr}{$pos} $stockCount\n";
			next unless $snpCount{$chr}{$pos} == $stockCount - 1;

			my($prefix) = $pos =~ /^(\d\d)/;
			$prefix = "a" if $pos < 1000000;

			my $stock;
			foreach my $tmpstock (keys %stocks) {
				next if $vcfData{$chr}{$tmpstock}{$prefix}{$pos};
				$stock = $tmpstock;
			}

			next unless $cnsData{$stock}{$chr}{$pos} =~ /(A|G|C|T)/;

			#print "SNP missing in $stock at $chr $pos\n";
			$vcfData{$chr}{$stock}{$prefix}{$pos} = "- $cnsData{$stock}{$chr}{$pos}";
			$allSNPs{$chr}{$prefix}{$pos}{$cnsData{$stock}{$chr}{$pos}}++;
			$vcfScores{$chr}{$stock}{$prefix}{$pos} = 220;
		}

		foreach my $stock (sort keys %stocks) {
			print "[".localtime(time)."]      - $stock\n";
			my($stockName) = $stock =~ /(chr\d\-\d+)/;

			my $maxPos = $chromosomeSizes{$chr};

			my $currMin 	= 0;
			my $step 	= 5000; # in kb
			my $currMax 	= 5000;
			while ($currMax < $maxPos) {
				my($prefix) = $currMin =~ /^(\d\d)/;
				$prefix = "a" if $currMin < 1000000;
				my $count = 0;
				foreach my $id (sort {$a<=>$b} keys %{$vcfData{$chr}{$stock}{$prefix}}) {
					#$vcfData{$stock}{$id} = undef if $id < $currMin; 

					next unless $id >= $currMin && $id <= $currMax;

					my $snp = $vcfData{$chr}{$stock}{$prefix}{$id};
					next unless $allSNPs{$chr}{$prefix}{$id}{$snp} <= 1;
					next unless $vcfScores{$chr}{$stock}{$prefix}{$id} > 150;
					$count++;
					print "$stock\t$chr\:$id\t$snp\t$vcfScores{$stock}{$prefix}{$id}\n";
				}
				$count = 1 if $count > 0;
				$output .= "$stockName,$chr,$currMax,$count\n";

				$currMin = $currMax + 1;
				$currMax += $step;
			}
		}
	}
}

open OUTF,">$outputFile";
print OUTF $output;
close OUTF;
