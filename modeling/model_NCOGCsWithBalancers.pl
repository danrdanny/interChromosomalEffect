#!/usr/bin/perl

use strict;

# We need to estimate 3 events:
# 1. NCOGCs onto CyO
# 2. NCOGCs onto TM6B
# 3. NCOGCs onto non-TM6B

my $trials = 100000;

# Assumptions
my $dsbmax = 24;
my $dsbmin = 16;
my $NCOGCSize = 450; # From Miller et al., 2016 Genetics

# Expected based on 2016 data # unnecessary for NCOGCs onto balancers
#my $COPerArm = 0.55;  # based on 541 COs / 980 arms
#my $NCOPerArm = 0.56; # based on 549 NCO-gcs / 980 arms

my %vcfLoc;
   $vcfLoc{"CyO"}      = "../vcf-merge.CyO.vcf.gz"; # merge of all CyO/ISO-1
   $vcfLoc{"TM6B"}     = "../vcf-merge.TM6B.vcf.gz"; # merge of all TM6B/ISO-1
   $vcfLoc{"non-TM6B"} = "../vcf-merge.nonTM6B.vcf.gz"; # merge of all non-TM6B/ISO-1

my %activeChrs;
   $activeChrs{"chr2L"} = "CyO";
   $activeChrs{"chr2R"} = "CyO";
   $activeChrs{"chr3L"} = "TM6B";
   $activeChrs{"chr3R"} = "TM6B";

my %indvCount;
   $indvCount{'CyO'}{"chr2L"} = 49;
   $indvCount{'CyO'}{"chr2R"} = 49;
   $indvCount{'TM6B'}{"chr3L"} = 57;
   $indvCount{'TM6B'}{"chr3R"} = 57;
   $indvCount{'non-TM6B'}{"chr3L"} = 42;
   $indvCount{'non-TM6B'}{"chr3R"} = 42;

my %chr;
   $chr{'CyO'}[0] = "chr2L";
   $chr{'CyO'}[1] = "chr2R";
   $chr{'TM6B'}[0] = "chr3L";
   $chr{'TM6B'}[1] = "chr3R";
   $chr{'non-TM6B'}[0] = "chr3L";
   $chr{'non-TM6B'}[1] = "chr3R";

my %chrs;
   $chrs{1} = "chrX";
   $chrs{2} = "chr2L";
   $chrs{3} = "chr2R";
   $chrs{4} = "chr3L";
   $chrs{5} = "chr3R";
   $chrs{6} = "chr4";

my %chromosomeSizes;
   $chromosomeSizes{"chrX"} = 23542271;
   $chromosomeSizes{"chr3R"} = 32079331;
   $chromosomeSizes{"chr3L"} = 28110227;
   $chromosomeSizes{"chr2L"} = 23513712;
   $chromosomeSizes{"chr2R"} = 25286936;
   $chromosomeSizes{"chr4"} = 1348131;

my %chromosomeFraction;
my $genomeSize = 23542271 + 32079331 + 28110227 + 23513712 + 25286936 + 1348131;
my($low,$newHigh) = (0,0);
foreach my $chr (keys %chromosomeSizes) {
	my $fraction = sprintf("%0.3f", $chromosomeSizes{$chr} / $genomeSize);
	$fraction *= 1000;
	$newHigh += $fraction;
	foreach ($low..$newHigh) {
		$chromosomeFraction{$_} = $chr;
	}
	#print "$chr\t$low - $newHigh\n";
	$low = $newHigh + 1;
}

foreach my $class ("CyO","TM6B","non-TM6B") {
	my(%snpPositions,%count);
	print "$class\t";

	foreach my $currChr (@{$chr{$class}}) {
		open INF,"gunzip -c $vcfLoc{$class} |" or die "Can't open $vcfLoc{$class}: $!";
		while (<INF>) {
			my(@F) = split /\t/, $_;
			my($chr,$id,$vcfScore,$ref,$alt,$det) = ($F[0],$F[1],$F[5],$F[3],$F[4],$F[7]);
			next unless $chr eq $currChr;
			next unless $vcfScore >= 220;
			next unless $ref =~ /^\w$/;
			next unless $alt =~ /^\w$/;
			next if $_ =~ /INDEL/;
			next if $alt =~ /\,/;

			my($mergeCount) = $det =~ /^AC\=(\d+)\;/;
			#print "$mergeCount " if $mergeCount != $indvCount{$class}{$chr};
			next unless $mergeCount >= $indvCount{$class}{$chr};

			my $type = "het";
			   $type = "hom" if $F[9] =~ /1\/1/;

			next if $type =~ /hom/;

			$snpPositions{$chr}{$id} = 1;
			$count{$chr}++;
		}
		print "$currChr SNPs: $count{$currChr}\t";
	}
	print "\n";

	my(%visibleNCO,%nonvisibleNCO,%NCOGCs);
	foreach my $trial (1..$trials) {
		#print "$trial\n" if $trial =~ /000/;

		my %currData;
		my $dsbs = int(rand($dsbmax - $dsbmin)) + $dsbmin;

		# Distribute DSBs randomly
		foreach my $dsbNum (1..$dsbs) {
			## Will this break be repaired as a NCO or CO?
			 # Assume 5-6 DSBs repaired as COs
			my $isCO = int(rand($dsbs)) + 1;
			my $COcount = int(rand(2)) + 5;
			#next if $isCO <= $COcount;

			## Which chromatid is this going on? 
			 # we only recover 1 of 4 chromatids with NCOGCs
			my $chromatid = int(rand(4)) + 1;
			next unless $chromatid == 1; 

			my $chrNum = int(rand(1000)) + 1;
			my $chrName = $chromosomeFraction{$chrNum};

			next unless $activeChrs{$chrName} =~ /$class/ || ($chrName =~ /3/ && $class =~ /non-TM6B/);

			my $maxbp = $chromosomeSizes{$chrName};
			my $pos = int(rand($maxbp)) + 1;

			my $halfNCOGCSize = $NCOGCSize / 2;
			my $lowPos = $pos  - $halfNCOGCSize;
			my $highPos = $pos + $halfNCOGCSize;

			my $visible = 0;
			foreach my $loc ($lowPos..$highPos) {
				$visible = 1 if $snpPositions{$chrName}{$loc} == 1;
			}

			$visibleNCO{$chrName}++    if $visible == 1;
			$nonvisibleNCO{$chrName}++ if $visible == 0;
		}
	}

	foreach my $chr (@{$chr{$class}}) {
		my $visibleFreq = $visibleNCO{$chr} / $trials;
		my $nonvisibleFreq = $nonvisibleNCO{$chr} / $trials;

		print "\t$chr\n";
		print "\tVisible NCOs: $visibleNCO{$chr}\t$visibleFreq per individual\n";
		print "\tNon-Visible NCOs: $nonvisibleNCO{$chr}\t$nonvisibleFreq per individual\n";
	}


#	foreach my $stock (1..$stocks) {
#		foreach my $chr (@chrList) {
#			if ($modeledCOs{$stock}{$chr} =~ /\d+\,\d+\,\d+\,\d+\,/) {
#				$coCounts{"QCO"}++; 
#			} elsif ($modeledCOs{$stock}{$chr} =~ /\d+\,\d+\,\d+\,/) {
#				$coCounts{"TCO"}++; 
#			} elsif ($modeledCOs{$stock}{$chr} =~ /\d+\,\d+\,/) {
#				$coCounts{"DCO"}++; 
#			} elsif ($modeledCOs{$stock}{$chr} =~ /\d+\,/) {
#				$coCounts{"SCO"}++; 
#			} else {
#				$coCounts{"NCO"}++; 
#			}
#		}
#	}
#	
#	### Distribute NCO events
#	my $countVisibleNCOs;
#	foreach my $ncoNum (1..$totalNCOs) {
#		my $stock = int(rand($stocks)) + 1;
#		my $chrNum = int(rand($armsPerStock)) + 1;
#		my $chrName = $chrs{$chrNum};
#		my $maxbp = $chromosomeSizes{$chrName};
#		my $pos = int(rand($maxbp)) + 1;
#
#		#print "$stock\t$chrName\t$pos\n";
#
#		# Is the NCO visible?
#		my $span = int(rand(1000))+450;
#
#		my $ncoStart = $pos - $span;
#		my $ncoEnd = $pos + $span;
#		my $visible = 0;
#		foreach my $snppos ($ncoStart..$ncoEnd) {
#			$visible = 1 if $snpDensity{$chrName}{$snppos} == 1;
#		}
#
#		$modeledNCOs{$stock}{$chrName} .= "$pos," if $visible == 1;
#		$countVisibleNCOs++ if $visible == 1;
#	}
#
#	foreach my $stock (1..$stocks) {
#		foreach my $chr (@chrList) {
#			if ($modeledNCOs{$stock}{$chr} =~ /\d+\,\d+\,\d+\,\d+\,/) {
#				$ncoCounts{"Q_NCO"}++; 
#			} elsif ($modeledNCOs{$stock}{$chr} =~ /\d+\,\d+\,\d+\,/) {
#				$ncoCounts{"T_NCO"}++; 
#			} elsif ($modeledNCOs{$stock}{$chr} =~ /\d+\,\d+\,/) {
#				$ncoCounts{"D_NCO"}++; 
#			} elsif ($modeledNCOs{$stock}{$chr} =~ /\d+\,/) {
#				$ncoCounts{"S_NCO"}++; 
#			} else {
#				$ncoCounts{"N_NCO"}++; 
#			}
#		}
#	}
#	$ncoCounts{"visibleCount"}{$countVisibleNCOs}++;
}

#open OUTF,">./$outputFile";
#print OUTF $output;
#close OUTF;

#print "Crossover Counts:\n";
#foreach my $event ("NCO","SCO","DCO","TCO","QCO") {
#	my $ave = sprintf("%0.1f", $coCounts{$event} / $trials);
#	print "$event\t$ave\n";
#}

#print "\n";
#print "Non-crossover Counts:\n";
#foreach my $event ("N_NCO","S_NCO","D_NCO","T_NCO","Q_NCO") {
#	my $ave = sprintf("%0.1f", $ncoCounts{$event} / $trials);
#	print "$event\t$ave\n";
#}

#print "\n";
#print "Visible NCOs per trial\n";
#foreach my $visible (sort keys %{$ncoCounts{"visibleCount"}}) {
#	print "$visible\t$ncoCounts{\"visibleCount\"}{$visible}\n";
#}
