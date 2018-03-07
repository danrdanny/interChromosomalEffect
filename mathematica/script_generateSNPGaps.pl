#!/usr/bin/perl

use strict;

foreach my $type ("WT","IC") {
	foreach my $currChr ("chrX","chr2L","chr2R","chr3","all") {
		print "$type\t$currChr\n";
		my(%finalOutput,$lastSNP);

		next if $type eq "IC" and $currChr eq "all";

		my $vcfCount;
		   $vcfCount = 50 if $type eq "IC" and $currChr eq "chrX";
		   $vcfCount = 52 if $type eq "IC" and $currChr =~ /chr2/;
		   $vcfCount = 196 if $type eq "WT" and $currChr eq "chrX";
		   $vcfCount = 196 if $type eq "WT" and $currChr =~ /chr2/;
		   $vcfCount = 196 if $type eq "WT" and $currChr eq "chr3";
		   #$vcfCount = 98 if $type eq "WT" and $currChr eq "all";

		my $dataFile;
		   $dataFile = "../out_uniqueParentalVariants.chrX.tsv" if $type eq "IC" and $currChr eq "chrX";
		   $dataFile = "../out_uniqueParentalVariants.chr2.tsv" if $type eq "IC" and $currChr =~ /chr2/;
		   $dataFile = "../out_uniqueParentalVariants.wt2016.tsv" if $type eq "WT";
		   #$dataFile = "/Dropbox/_Manuscripts/__Completed/CO\ and\ GC\ Paper/resubmissionDataDir/snpDensity/out_uniqueParentalVariants.tsv" if $type eq "WT";

		open INF,"$dataFile" or die "Can't open $dataFile: $!";
		while (<INF>) {
        		my(@F) = split /\t/, $_;
        		next unless $F[1] =~ /[0-9]/;

        		my($chr,$id,$ref,$acns,$avcf,$bcns,$bvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

        		next if $avcf ne $acns && $avcf !~ /\./;
        		next if $bvcf ne $bcns && $bvcf !~ /\./;
        		next if $acns !~ /[A|G|C|T]/ || $bcns !~ /[A|G|C|T]/;
        		next if $acns eq $bcns;

			my $skip = 0;
			$skip = 1 if $chr ne $currChr;
			$skip = 0 if $currChr eq "all";
			$skip = 1 if $chr eq "chr4";
			#next unless $chr =~ /$currChr/;
			next if $chr eq "chr2R" && $id < 5000000;
			next if $chr eq "chr3R" && $id < 5000000;
			next if $skip == 1;
	
			if ($lastSNP) {
				my $gap = $id - $lastSNP;
				#print "$chr\t$id\t$lastSNP\t$gap\t$finalOutput{$gap}\t";
				$finalOutput{$gap} += $vcfCount if $gap <= 10000;
				#$finalOutput{$gap} += $vcfCount if $gap <= 10000 && $chr eq "chrX";
				#print "$finalOutput{$gap}\n";
			}
			$lastSNP = $id;		
		}

		my $output = "{";
		foreach my $gap (1..10000) {
			my $value = $finalOutput{$gap};
	   		$value = 0 if !$value;
			$output .= "{$gap, $value}, ";
		}
		$output =~ s/\,\s$//;
		$output .= "}\n";
		
		open OUTF,">data.SNPGaps.$type.$currChr.tsv";
		print OUTF "$output";
		close OUTF;
	}
}
