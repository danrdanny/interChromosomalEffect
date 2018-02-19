#!/usr/bin/perl

use strict;
use Getopt::Std;

my $minCNSBaseQuality	= 60;
my $minVCFScore		= 150; #only use SNPs with scores greater than or equal to this number
my $maxIndelLength	= 1; #ignore indels greater than or equal to this number
my $skipIndels		= 1; # 1 for yes, 0 for no
my $minParentalDepth 	= 20;
my $minChildDepth	= 8;

my $repeatMasker	= "dm6.rmsk.txt";
my $chromSizes		= "dm6.chrom.sizes";

## Gather system data
my $pwd = `pwd`;
chomp($pwd);
my $refLoc = "/Users/danny/projects/genomes/dmel";
my $parentLoc = "parentalData";

## Command-line options
my %opts;
getopts('a:b:ehols:c:', \%opts); # options as above. Values in %opts

## Usage Output
if ($opts{'h'} || !$opts{'a'} || !$opts{'b'} || !$opts{'c'}) {
        print "
	This script identifies CO and NCO events in child .vcf files from two known
	parental genotypes. Parental files must exist in a parental/ directory,
	chid data must exist in a offspring/ directory.

	Required:

		-a Parent A
		-b Parent B
		-c Chromosome: chrX, chr2L, chr2R, chr3L, chr3R, chr4, or A(ll) Default A(ll).

	Optional:

		-e If out_uniqueParentalVariants.tsv already exists you can use it
		   and not re-make it every time.

		-o If you only want to populate (or re-populate out_uniqueParentalVariants.tsv

		-l Populate or add to the cogcPositionsToSkip.tsv file.

		-s Specific stock you want to check.

		-h This helpful help.
        \n";
        exit 0;
}

exit "No parentA passed!\n" if !$opts{'a'};
exit "No parentB passed!\n" if !$opts{'b'};
exit "No chr passed!\n" if !$opts{'c'};

my $parentA     = $opts{'a'};
my $parentB     = $opts{'b'};

my @chrList = qw/ chrX chr2L chr2R chr3L chr3R chr4 /;
my %chrListToCheck;
$chrListToCheck{"chrX"} = 1;
$chrListToCheck{"chr2L"} = 1;
$chrListToCheck{"chr2R"} = 1;
$chrListToCheck{"chr3L"} = 1;
$chrListToCheck{"chr3R"} = 1;
$chrListToCheck{"chr4"} = 1;

if ($opts{'c'} =~ /chr/) {
	undef @chrList;
	undef %chrListToCheck;
	if ($opts{'c'} =~ /\,/) {
		foreach my $inChr (split /\,/, $opts{'c'}) {
			$inChr =~ s/\,//g;
			next unless $inChr =~ /chr/i;
			chomp($inChr);
			push(@chrList,$inChr);
			$chrListToCheck{$inChr} = 1;
		}
	} else {
		push(@chrList,$opts{'c'});
		$chrListToCheck{$opts{'c'}} = 1;
	}
}

print "Checking chromosomes: ";
foreach my $chr (@chrList) { 
	print "$chr ";
}
print "\n";

## Subroutine to check memory usage - Doesn't work on every system
sub memUsage {
	return "unknown";
        my $gb;
        if (-e "/proc/$$/status") {
                open DATA, "< /proc/$$/status" or die "Unable to read /proc/$$/status: $!\n";
                local $/;
                <DATA> =~ m/^VmSize:\s+(.*?)$/m;
                my($kb) = $1 =~ /(\d+)/;
                $gb = sprintf("%0.1f", $kb / 1000 / 1000);
                $gb .= " Gigabytes";
        } else {
                $gb = "unable to check on this OS\n";
        }
        return $gb;
}

my %qualityScores;
foreach (0..90) {
        $_ += 33;
        my $ascii = chr($_);
        $_ -= 33;
        $qualityScores{$ascii} = $_;
}

## Print out all variables for the user
print "[".localtime(time)."] \n";
print "[".localtime(time)."]       Script begin. Variables:\n";
printf "[".localtime(time)."] %30s %8d\n", "minimum cns base quality", $minCNSBaseQuality; 
printf "[".localtime(time)."] %30s %8d\n", "minimum VCF score", $minVCFScore; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for parents", $minParentalDepth; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for child", $minChildDepth; 
if ($skipIndels == 0) {
	printf "[".localtime(time)."] %30s %8d\n", "maximum INDEL length", $maxIndelLength; 
} else {
	printf "[".localtime(time)."] %30s\n", "Skipping all indels"; 
}
print "[".localtime(time)."] \n";

## Grab chromosome names and sizes from the chrom.sizes file
my %chromosomeSizes;
open INF, "$refLoc/$chromSizes" or die "Can't open $chromSizes: $!";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	# format is: <chromosome>\t<size in nt>
	$chromosomeSizes{$F[0]} = $F[1];
}
close INF;

if (!$opts{'e'}) {  # -e flag is to use existing out_uniqueParentalVariants.tsv file
	## Open repeatmasker file
	my %repeats;
	print "[".localtime(time)."] Getting repeats from $repeatMasker.\n";
	open INF,"$refLoc/$repeatMasker";
	while (<INF>) {
		next if $_ =~ /Simple_repeat/;
		my(@F) = split /\t/, $_;
		foreach my $id ($F[6]..$F[7]) {
			$repeats{$F[5]}{$id} = 1;
			my $chr = $F[5];
			#$chr =~ s/chr//;
			$repeats{$chr}{$id} = 1;
		}
	}
	close INF;
	print "[".localtime(time)."] Repeats gathered.\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";

	## Grab parental VCF data
	#my $parentCount = 0;
	my %parentalSNPs;
	my %variantPositions;
	print "[".localtime(time)."] Getting parental data.\n";
	foreach my $parent ($parentA, $parentB) {
		my($SNPcount,$skipRepeat,$skipScore,$skipHet,$skipAlt,$skipINDEL,$skipDepth,$indelCount) = (0,0,0,0,0,0,0,0);
		print "[".localtime(time)."] - Getting vcf data from $parent.\n";
		open INF,"gunzip -c $parentLoc/$parent.vcf.gz |" or die "Can't open $parentLoc/$parent.vcf.gz: $!";
		while (<INF>) {
			my(@F) = split /\t/, $_;
			my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);

			next unless $chrListToCheck{$chr} == 1;
			$skipRepeat++ if $repeats{$chr}{$id} == 1;	# skip if in repetative region
			next if $repeats{$chr}{$id} == 1;	
			$skipScore++ if $vcfScore < $minVCFScore;	# skip based on vcf score
			next unless $vcfScore >= $minVCFScore;	
                	$skipHet++ if $F[9] =~ /0\/1/; 			# skip it if it's het when looking at parents
                	next if $F[9] =~ /0\/1/; 	
                	$skipAlt++ if $alt =~ /\,/; 			# skip if the alt allele has two bases
                	next if $alt =~ /\,/; 		

			my $isIndel = 0;
		   	$isIndel = 1 if $_ =~ /INDEL/;
			next if $isIndel == 1 && $skipIndels == 1;

			my $indelLength = 0;
			if ($isIndel == 1) {
				my $refLength = length($ref);
				my $altLength = length($alt);

				$indelLength = abs($refLength - $altLength);
			}
			$skipINDEL++ if $indelLength >= $maxIndelLength;
			next if $indelLength >= $maxIndelLength;

			# calculate depth
			$_ =~ /DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;/;
			my $refDepth = $1 + $2; # this should actually always be 0 for the parents
			my $altDepth = $3 + $4;
			my $totalDepth = $refDepth + $altDepth;
			$skipDepth++ if $totalDepth <= $minParentalDepth;
			next if $totalDepth <= $minParentalDepth;

			#print "$chr\t$id\t$vcfScore\t$refDepth\t$altDepth\t$totalDepth\t$ref\t$alt\n";

			$parentalSNPs{$parent}{$chr}{$id} = "$ref|$alt|$vcfScore|$totalDepth";
			$variantPositions{$chr}{$id} = 1;
			++$SNPcount unless $_ =~ /INDEL/;
			++$indelCount if $_ =~ /INDEL/;
		}
		close INF;

		my $total = $SNPcount + $indelCount;

		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t repeat", $skipRepeat; 
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t score < $minVCFScore", $skipScore;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t het SNP", $skipHet;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t two alt alleles", $skipAlt;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t INDEL >= $maxIndelLength", $skipINDEL;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t depth <= $minParentalDepth", $skipDepth;
		printf "[".localtime(time)."] %30s %8d\n", "Total INDELs", $indelCount;
		printf "[".localtime(time)."] %30s %8d\n", "Total SNPs", $SNPcount;
		printf "[".localtime(time)."] %30s %8d\n", "Total variants", $total;
		print "[".localtime(time)."] \n";
	}
	print "[".localtime(time)."] Done gathering parental data.\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
	print "[".localtime(time)."] \n";

	## Grab parental consensus sequence data
	my %cns;
	print "[".localtime(time)."] Opening consensus file for each parental line - can be slow.\n";
	foreach my $parent ($parentA, $parentB) {
		print "[".localtime(time)."] - Opening consensus files for $parent.\n";
		foreach my $chr (@chrList) {
			my $cns = "$parentLoc/$parent.cns.$chr.gz";
        		open INF,"gunzip -c $cns |" or die "Can't open file $cns: $!\n";
			while (<INF>) {
				chomp($_);
				my(@F) = split /\t/, $_;
			
				$F[1] =~ tr/[a-z]/[A-Z]/;
                		$cns{$parent}{$chr}{$F[0]} = "$F[1]|$F[2]";
                	}
			close INF;
        	}
	}
	print "[".localtime(time)."] Done collecting consensus sequence. Memory usage ".memUsage()."\n";
	print "[".localtime(time)."] \n";

	## Identify positions that differ between $parentA and $parentB
	my %skip;
	my %diffVariants;
	my %diffVariantCounts;
	my $output = "Chr\tID\tRef\t$parentA\_cns\tScore\t$parentA\_vcf\tScore\tDepth\t$parentB\_cns\tScore\t$parentB\_vcf\tScore\tDepth\n";
	foreach my $chr (keys %variantPositions) {
		foreach my $id (sort {$a<=>$b} keys %{$variantPositions{$chr}}) {
			my($parentAcns,$parentAcnsScore) = $cns{$parentA}{$chr}{$id} =~ /(\w+)\|(\d+)/;
			my($parentBcns,$parentBcnsScore) = $cns{$parentB}{$chr}{$id} =~ /(\w+)\|(\d+)/;

			$skip{'missingCNS'}++ if !$parentAcns || !$parentBcns;
			next if !$parentAcns || !$parentBcns;

			$skip{'LowCNSScore'}++ if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;
			next if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;

			my($refA,$parentAalt,$AvcfScore,$AvcfDepth) = $parentalSNPs{$parentA}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;
			my($refB,$parentBalt,$BvcfScore,$BvcfDepth) = $parentalSNPs{$parentB}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;

			my $ref = $refA;
		   	$ref = $refB if !$ref;

			$parentAalt 	= "."  if !$parentAalt; 
			$AvcfScore 	= "..." if $parentAalt eq '.'; 
			$AvcfDepth	= ".." if $parentAalt eq '.'; 
			$parentBalt 	= "."  if !$parentBalt; 
			$BvcfScore 	= "..." if $parentBalt eq '.';
			$BvcfDepth 	= ".." if $parentBalt eq '.';

			$skip{'IdenticalAlt'}++ if $parentAalt eq $parentBalt;
			next if $parentAalt eq $parentBalt;

			$diffVariants{$chr}{$id} = "$parentAcns\t$parentAcnsScore\t$parentAalt\t$AvcfScore\t$AvcfDepth\t$parentBcns\t$parentBcnsScore\t$parentBalt\t$BvcfScore\t$BvcfDepth";
			$output .= "$chr\t$id\t$ref\t$diffVariants{$chr}{$id}\n";
			$diffVariantCounts{$chr}++;
	
			# need to add depth of coverage filter at the cns positions
		}
	}

	# deleting structures not needed. See this thread: http://www.perlmonks.org/?node_id=182343
	undef %parentalSNPs;
	undef %cns;

	open OUTF,">$pwd/out_uniqueParentalVariants.chr2.tsv";
	print OUTF $output;
	close OUTF;

	undef $output;

	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t missing base", $skip{'missingCNS'}; 
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t low CNS score", $skip{'LowCNSScore'}; 
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t identical variants", $skip{'IdenticalAlt'};

	printf "[".localtime(time)."] \n";
	printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
	my $totalCount;
	foreach my $chr (keys %diffVariantCounts) {
		#my $tmpChr = "chr".$chr;
		my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $diffVariantCounts{$chr});
		printf "[".localtime(time)."] %15s %8d - 1 snp / %4d bp\n", $chr, $diffVariantCounts{$chr}, $snpPerBP; 
		$totalCount += $diffVariantCounts{$chr};
	}
	printf "[".localtime(time)."] %30s %8d\n", "Total Count", $totalCount;
	print "[".localtime(time)."] \n";
	print "[".localtime(time)."] Unique variants saved in out_uniqueParentalVariants.tsv\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
	print "[".localtime(time)."] \n";
} else {
	die "Error: out_uniqueParentalVariants.csv doesn't exist! Quitting.\n" if !-e "$pwd/out_uniqueParentalVariants.tsv";
	print "[".localtime(time)."] -e flag passed, out_uniqueParentalVariants.tsv exists.\n";
}

exit 0 if $opts{'o'};

# open out_uniqueParentalVariants.csv
my(%parental,%SNP,$countVariants,%countSNPs);
print "[".localtime(time)."] Opening out_uniqueParentalVariants.tsv.\n";
open INF,"$pwd/out_uniqueParentalVariants.tsv" or die "Can't open out_uniqueParentalVariants.tsv: $!";
while (<INF>) {
	my(@F) = split /\t/, $_;
	next unless $F[1] =~ /[0-9]/;

	my($chr,$id,$ref,$wcns,$wvcf,$cscns,$csvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

	next if $wvcf ne $wcns && $wvcf !~ /\./;
	next if $csvcf ne $cscns && $csvcf !~ /\./;
	next if $wcns !~ /[A|G|C|T]/ || $cscns !~ /[A|G|C|T]/;
	next if $wcns eq $cscns;

	$SNP{$chr}{$id} = $ref;
	$parental{$parentA}{$chr}{$id} = $wcns;
	$parental{$parentB}{$chr}{$id} = $cscns;

	#print "$parentA\t$chr\t$id\t$wcns\t\t$parentB\t$cscns\n";

	++$countVariants;
	$countSNPs{$chr}++;
}
print "[".localtime(time)."] $countVariants total variants identified.\n";
printf "[".localtime(time)."] \n";
printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
foreach my $chr (@chrList) {
	my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $countSNPs{$chr});
	printf "[".localtime(time)."] %17s %10d - 1 snp / %4d bp\n", $chr, $countSNPs{$chr}, $snpPerBP; 
}
print "[".localtime(time)."] \n";
print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
print "[".localtime(time)."] \n";


######
if (-e "$pwd/out_cogcPositionsToSkip.tsv") {
	open INF,"$pwd/out_cogcPositionsToSkip.tsv";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		my($chr,$count) = ($F[0],$F[2]);
		my($id1,$id2) = $F[1] =~ /(\d+)\-(\d+)/;

		next if $count <= 2;
	
		foreach ($id1..$id2) {
			#print "Skipping: $chr\t$_\n";
			$parental{$parentA}{$chr}{$_} = undef;
			$parental{$parentB}{$chr}{$_} = undef;
		}
	}
	close INF;
}


my($sep,$stockCount,%eventCount,%eventStocks,$snpsAnalyzedOutput);

$sep= "+-----+------++-----------+-----+----++-----------+-----+----+-----------+-----++-----------+-----+----++----------+----------++-------------+-------------+";
#print "$sep\n";
#print "| Stk |  Chr ||   Last ID | Bas | Pa ||        ID | Bas | Pa |        ID | Bas ||   Next ID | Bas | Pa ||     fGap |     bGap || Lst->CurSNP | SNPsInEvent | \n";
#print "Chr\tStock\tLastParent\tRange\t\n";


my @files = `ls -1 $pwd`;
foreach my $stock (@files) {
	#next unless $stock =~ /^chr2/i;
	next unless $stock =~ /^chrX/i;

	chomp($stock);
	#print "[".localtime(time)."] $stock\n";
	print "\n";
	#print "$stock\n";
	$stockCount++;

	my(%vcfData,%allVCFData,%totalSNPsAnalyzed);

	# open VCF file, store data in %vcf
	my $vcfFile = "$stock.vcf.gz";
	next if !-e "$stock/$vcfFile";
	open INF,"gunzip -c $pwd/$stock/$vcfFile |" or die "Can't open $pwd/$stock/$vcfFile: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);
		$allVCFData{$chr}{$id} = 1;

		next unless $chrListToCheck{$chr} == 1;
		next unless $parental{$parentA}{$chr}{$id} && $parental{$parentB}{$chr}{$id};
		next unless $vcfScore >= $minVCFScore;	
               	next if $alt =~ /\,/; 		

		my $isIndel = 0;
	   	$isIndel = 1 if $_ =~ /INDEL/;
		next if $isIndel == 1 && $skipIndels == 1;

		my $type = "het";
		   $type = "hom" if $F[9] =~ /1\/1/;

               	#next if $type eq "het" && $chr eq "chrX"; 	# skip X chr calls that are het

		my $parent = "neither";
		
		# these are over ISO-1
		if ($type eq "hom" ) {
			$parent = $parentA if $alt eq $parental{$parentA}{$chr}{$id};
			$parent = $parentB if $alt eq $parental{$parentB}{$chr}{$id};
		} elsif ($type eq "het") {
			#$parent = "both" if $alt eq $parental{$parentA}{$chr}{$id} && $ref eq $parental{$parentB}{$chr}{$id};
			#$parent = "both" if $ref eq $parental{$parentA}{$chr}{$id} && $alt eq $parental{$parentB}{$chr}{$id};
			$parent = $parentA if $alt eq $parental{$parentA}{$chr}{$id};
			$parent = $parentB if $alt eq $parental{$parentB}{$chr}{$id};
		}

		$vcfData{$chr}{$id} = "$parent|$type|$ref|$alt";
		$totalSNPsAnalyzed{$chr}++;
		#print "$vcfData{$chr}{$id}\n" if $id < 1000000 && $chr eq "chr2L";
	}

	foreach my $chr (@chrList) {
		foreach my $id (keys %{$SNP{$chr}}) {
			next if $vcfData{$chr}{$id};
			next if $allVCFData{$chr}{$id} == 1;

			my $parent = "neither";
			$parent = $parentA if $SNP{$chr}{$id} eq $parental{$parentA}{$chr}{$id};
			$parent = $parentB if $SNP{$chr}{$id} eq $parental{$parentB}{$chr}{$id};

			$vcfData{$chr}{$id} = "$parent|ref|$SNP{$chr}{$id}|N";
			$totalSNPsAnalyzed{$chr}++;
			#print "$parent|ref|$SNP{$chr}{$id}|N\n" if $id == 15898686;
		}


		#print "[".localtime(time)."] $stock | $chr \n";
		my($lastParent,$lastParentChange,$lastID,$snpCount,$changeList);

		foreach my $id (sort {$a<=>$b} keys %{$vcfData{$chr}}) {
			my($parent,$type,$ref,$alt) = split /\|/, $vcfData{$chr}{$id};
			#print "$vcfData{$chr}{$id}\n" if $id == 15898686;
			#print "$vcfData{$chr}{$id}\n" if $id == 15898498;

			next if $parent eq "neither";

			if (!$lastParent) {
				$lastParent = $parent;
				$lastParentChange = $id;
				$snpCount = 1;
				$lastID = $id;

				print "\n";
				print "$stock\tBEGIN_CHR\t$parent\t$chr:$id\t$id\t$snpCount\n";

			} elsif ($parent eq $lastParent) {
				$snpCount++;
				$lastID = $id;
			} elsif ($parent ne $lastParent) {
				my($depth,$allele,$alleles) = getAllele($stock,$chr,$id);

				my $tmpParent = "neither";
				   $tmpParent = $parentA if $allele eq $parental{$parentA}{$chr}{$id} && $alleles !~ /\|\w\|/;
				   $tmpParent = $parentB if $allele eq $parental{$parentB}{$chr}{$id} && $alleles !~ /\|\w\|/;
				#print "$allele\t$alleles\t$depth\t$parent\t$tmpParent\n";

				if ($alleles =~ /(\w)\|(\w)\|/) {
					my($a,$b) = ($1,$2);
					$tmpParent = "both" if $a eq $parental{$parentA}{$chr}{$id} && $b eq $parental{$parentB}{$chr}{$id};
					$tmpParent = "both" if $b eq $parental{$parentA}{$chr}{$id} && $a eq $parental{$parentB}{$chr}{$id};
				}

				#print "$stock\t$chr\t$id\t$depth\t$allele\t$tmpParent\t$parent\n" if $id == 24865481;
				if ($depth > $minChildDepth && $tmpParent eq $parent) {
					#if ($changeList) {
						my $gap = $id - $lastID;
						my $lastGap = $lastID - $lastParentChange;
						#print "$changeList-$lastID\t$lastGap\t$snpCount\t->\t$parent\t$id\t$gap\n";
						print "$stock\tCHANGE\t$lastParent\t$chr:$lastID\t$lastID\t$snpCount\t$tmpParent\t$chr:$id\t$id\t$chr:$lastID-$id\t$gap\n";
						#$changeList = undef;

						if ($snpCount < 100) {
							$eventCount{$chr}{"$lastParentChange-$lastID"}++;
							$eventStocks{$chr}{"$lastParentChange-$lastID"} .= "$stock,";
						}
						
					#} else {
						#my $gap = $id - $lastID;
						#my $lastGap = $lastID - $lastParentChange;
						#$changeList = "$chr\t$stock\t$lastParent\t$lastParentChange-$lastID\t$snpCount\t$lastGap\t->\t$gap\t$parent\t$chr:$id";
					#}

					$lastID = $id;
					$lastParent = $parent;
					$lastParentChange = $id;
					$snpCount = 1;
				}
			}
		}
		#if ($changeList =~ /$chr/) {
			#my $lastGap = $lastID - $lastParentChange;
			#print "$changeList-$lastID\t$lastGap\t$snpCount\t-> END_OF_CHR\n";
		#}
		print "$stock\tEND_CHR\t$lastParent\t$chr:$lastID\t$lastID\t$snpCount\n";
	}

	foreach my $chr (sort keys %totalSNPsAnalyzed) {
		$snpsAnalyzedOutput .= "$stock\tSNPs\t$chr\t$totalSNPsAnalyzed{$chr}\n";
	}
}

if ($opts{'l'}) {
	open OUTF,">$pwd/out_cogcPositionsToSkip.tsv";
	print OUTF "Chr\tRange\tCount\tStocks Seen In (stocks checked: $stockCount)\n";
	foreach my $chr (keys %eventCount) {
		foreach my $id (sort {$a<=>$b} keys %{$eventCount{$chr}}) {
			my $count = $eventCount{$chr}{$id};
			#next unless $count > 3;
			my($id1,$id2) = $id =~ /(\d+)\-(\d+)/;
			my $gap = $id2 - $id1;
			print OUTF "$chr\t$id\t$count\t$gap\t$eventStocks{$chr}{$id}\n";
		}
	}
	close OUTF;
}

#####################

sub getAllele {
	my($stock,$chr,$id) = ($_[0],$_[1],$_[2]);

	my %alleleCount;

	foreach my $chunk (split /\n/, `samtools view $pwd/$stock/$stock.bam $chr:$id-$id`) {
		my(@F) = split /\t/, $chunk;
		next if $F[4] < 60;
		my($splitInfo,$originalSequence,$originalQuality) = ($F[5],$F[9],$F[10]);
	
       		my @seq;
       		foreach (split //, $originalSequence) {
       			push(@seq,$_);
       		}
       		$originalSequence = undef;
       		my $count = 0;
       		foreach (split //, $originalQuality) {
       			my $score = $qualityScores{$_};
       			my $base = $seq[$count];
       			$base = "n" if $qualityScores{$_} < 15;
       			$originalSequence .= $base;
       			++$count;
       		}

        	my $unknownAction = 0;
		my $seq;
        	foreach (split /([0-9]+[A-Z])/, $splitInfo) {
        		my($len,$action) = $_ =~ /([0-9]+)([A-Z])/;
        		next unless $len =~ /[0-9]/;

        		if ($action =~ /M/) {
        			$originalSequence =~ /^(.{$len})/;
        			$seq .= $1;
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /D/) {
        			foreach (1..$len) {
        				$seq .= ".";
        			}
        		} elsif ($action =~ /I/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /S/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} else {
				$unknownAction = 1;
			}
		}

		next if !$seq;
		my $length = length($seq);
		my $startID = $F[3];
		my $endID = $F[3] + $length - 1;

		$startID = $endID if $startID > $id;

		my $currID = $startID;
		foreach (split //, $seq) {
			#print "$currID\t$_\n" if $currID == 7023293;
			$alleleCount{$_}++ if $currID == $id;
			$currID++;
		}

	}

	my $readCount = $alleleCount{'A'} + $alleleCount{'G'} + $alleleCount{'C'} + $alleleCount{'T'};
	my($allele,$alleles);
	if ($readCount > 0) {
		my $tmpReadCount = $readCount / 2;
		$allele = "A" if $alleleCount{'A'} >= $tmpReadCount;
		$allele = "G" if $alleleCount{'G'} >= $tmpReadCount;
		$allele = "C" if $alleleCount{'C'} >= $tmpReadCount;
		$allele = "T" if $alleleCount{'T'} >= $tmpReadCount;

		my $alleleFreq = sprintf("%0.0f", (($alleleCount{$allele} / $readCount ) * 100));

		if ($alleleFreq >= 30 && $alleleFreq <= 70) {
			$alleles .= "A|" if $alleleCount{'A'} >= 2;
			$alleles .= "G|" if $alleleCount{'G'} >= 2;
			$alleles .= "C|" if $alleleCount{'C'} >= 2;
			$alleles .= "T|" if $alleleCount{'T'} >= 2;
		} elsif ($alleleFreq >= 5 && $alleleFreq <= 95) {
			$allele = "N";
		}

		#print "$alleleFreq\t$allele\t$alleles\t$alleleCount{'A'}, $alleleCount{'G'}, $alleleCount{'C'}, $alleleCount{'T'}\n" if $id == 24865481;
	}

	return($readCount,$allele,$alleles);
}

open OUTF,">out_SNPsAnalyzed.tsv";
print OUTF $snpsAnalyzedOutput;
close OUTF;
