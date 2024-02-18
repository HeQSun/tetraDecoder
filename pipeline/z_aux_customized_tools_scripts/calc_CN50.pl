#! /usr/bin/perl
use POSIX;
use strict;
use warnings;

my $SMALL_CONTIG_SIZE = 10000;

my $usage = "$0 fasta [genomesize in bp] [chromosome number]\n";
my $FILE = shift or die $usage;
my $SIZE = shift or die $usage;
my $CHR_NUM = shift or die $usage;

###########################################################
# Read in file
open FILE, $FILE or die $usage;

my @CTG_LGTH = ();
my $NUM_N = 0;
my $CUM_LGTH = 0;
my $CTG_C = 0;
my $CTG_SMALL_C = 0;
my $CTG_SMALL_LGTH = 0;

my $seq = "";

while (my $line = <FILE>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		if ($seq ne "") {
			push @CTG_LGTH, length($seq);
			$CUM_LGTH += length($seq);
			for (my $i=0; $i<length($seq); $i++) { $NUM_N++ if (uc(substr($seq, $i, 1)) eq "N"); }
			if (length($seq) <= $SMALL_CONTIG_SIZE) {
				$CTG_SMALL_C++;
				$CTG_SMALL_LGTH+=length($seq);
			}
		}
		$seq = "";
		$CTG_C++;
	}
	else {
		$seq .= $line;
	}
}
if ($seq ne "") {
	push @CTG_LGTH, length($seq);
	$CUM_LGTH += length($seq);
	for (my $i=0; $i<length($seq); $i++) { $NUM_N++ if (uc(substr($seq, $i, 1)) eq "N"); }
	if (length($seq) <= $SMALL_CONTIG_SIZE) {
        	$CTG_SMALL_C++;
                $CTG_SMALL_LGTH+=length($seq);
        }
}
close FILE or die "cannot close file\n";
@CTG_LGTH = sort {$b <=> $a} @CTG_LGTH;
#print "A:", join(",", @CTG_LGTH), "\n";

###########################################################
# Analyze lengths and number

my ($N_ref, $L_ref, $NG_ref, $LG_ref, $MIN, $MAX) = calc_N50(\@CTG_LGTH);
my ($CN50, $CL50) = calc_CN50(\@CTG_LGTH);

print "EstGenomeSize:\t\t$SIZE\n";
print "ContigCount:\t\t$CTG_C\n";
print "ContigLength:\t\t$CUM_LGTH\n";
print "AmbiguousBases:\t\t$NUM_N\n";
print "N25:\t\t${$N_ref}[0]\n";
print "N50:\t\t${$N_ref}[1]\n";
print "N75:\t\t${$N_ref}[2]\n";
print "N90:\t\t${$N_ref}[3]\n";
print "L25:\t\t${$L_ref}[0]\n";
print "L50:\t\t${$L_ref}[1]\n";
print "L75:\t\t${$L_ref}[2]\n";
print "L90:\t\t${$L_ref}[3]\n";
print "NG25:\t\t${$NG_ref}[0]\n";
print "NG50:\t\t${$NG_ref}[1]\n";
print "NG75:\t\t${$NG_ref}[2]\n";
print "NG90:\t\t${$NG_ref}[3]\n";
print "LG25:\t\t${$LG_ref}[0]\n";
print "LG50:\t\t${$LG_ref}[1]\n";
print "LG75:\t\t${$LG_ref}[2]\n";
print "LG90:\t\t${$LG_ref}[3]\n";
print "Min:\t\t$MIN\n";
print "Max:\t\t$MAX\n";
print "TotalSmallContigs(<$SMALL_CONTIG_SIZE):\t$CTG_SMALL_C\n";
print "SmallContigsLength(<$SMALL_CONTIG_SIZE):\t$CTG_SMALL_LGTH\n";
print "CN50:\t\t$CN50\n";
print "CL50:\t\t$CL50\n";

sub calc_CN50 {
	my ($ctg_lgth_ref) = @_;

	my @N50 = ();
	my @L50 = ();
	my @chr_ctgs = ();
	my $count = 0;

	# set up contig-length array for each chromosome
	for (my $c = 0; $c < $CHR_NUM; $c++) {
		my @n = ();
		push @chr_ctgs, \@n;
	}

	# build up batches of contigs represent average subassemblies of chromosomes
	while ($count < $CTG_C) {
		for (my $c = 0; $c < $CHR_NUM and $count < $CTG_C; $c++) {
			push @{$chr_ctgs[$c]}, ${$ctg_lgth_ref}[$count];
			$count++;
		}

		for (my $c = $CHR_NUM-1; $c >= 0 and $count < $CTG_C; $c--) {
                        push @{$chr_ctgs[$c]}, ${$ctg_lgth_ref}[$count];
                        $count++;
                }
	}

	# ... and calculate N50 for each of them.
	for (my $c = 0; $c < $CHR_NUM; $c++) {
#print "$c:", join(",", @{$chr_ctgs[$c]}), "\n";
		my ($N_ref, $L_ref, $NG_ref, $LG_ref, $MIN, $MAX) = calc_N50($chr_ctgs[$c]);
		push @N50, ${$N_ref}[1];
		push @L50, ${$L_ref}[1];
	}

	# ... and select median ctg as CN50 contig.
	@N50 = sort {$b <=> $a} @N50;
	@L50 = sort {$b <=> $a} @L50;
	my $median_i = ceil($CHR_NUM / 2) - 1;	

	return ($N50[$median_i], $L50[$median_i]);

}

sub calc_N50 {
	my ($ctg_lgth_ref) = @_;

	my $ctg_count = $#{$ctg_lgth_ref}+1;
	my $cum_lgth = 0;
	my $max = 0;
        my $min = -1;

	for (my $i = 0; $i < $ctg_count; $i++) {
		my $length = ${$ctg_lgth_ref}[$i];
		$cum_lgth += $length;
		# set max length
                if ($max < $length) {
                	$max = $length;
                }
                if ($min == -1 || $min > $length) {
                	$min = $length;
                }
	}

	my @N = (); # N25: $N[0], N50: $N[1], N75: $N[2], N90: $N[3] 
	my @L = (); # L25: $L[0], L50: $L[1], L75: $L[2], L90: $L[3] 
	my @NG = (); # N25: $N[0], N50: $N[1], N75: $N[2], N90: $N[3] 
        my @LG = (); # L25: $L[0], L50: $L[1], L75: $L[2], L90: $L[3]

	my $add_lgth = 0;
	my $ctg_c = 0;

	for (my $i = 0; $i < $ctg_count; $i++) {
		my $ctg_lgth = ${$ctg_lgth_ref}[$i];
		$add_lgth += $ctg_lgth;
        	$ctg_c += 1;

        	# set N values
		if  (not defined($N[0]) and $add_lgth >= $cum_lgth*0.25) {
                	$L[0] = $ctg_c;
                        $N[0] = $ctg_lgth;
                }
		if  (not defined($N[1]) and $add_lgth >= $cum_lgth*0.5) {
                        $L[1] = $ctg_c;
                        $N[1] = $ctg_lgth;
                }
		if  (not defined($N[2]) and $add_lgth >= $cum_lgth*0.75) {
                        $L[2] = $ctg_c;
                        $N[2] = $ctg_lgth;
                }
		if  (not defined($N[3]) and $add_lgth >= $cum_lgth*0.9) {
                        $L[3] = $ctg_c;
                        $N[3] = $ctg_lgth;
                }
			
		# set NG values
		if  (not defined($NG[0]) and $add_lgth >= $SIZE*0.25) {
                        $LG[0] = $ctg_c;
                        $NG[0] = $ctg_lgth;
                }
                if  (not defined($NG[1]) and $add_lgth >= $SIZE*0.5) {
                        $LG[1] = $ctg_c;
                        $NG[1] = $ctg_lgth;
                }
                if  (not defined($NG[2]) and $add_lgth >= $SIZE*0.75) {
                        $LG[2] = $ctg_c;
                        $NG[2] = $ctg_lgth;
                }
                if  (not defined($NG[3]) and $add_lgth >= $SIZE*0.9) {
                        $LG[3] = $ctg_c;
                        $NG[3] = $ctg_lgth;
                }	
	}

	for (my $i = 0; $i < 4; $i++) {
		if  (not defined($N[$i])) {
			$L[$i] = 0;
                	$N[$i] = 0;
		}
		if  (not defined($NG[$i])) {
                	$LG[$i] = 0;
                	$NG[$i] = 0;
		}
        }	

	return (\@N, \@L, \@NG, \@LG, $min, $max);

}	


