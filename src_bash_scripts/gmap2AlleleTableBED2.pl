#!/usr/bin/perl -w

die "Usage: perl $0 ref.bed\n" if(!defined ($ARGV[0]));
my $refGFF = $ARGV[0];
open(IN, "grep 'gene' gmap.gff3 |") or die"";
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    my $gene = $1 if(/Name=(\S+)/);
    $infordb{$gene} .= $data[0]."\t";
    }   
close IN; 

open(OUT, "> Allele.ctg.table") or die"";
open(IN, $refGFF) or die"";
while(<IN>){
    chomp;
    my @data = split(/\s+/,$_);
    my $gene = $data[3];
       $gene =~ s/;.*//g;
    next if(!exists($infordb{$gene}));
    my @tdb = split(/\s+/,$infordb{$gene});
    my %tmpdb = (); 
    map {$tmpdb{$_}++} @tdb;
    print OUT $data[0]."\t".$data[1]."\t"; # shq: get chr \t position \t a1 \t a2 \t a3 \t a4
    map {print OUT $_."\t"} keys %tmpdb;
    print OUT "\n";
    }   
close IN; 
close OUT;
