#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;

open (A, "<$file")||die "could not open $file\n";
open (B, ">kegg.job")||die "could not open kegg.job\n";
open (C, ">stat2KO.job")||die "could not open stat2KO.job\n";
open (D, ">enrichment.job")||die "could not open enrichment.job\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
 

    print B "perl ../src/test_kegg_compound_new.pl ../ann/Ecoli_BW25113.koid ../stat/$a ../ann/Ecoli_BW25113.name ../ann/Ecoli_BW25113.symbol ../ann/Ecoli_BW25113.geneid -l ../ann/Ecoli_BW25113.loc 1>$a.out 2>$a.err\n";
    
    print C "perl ../src/stat_2_KO.pl  ../stat/$a\n";
    print D "Rscript ../src/enrichment_test.R ../stat/$a.KO\n";

}
close A;
close B;
close C;
close D;
system("chmod 777 *job\n");
