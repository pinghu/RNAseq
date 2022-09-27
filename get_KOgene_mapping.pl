#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### gene list file will be the [uniref] [uniref90]
##########################################################################
use strict;
my $usage="$0  [listfile][col] first line is title\n"; 
die $usage unless @ARGV ==2;
my ($file, $N) = @ARGV;
my %data;
open (A, "</mnt/G6_2D/db/ko/ko_genes.list")||die "could not open ko_genes.list \n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
   
    $data{$tmp[1]}=$tmp[0];
    
   
}
close A;
open (A, "<$file")||die "could not open $file\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    my $id=$tmp[$N-1];
    if(defined $data{$id}){
	print $a, "\t", $data{$id}, "\n";
    }else{
	print $a, "\tNA\n";
	print STDERR $tmp[$N-1], "\n";
    }
   
}
close A;

