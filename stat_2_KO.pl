#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### gene list file will be the [uniref] [uniref90]
##########################################################################
use strict;
my $usage="$0  [listfile][category ann file] first line is title\n"; 
die $usage unless @ARGV >=1;
my ($file, $KOMAP) = @ARGV;
if(! defined $KOMAP){
    $KOMAP="/Disk2/project/OralTrans/starRsen/ID2KO";
}
my %data;
open (A, "<$KOMAP")||die "could not open \n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    if(!defined $tmp[1]){
	#print STDERR $a, "\n";
	next;
    }
    if($tmp[1] =~/ko\:(K\d+)/){
	$tmp[1]=$1;
    }
    $data{$tmp[0]}=$tmp[1];
}
close A;
open (A, "<$file")||die "could not open $file\n";
open (B, ">$file.KO")||die "could not open $file.KO\n";
my $tt=<A>;
print B "KOID\tFC\tPvalue\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    my $id=$tmp[0];
    if(defined $data{$id}){
	print  B $data{$id},"\t", $tmp[2], "\t", $tmp[1], "\n";
    }else{
	print STDERR $tmp[0], " no mapping\n";
    }
   
}
close A;
close B;

