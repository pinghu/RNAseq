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
open (A, "</mnt/G6_2D/db/ko/ko_name.list")||die "could not open ko.name\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
   
    $data{$tmp[0]}=$a;
    
   
}
close A;

my %syn;
open (A, "</mnt/G6_2D/db/ko/ko.symbol")||die "could not open ko.symbol\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
   
    $syn{$tmp[0]}=$tmp[1];
    
   
}
close A;
open (A, "<$file")||die "could not open $file\n";
#my $tt=<A>;
#chomp $tt;
#for($tt){s/\r//gi;}
#print "$tt\tKOSyn\tKOAnn\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    my $id=$tmp[$N-1];
    if($id =~ /ko\:(\S+)/){$id=$1;}
    my $ann="NA";
    if(defined $data{$id}){
	$ann=$data{$id};
    }else{
	print STDERR $tmp[$N-1], " ", $id,  " no ann \n";
    }
    my $sss="NA";
    if(defined $syn{$id}){
	$sss=$syn{$id};
    }else{
	print STDERR $tmp[$N-1], , " ", $id," no symbol\n";
    }
    print $a, "\t", $sss, "\t", $ann, "\n";
}
close A;

