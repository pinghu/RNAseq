#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  <badfile> <file> \n"; 
die $usage unless @ARGV == 2;
my ($badfile,  $file) = @ARGV;
my %badid;
open (A, "<$badfile")||die "could not open $badfile\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, lc($a);	
    $badid{$tmp[0]}=0;
}
close A;

open (A, "<$file")||die "could not open $file\n";
my $tt=<A>;
print $tt;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, lc($a);
    if (defined $badid{$tmp[0]}){
	    print STDERR "remove ", $tmp[0], "\n";
    }else{
	print $a, "\n";	    	
    }
}
close A;

