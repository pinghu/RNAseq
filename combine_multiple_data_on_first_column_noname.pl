#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;

my (@files) = @ARGV;

my %data;
my %title;
my %genes;

for (my $i =0; $i<@files; $i++){
    open (C, "<$files[$i]")|| die "could not open $files[$i]\n";
    my $t=<C>;
    chomp $t;for($t){s/\r//g;}
    my @ss=split /\t/, $t;
    for(my $j=0; $j<@ss; $j++){
	$ss[$j]=$ss[$j];
    }
    $t=join("\t", @ss);
    $title{$files[$i]}=$t;
    while (my $M =<C>){
	chomp $M;
	for($M){s/\r//g;}
	my @tmp=split /\t/, $M;	
	if ($tmp[0] =~/^\s*$/){next;}
	$data{lc($tmp[0])}{$files[$i]}=$M;
	$genes{lc($tmp[0])}=$tmp[0];
    }
    close C;
}

print "gene";
for (my $i =0; $i<@files; $i++){
    print "\t",$title{$files[$i]};
}
print "\n";
foreach my $j (keys %genes){
    print $j;
    for (my $i =0; $i<@files; $i++){
	print "\t";
	if (defined $data{$j}{$files[$i]}){
	    my $xx= $data{$j}{$files[$i]};
	    print $xx;
	    
	}else{
	    my @x=split /\t/,  $title{$files[$i]};
	    for (my $i=0; $i<@x; $i++){
		if ($i == @x-1){
		    print "NA";
		}else{
		    print "NA", "\t";
		}
	    }
	}
    }
    print "\n";
}
