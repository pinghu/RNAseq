#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
###ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
###ncbi gene data
####ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz 
#### gene list file will be the [uniref] [uniref90]
####/mnt/G6_2D/project/oralKF/refseq/idmapping.GO
####eggnog: http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz
##########################################################################
use strict;
my $usage="$0  [listfile][col][category ann file] first line is title\n"; 
die $usage unless @ARGV ==3;
my ($file, $N, $file2) = @ARGV;
my %data;
open (A, "<$file2")||die "could not open \n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
   
    $data{$tmp[0]}=$tmp[2];
    
   
}
close A;
open (A, "<$file")||die "could not open $file\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    my $id=$tmp[$N-1];
    if($id=~/UniRef\d+\_(\S+)/){$id=$1; print STDERR $1, "\n";}
    if(defined $data{$id}){
	print $a, "\t", $data{$id}, "\n";
    }else{
	print $a, "\tNA\n";
	#print STDERR $tmp[$N-1], "\n";
    }
   
}
close A;

