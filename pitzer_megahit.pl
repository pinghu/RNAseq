#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### https://www.osc.edu/resources/available_software/software_list/python
#https://www.osc.edu/resources/getting_started/howto/howto_locally_installing_software
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;
my $rstM="/fs/scratch/PYS0226/pg123/GSS2910_11Thymol/merge/";
my $rstS="/fs/scratch/PYS0226/pg123/GSS2910_11Thymol/megahit/";


open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    open (B, ">OSC.megahit.$a")||die "could not open OSC.megahit.count\n";
    print B "#!/bin/bash\n";
    print B "#SBATCH --time=12:00:00\n";
    print B "#SBATCH --nodes=1 --ntasks=28\n";
    print B "#SBATCH --account=PYS0226\n";
    print B "\n";
    print B "source /users/PYS0226/pg123/.bashrc\n";
   
    print B "module load python\n";
    print B "module load bwa\n"; 
    #print B "gunzip $rstM/$a.R1.fastq.gz\n";
    #print B "gunzip $rstM/$a.R2.fastq.gz\n";
    print B "time /users/PYS0226/pg123/pkg/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit  -1 $rstM/$a.R1.fastq -2 $rstM/$a.R2.fastq -o $rstS/$a \n";
    close B;
    print "sbatch OSC.megahit.$a\n";
}
close A;
