#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### https://www.osc.edu/resources/available_software/software_list/python
#https://www.osc.edu/resources/getting_started/howto/howto_locally_installing_software
##software located at ~/.local/bin
#pip install humann2
#humann2_databases --download chocophlan full humann2_database_downloads
#humann2_databases --download uniref uniref90_ec_filtered_diamond humann2_database_downloads
#humann2_databases --download utility_mapping full humann2_database_downloads
#humann2_databases --download utility_mapping full humann2_database_downloads
#pip install biom-format --user
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;

my $orig="/fs/scratch/PYS0226/pg123/GSS2910_11Thymol/org/";
my $dir="/fs/scratch/PYS0226/pg123/GSS2910_11Thymol/merge/";


open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    
    $count++;
    open (B, ">OSC.HN.$count")||die "could not open OSC.HN.$count\n";
    print B "#!/bin/bash\n";
    print B "#SBATCH --time=12:00:00\n";
    print B "#SBATCH --nodes=1 --ntasks-per-node=1\n";
    print B "#SBATCH --account=PYS0226\n";
    print B "\n";
    print B "source /users/PYS0226/pg123/.bashrc\n";
   
    print B "cat $orig/$a","_L001_R1_001.fastq.gz $orig/$a","_L002_R1_001.fastq.gz  $orig/$a","_L003_R1_001.fastq.gz $orig/$a","_L004_R1_001.fastq.gz  >$dir/$a.R1.fastq.gz\n";
    print B "cat $orig/$a","_L001_R2_001.fastq.gz $orig/$a","_L002_R2_001.fastq.gz  $orig/$a", "_L003_R2_001.fastq.gz $orig/$a","_L004_R2_001.fastq.gz  >$dir/$a.R2.fastq.gz\n";  
    close B;
    print "qsub OSC.HN.$count\n";
}
close A;
