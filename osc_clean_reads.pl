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

my $dir="/fs/scratch/ar2767/OralBGIRNASeq092019/orig/";
my $rstD="/fs/scratch/ar2767/OralBGIRNASeq092019/clean/";


open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    
    $count++;
    open (B, ">OSC.HN.$count")||die "could not open OSC.HN.$count\n";
    print B "#PBS -l walltime=2:00:00\n";
    print B "#PBS -l nodes=1:ppn=28\n";
    print B "#PBS -S /bin/bash\n";
    #print B "module load blast\n";
    #print B "module load bowtie2\n";
    #print B "module load samtools\n";
    #print B "module load bwa\n";
    print B "source /users/PYS0226/ar2767/.bash_profile\n";
    print B "gunzip $dir/$a*.gz\n";
    
    #print B "/users/PYS0226/ar2767/pkg/hisat2-2.1.0/hisat2 -p 28 --dta -x /users/PYS0226/ar2767/project/RnR/71_3/TE31_71_3-2011.tran -1 $dir/$a"."_1.fq -2 $dir/$a"."_2.fq -S $rstD/$a.sam\n";
    #print B "samtools sort -@ 28 -o $rstD/$a.bam $rstD/$a.sam\n";
  print B "/users/PYS0226/ar2767/.local/bin/cutadapt -q 30 --minimum-length 60 -o $rstD/$a"."_1.fq -p $rstD/$a"."_2.fq $dir/$a"."_1.fq $dir/$a"."_2.fq 2>$rstD/$a.err\n";  
#print B "cd $rstD\n";
 #   print B "/users/PYS0226/ar2767/pkg/stringtie-1.3.3b.Linux_x86_64/stringtie -p 28 -G /users/PYS0226/ar2767/project/RnR/71_3/TE31_71_3-2011.gtf -o $rstD/$a.gtf -l $rstD/$a  $rstD/$a.bam\n";
    close B;
    print "qsub OSC.HN.$count\n";
}
close A;

