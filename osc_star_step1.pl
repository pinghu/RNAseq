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

my $dir="/fs/scratch/ar2767/OralBGIRNASeq092019/clean/";
my $rstD="/fs/scratch/ar2767/OralBGIRNASeq092019/star_step1/";
my $REFindex="/fs/scratch/ar2767/OralBGIRNASeq092019/Ann/STAR.combine6DB/";
my $REFgtf="/fs/scratch/ar2767/OralBGIRNASeq092019/Ann/oral6.gtf";
my $REFlen=149;

open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    
    $count++;
    open (B, ">OSC.HN.$count")||die "could not open OSC.HN.$count\n";
    print B "#PBS -l walltime=24:00:00\n";
    print B "#PBS -l nodes=1:ppn=28\n";
    print B "#PBS -S /bin/bash\n";
    print B "module load blast\n";
    print B "module load bowtie2\n";
    print B "module load samtools\n";
    print B "module load bwa/0.7.17-r1198\n";
    #print B "module load start\n";
    print B "source /users/PYS0226/ar2767/.bash_profile\n";
    
    print B "/users/PYS0226/ar2767/pkg/STAR-2.6.0a/bin/Linux_x86_64/STAR --runThreadN 28 --genomeDir $REFindex --sjdbGTFfile $REFgtf --sjdbOverhang $REFlen --readFilesIn $dir/$a"."_1.fq $dir/$a"."_2.fq  --outSAMtype BAM Unsorted --outFileNamePrefix $rstD/$a --quantMode TranscriptomeSAM GeneCounts --chimOutType WithinBAM --chimSegmentMin 20 --outSAMstrandField introMotif --genomeLoad NoSharedMemory --twopassMode Basic\n";
    
    close B;
    print "qsub OSC.HN.$count\n";
}
close A;
