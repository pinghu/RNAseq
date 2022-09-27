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
my $indir="/fs/scratch/ar2767/OralBGIRNASeq092019/star_step1/";
my $REFindex="/fs/scratch/ar2767/OralBGIRNASeq092019/Ann/RSEM.combine6DB/oral6";
my $outdir="/fs/scratch/ar2767/OralBGIRNASeq092019/rsem_step2/";

open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    
    $count++;
    open (B, ">OSC.rsem.$count")||die "could not open OSC.HN.$count\n";
    print B "#PBS -l walltime=24:00:00\n";
    print B "#PBS -l nodes=1:ppn=28\n";
    print B "#PBS -S /bin/bash\n";
    print B "module load blast\n";
    print B "module load bowtie2\n";
    print B "module load bowtie1\n";
    print B "module load samtools\n";
    print B "module load bwa/0.7.17-r1198\n";
    print B "module load star\n";
    print B "source /users/PYS0226/ar2767/.bash_profile\n";
    
    print B "/users/PYS0226/ar2767/pkg/RSEM-1.3.1/rsem-calculate-expression -p 12 --paired-end --bam --estimate-rspd --append-names --output-genome-bam $indir/$a"."Aligned.toTranscriptome.out.bam $REFindex $outdir/$a\n"; 
###STAR --runThreadN 28 --genomeDir $REFindex --sjdbGTFfile $REFgtf --sjdbOverhang $REFlen --readFilesIn $dir/$a.R1.fastq $dir/$a.R2.fastq  --outSAMtype BAM Unsorted  --outFileNamePrefix $rstD/$a --quantMode TranscriptomeSAM GeneCounts --chimOutType WithinBAM --chimSegmentMin 20 --outSAMstrandField introMotif --genomeLoad NoSharedMemory --twopassMode Basic\n";
    
    close B;
    print "qsub OSC.rsem.$count\n";
}
close A;
