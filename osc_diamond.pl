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
die $usage unless @ARGV >=1;
my ($file) = @ARGV;

my $dir="/fs/scratch/ar2767/OralBGIRNASeq092019/Ann/protein/";
my $rstD="/fs/scratch/ar2767/OralBGIRNASeq092019/Ann/protein/";


open (A, "<$file")||die "could not open $file\n";
my $count=0;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    
    $count++;
    open (B, ">OSC.HN.$count")||die "could not open OSC.HN.$count\n";
    print B "#PBS -l walltime=120:00:00\n";
    print B "#PBS -l nodes=1:ppn=28\n";
    print B "#PBS -S /bin/bash\n";
    print B "source /users/PYS0226/ar2767/.bash_profile\n";
    
   
    print B "/users/PYS0226/ar2767/bin/diamond blastp -d /fs/project/PYS0226/ar2767/diamond_database/KEGG_prokaryotes_pep -q $dir/$a --threads 28 -o $rstD/$a.keggPro.tbl 2>$rstD/$a.keggPro.err\n";
    print B "/users/PYS0226/ar2767/bin/diamond blastp -d /fs/project/PYS0226/ar2767/diamond_database/VFDB_all_proteins.dmnd -q $dir/$a --threads 28 -o $rstD/$a.vfdb.tbl 2>$rstD/$a.vfdb.err\n";
     print B "/users/PYS0226/ar2767/bin/diamond blastp -d /fs/project/PYS0226/ar2767/diamond_database/T3DB_protein.dmnd -q $dir/$a --threads 28 -o $rstD/$a.T3DB.tbl 2>$rstD/$a.T3DB.err\n";
    print B "/users/PYS0226/ar2767/bin/diamond blastp -d /fs/project/PYS0226/ar2767/diamond_database/uniref100 -q $dir/$a --threads 28 -o $rstD/$a.uniref100.tbl 2>$rstD/$a.uniref100.err\n";
    print B "/users/PYS0226/ar2767/bin/diamond blastp -d /fs/project/PYS0226/ar2767/diamond_database/nr -q $dir/$a --threads 28 -o $rstD/$a.nr.tbl 2>$rstD/$a.nr.err\n";
   
    close B;
    print "qsub OSC.HN.$count\n";
}
close A;

