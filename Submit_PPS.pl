#!/usr/bin/perl
use strict;
use warnings;


my $shelldir = $ARGV[0];
my $startfeat = $ARGV[1];
my $lastfeat = $ARGV[2];
my $featsize =$ARGV[3];




for( my $item = $startfeat; $item<$lastfeat; $item+=$featsize){
    
    my $tstop = $item+$featsize;
    if($tstop> $lastfeat){
	$tstop = $lastfeat;
    }
    my $shellbase = "${shelldir}/FGEM_${item}_${tstop}";
    my $shellout = "${shellbase}.out";
    my $shellerr = "${shellbase}.err";
    my $shellfile = "${shellbase}.sh";
    print $shellfile."\n";

    unless (-e $shellfile){
	open(SCRIPT,">${shellfile}");
	print SCRIPT "#!/bin/sh\n";
	print SCRIPT "#\$ -o ${shellout}\n";
	print SCRIPT "#\$ -e ${shellerr}\n";
	print SCRIPT "#\$ -l h_vmem=1g\n";
	print SCRIPT "/mnt/gluster/home/nwknoblauch/bin/Rscript /mnt/gluster/home/nwknoblauch/Gene_Annotation/FGEM/Bayes_EM_GWAS.R /mnt/gluster/home/nwknoblauch/Gene_Annotation/GOmat.RDS /mnt/gluster/home/nwknoblauch/Gene_Annotation/TADA_ASC_SSC_results_Dec23.csv ${item} ${featsize} /mnt/gluster/home/nwknoblauch/Gene_Annotation/ \n";
	close(SCRIPT);
	system("qsub ${shellfile}");
	print "${item}\n";
	sleep(1);

    
    }
}	
