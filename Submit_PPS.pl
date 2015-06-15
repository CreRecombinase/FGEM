#!/usr/bin/perl
use strict;
use warnings;

my $homedir = $ARGV[0];
my $shelldir = $ARGV[1];
my $GOfile = $ARGV[2];
my $BFfile= $ARGV[3];
my $startfeat = $ARGV[4];
my $lastfeat = $ARGV[5];
my $featsize =$ARGV[6];

my @array = ($startfeat .. $lastfeat);

for{ my $item =$startfeat; $item<$lastfeat; $item+=$featsize){
    my $tstop = $item+$featsize;
    my $shellfile = "${shelldir}/FGEM_${item}_${tstop}";
    print $shellfile."\n";
}
	
