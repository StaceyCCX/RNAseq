#!/usr/bin/perl -w
use strict;


##this script is used to buqi the output of topgo
#
die '$ARGV[0] is the /psc/home/zhaocheng/species/GO.terms_and_ids, $ARGV[1] is the topgo output file' if @ARGV !=2;

open A,$ARGV[0] or die "$ARGV[0] can't be opened";
open B,$ARGV[1] or die "$ARGV[1] can't be opened";

my %GO;
while (<A>) {
	if ($_=~/^!/) {next};
	chomp;
	my @array=split "\t",$_;
	$GO{$array[0]}=$array[1];
}
close A;
while (<B>) {
	chomp;
	my @array=split "\t",$_;
	if (not exists $GO{$array[0]}) {$GO{$array[0]}=$array[1]};
	$array[1]=$GO{$array[0]};
	my $str=join "\t",@array;
	print "$str\n";
}
close B;


