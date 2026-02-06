#!/usr/bin/perl

use strict;
use warnings;

my $sample = shift or die $!;
my $dir = shift or die $!;

my %hash;
my @sam;
open IN, "$sample" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $id = $tmp[0];
	push @sam, $id;
	open INI, "$dir/$id.count " or die $!;
	while(<INI>){
		chomp;
#		next if (/^__/);
		next if (/^#/);
		next if (/^Geneid/);
		my @tmp = split /\t/;
		my $esgn = $tmp[0];
		my $value = $tmp[-1];
		$hash{$esgn}{$id} = $value;
	}
	close INI;
}
close IN;

my $header = join("\t", @sam);
print "GID\t$header\n";

for my $i (sort {$a cmp $b} keys %hash){
	print "$i";
	for my $j (@sam){
		print "\t$hash{$i}{$j}";
	}
	print "\n";
}

