#!/usr/bin/perl

use strict;
use warnings;

my $sample = shift or die $!;
my $dir = shift or die $!;

my %hash;
my @sam;
open IN, "$sample" or die $!;
my $header = <IN>;
die "Empty samplesheet: $sample\n" unless defined $header;
my $sep = ($header =~ /\t/ && $header !~ /,/) ? "\t" : ",";
while(<IN>){
	chomp;
	next unless /\S/;
	my @tmp = split /\Q$sep\E/;
	my $id = $tmp[0];
	$id =~ s/^\s+|\s+$//g;
	push @sam, $id;
	my $count_file = "$dir/$id.count";
	open INI, $count_file or die "Cannot open count file '$count_file': $!";
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
