#!/usr/bin/perl

use strict;
use warnings;


my $count = shift or die $!;
my $genelen = shift or die $!;
my $prefix = shift or die $!;

#my $genelen = "/dssg/home/acct-medty/medty-c/database/GENCODE/gene_len.v43.txt";

my %ensg_len;
open IN, "$genelen" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $id = (split (/\./, $tmp[0]))[0];
	$ensg_len{$id}{'length'} = $tmp[3];
	$ensg_len{$id}{'symbol'} = $tmp[1];
	$ensg_len{$id}{'type'} = $tmp[2];

}
close IN;

my %count2tpm;
#my @all_gene;
my %all_gene;
#####################Step1: cal transcript = gene reads / gene length
open IN, "$count" or die $!;
chomp(my $sample = <IN>);
my @all_sample = split (/\t/, $sample);
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $id = (split (/\./, $tmp[0]))[0];
	next unless (defined $ensg_len{$id});
#	push @all_gene, $id;
	(defined $all_gene{$id}) ? next : ($all_gene{$id} = 1);
	for(my $i=1;$i<@tmp;$i++){
		$count2tpm{$all_sample[$i]}{$id} = $tmp[$i]/$ensg_len{$id}{'length'};
	}
}

my %sum_trans;
my %sum_count;
#####################Step2: sum all genes transcript per sample 
shift @all_sample;
for my $k (@all_sample){
	for my $j (keys %{$count2tpm{$k}}){
		$sum_trans{$k} += $count2tpm{$k}{$j};
		$sum_count{$k} += $count2tpm{$k}{$j} * $ensg_len{$j}{'length'};
	}
}


#####################Step3: cal TPM = gene transcript * 10^6 / sum trinscript per sample
open OUT_TPM, ">$prefix.tpm.txt" or die $!;
open OUT_FPKM, ">$prefix.fpkm.txt" or die $!;

print OUT_TPM "GID\tSymbol\tLength\tType\t" . join ("\t",@all_sample) . "\n";
print OUT_FPKM "GID\tSymbol\tLength\tType\t" . join ("\t",@all_sample) . "\n";
#for my $g (@all_gene){
for my $g (sort {$a cmp $b} keys %all_gene){
	print OUT_TPM "$g\t$ensg_len{$g}{'symbol'}\t$ensg_len{$g}{'length'}\t$ensg_len{$g}{'type'}";
	print OUT_FPKM "$g\t$ensg_len{$g}{'symbol'}\t$ensg_len{$g}{'length'}\t$ensg_len{$g}{'type'}";
	for my $s(@all_sample){
		my $tpm = sprintf("%.4f",$count2tpm{$s}{$g} * 1000000 / $sum_trans{$s});
		my $fpkm = sprintf("%.4f",$count2tpm{$s}{$g} * 1000000000 /  $sum_count{$s});
		print OUT_TPM "\t$tpm";
		print OUT_FPKM "\t$fpkm";
	}
	print OUT_TPM "\n";
	print OUT_FPKM "\n";
}
close OUT_TPM;
close OUT_FPKM;
