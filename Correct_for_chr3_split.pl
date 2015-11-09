#!/usr/bin/perl
use strict;
my $line;
my @temp;
use Bio::SeqIO;
use List::MoreUtils qw(uniq);
use Math::Round qw(nearest);
my $new_info;
my $counter=0;

open (INPUT, $ARGV[0]);			#GATK vcf file
while (my $liney=<INPUT>){
chomp $liney;                           
my @array=split(/\t/,$liney);
if($liney =~ /^\#\#/ && $liney !~ /ID=3B/){
print "$liney\n";
}
elsif($liney =~ /^\#\#/ && $liney =~ /ID=3B/ && $counter == 0){
print "##contig=<ID=3B,length=774434471>\n";
$counter++;
}
elsif($liney !~ /^\#\#/){
if($array[0] =~ /3B1/){
print "3B\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]\n";
}
elsif($array[0] =~ /3B2/){
my $new_pos=($array[1] + 480000000);
print "3B\t$new_pos\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]\n";
}
else{
print "$liney\n";
#if(!exists($hash{$array[0]})){
#print "$liney\n";
#$hash{$array[0]}="";
}
}
}
close INPUT;



