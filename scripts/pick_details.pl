#!usr/bin/perl
use strict;

use Data::Dumper;
my $hash={}; my $hash1={};
#-Open Bacteria Indval taxonomy file
my $file = $ARGV[0];
open (IF, $file);
while (chomp (my $line=<IF>))
{
#ASV	New_names	Group	Group	Indval	pvalue	freq	Taxonomy								
#74b7dbcd58dc83617d1a82d321245a1e	ASV_B762	6	1	0.001	11	k__Bacteria; p__Spirochaetes; c__Spirochaetes; o__Spirochaetales; f__Spirochaetaceae; g__Treponema; s__.

	my @arr = split ("\t", $line);
	$hash->{$arr[0]}=$arr[1]."\t".$arr[2]."\t".$arr[6];
}
print Dumper $hash;

#-Open Fungal Indval taxonomy file
my $file1 = $ARGV[1];
open (IF1, $file1);
while (chomp (my $line1=<IF1>))
{
#ASV	New_Names	Group	Group	Indval	pvalue	freq	Taxonomy
#473285b616b9b3435278f9a53a1fbf8a	ASV_F277	8	0.884615385	0.001	23	k__Fungi;p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified.

        my @arr1 = split ("\t", $line1);
        $hash1->{$arr1[0]}=$arr1[1]."\t".$arr1[2]."\t".$arr1[6];
}
print Dumper $hash1;

my $file2 = $ARGV[2];
open (IF2, $file2);
open (OF, ">$file2.Details.txt");
while (chomp (my $line2=<IF2>))
{
        my @arr2 = split ("\t", $line2);
	if (exists ($hash->{$arr2[0]}))
	{
		print OF $arr2[0]."\t"."$hash->{$arr2[0]}"."\t";
	}
	if (exists ($hash1->{$arr2[1]}))
        {
                print OF $arr2[1]."\t"."$hash1->{$arr2[1]}"."\t"."$arr2[2]\t$arr2[3]\n";
        }
	
}
