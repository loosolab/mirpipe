#!/usr/bin/perl -w
# open file for reading
#
use strict;
my($file,$outfile,$outfile_long,$fh1,$fh2,$newline,$name,$sequence,$entry,$cnt1,$cnt2,%short,%middle,%long,$maxlen);
$file=$ARGV[0];
$outfile=$ARGV[1];
$maxlen=$ARGV[2];
$outfile_long="$outfile.long";
$entry=1;
$cnt1=0;
$cnt2=0;

open(INFILE, $file) || die ("could not read file $file\n");
open($fh1, ">", $outfile) or die "Couldn't open $outfile: $!";
open($fh2, ">", $outfile_long) or die "Couldn't open $outfile_long: $!";
while(<INFILE>) {
	$newline=$_;
	if ($entry){
		$name=$newline;
		$entry=0;
	} else{
		$sequence=$newline;
		if((length($sequence)-1) <= $maxlen){
			print $fh1 $name,$sequence;
			$cnt1++;
		}elsif((length($sequence)-1) > $maxlen){
			print $fh2 $name,$sequence;			
			$cnt2++;
		}
		$entry=1;
	}
}
close (INFILE);		
close($fh1);
close($fh2);

print "Total reads: ".($cnt1+$cnt2)."\n";
print "Reads <= max bp ($maxlen): $cnt1\n";
print "Reads > max bp ($maxlen): $cnt2\n";
exit(0);
