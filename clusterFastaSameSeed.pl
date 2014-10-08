#!/usr/bin/perl -w
# open file for reading
#
# Clustert Reads mit gleichem Anfang / 5' Ende
# Kernsequnzen mit dem größten Count werden mitgeführt
# Sequenzen werden zusammengeführt:
#  a) Sequenz 1 Teil von 2 bzw. andersrum
#  b) nur die letzte Base von Seq2 falsch ist
use strict;
use Getopt::Std;

our($opt_f, $opt_o, $opt_m);

my $usage = "Cluster reads with identical 5';only keep sequences with a minimum number of counts\nsyntax: $0 -f <fasta> -o <out> -m <minreads>\n\n";
getopt('fom');

my($fasta,$outfileprefix,$readMin);

$fasta = $opt_f // die($usage);
$outfileprefix = $opt_o // $fasta;
$readMin = $opt_m // 5;

if($fasta eq $outfileprefix) {
	$fasta =~ s/\.fasta$/\.seedCluster/g;
	$outfileprefix = $fasta;
}

my($newline,$mirID,$count,$hit,@temp,%mirhash,%cluster,%mirseq,$input,$output);
$input=0;
$output=0;
open(_F, $fasta) || die ("Could not read file $fasta\n");
while(<_F>) {
	$newline = $_;
	if($newline =~ m/^\>/){
		$input+=1;
		$newline =~ s/\>//;
		chomp($newline);
		($mirID,$count)=split("-",$newline);
		if($count >= $readMin){
			$mirhash{$mirID}{"count"}=$count;
			$mirhash{$mirID}{"mirID"}=$mirID;
		}
	}else{
		if($count >= $readMin){
			chomp($newline);
			$mirhash{$mirID}{"seq"}=$newline;
		}
	}
}
close (_F);


#$fasta =~ s/\.fasta$/\.seedCluster/;
open(OUTF, ">".$outfileprefix) || die ("Could not read file $outfileprefix\n");

for $mirID ( sort { $mirhash{$b}{"count"} <=> $mirhash{$a}{"count"} } keys %mirhash ) {
	push(@temp,$mirhash{$mirID});
}

if(!@temp) {
	print("No miRNA's left for clustering... ");
	exit(0);
} 

$cluster{$temp[0]{"seq"}}{"count"} = $temp[0]{"count"};
$cluster{$temp[0]{"seq"}}{"mirID"} = $temp[0]{"mirID"};
$cluster{$temp[0]{"seq"}}{"mitglieder"} = ">".$temp[0]{"count"}."\t".$temp[0]{"seq"}."\n";
$cluster{$temp[0]{"seq"}}{"mitgliedzahl"} = 1;

for(my $i=0;$i <= $#temp;$i++ ){
	
	my $seq1 = $temp[$i]{"seq"};
	if ($temp[$i]{"count"} == 0){next;}
	
	for(my $j=$i+1;$j <= $#temp;$j++ ){
	
		my $seq2 = $temp[$j]{"seq"};
		if ($temp[$j]{"count"} == 0){
			next;
		}
		$hit=0;
		
		if(length($seq1)>length($seq2)){
			if(index($seq1,$seq2) == 0){
				$hit=1;
			}
		}else{
			if(index($seq2,$seq1) == 0 || index($seq1,substr($seq2,0,length($seq2)-1)) == 0){
				$hit=1;	
			} 
		}
		
		if($hit==1){
			$cluster{$seq1}{"count"} = $cluster{$seq1}{"count"} + $temp[$j]{"count"};
			$cluster{$seq1}{"mitgliedzahl"}++;
			$cluster{$seq1}{"mitglieder"} = $cluster{$seq1}{"mitglieder"}.">".$temp[$j]{"count"}."\t".$temp[$j]{"seq"}."\n";
			$temp[$j]{"count"}=0; # "deletes" element
			if(defined $cluster{$seq2}){delete $cluster{$seq2};}
		}else{
			$cluster{$seq2}{"count"} = $temp[$j]{"count"};
			$cluster{$seq2}{"mirID"} = $temp[$j]{"mirID"};
			$cluster{$seq2}{"mitglieder"} = ">".$temp[$j]{"count"}."\t".$temp[$j]{"seq"}."\n";
			$cluster{$seq2}{"mitgliedzahl"} = 1;
		}
		
		
	}

	
}

for my $seq(sort { $cluster{$a}{"count"} <=> $cluster{$b}{"count"} } keys %cluster ){
	print OUTF ">",$cluster{$seq}{"mirID"},":",$cluster{$seq}{"count"},"\n$seq\nCluster:\n",$cluster{$seq}{"mitglieder"},"\n";
	$output+=1;
}
close (OUTF);

#$fasta =~ s/\.seedCluster$/\.seedCluster\.fasta/;
open(OUTF, ">".$outfileprefix.".fasta") || die ("Could not read file $fasta\n");
for my $seq(sort { $cluster{$a}{"count"} <=> $cluster{$b}{"count"} } keys %cluster ){
	print OUTF ">",$cluster{$seq}{"mirID"},"-",$cluster{$seq}{"count"},"\n$seq\n"; 
}
close (OUTF);

print "Input reads for isomiR clustering: ".$input."\n";
print "Output reads after isomiR clustering and above minimum count: ".$output."\n";
