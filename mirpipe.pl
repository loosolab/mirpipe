#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Spec;
use File::Which;
#use File::Basename qw( dirname );

###############################################
# 
#  C o m m a n d   l i n e   o p t i o n s
#
###############################################

# help message / syntax
my $usage = "
MIRPIPE miRNA Quantification Pipeline 1.0

Usage: $0 -file <reads(.fastq/.fasta)> -ref <reference mirna(.fasta)> [options]
Example: $0 -file test.fastq -ref mirbase20/mature.fa -mismatches 3

MANDATORY:
-file <string>
Input file bearing reads (FASTA/FASTQ).

-ref <string>
FASTA file bearing reference miRNAs.

OPTIONS:
-adapter <string> []
Trim 3' adapter sequences.

-minlen <int> [18]
Minimum length of a read after trimming to be considered in the analysis.

-qval <int> [20]
Minimum phred quality for FASTQ data. Nucleotides with lower quality will be trimmed. 
This parameter is not used if FASTA formatted read data is supplied.

-maxlen <int> [28]
Maximum length of a raw read to be considered in the analysis.

-mincount <int> [5]
A read sequence must be present at least this number of times to be included.

-mismatches <int> [4]
Maximum number of mismatches allowed between reference miRNA and read sequence. This 
parameter controls the size of the miRNA clusters (more mismatches allowed = larger clusters).

-family <string> [Yes]
Cluster mirna names on the family level (mmu-miR-1-1-3p -> miR-1). Switch to \"No\" to turn off, 
e.g. in case the reference database of miRNA sequences does not obey to the naming convention of 
miRBase (<species>-miR-<#>-<suffix>).

-abundance <int> [5]
Remove read sequences from a cluster that are less than x% as abundant as the most abundant 
sequence. This is intended to suppress reads resulting from sequencing errors or biological 
miRNA variations that are expressed near the detection limit (controls the size 
of the miRNA clusters: lower minimum cluster abundance = larger clusters.
";

my($reads,$mirbase,$blastdb,$output,$output2,$output3,$tempfiles,$mismatches,$mincount,$print_log,$gzip_tempfiles,$remove_tempfiles,$nameswitch,$maxlen,$minlen,$qval,$adpt,$abundance);

GetOptions(	"file=s" => \$reads,
		"ref=s" => \$blastdb,
#		"mirbase=s" => \$mirbase,
#		"output=s" => \$output,
#		"output2=s" => \$output2,
#		"output3=s" => \$output3,
		"mismatches=i" => \$mismatches,
		"mincount=i" => \$mincount,
		"adapter=s" => \$adpt,
		"maxlen=i" => \$maxlen,
		"qval=i" => \$qval,
		"minlen=i" => \$minlen,
		"tempfiles=s" => \$tempfiles,
		"family=s" => \$nameswitch,
		"abundance=i" => \$abundance) or die("Error in command line arguments\n");


# set defaults
$mismatches = '4' if (!defined $mismatches);
$mincount = '5' if (!defined $mincount);
$adpt = 'No' if (!defined $adpt);
$maxlen = '28' if (!defined $maxlen);
$qval = '20' if (!defined $qval);
$minlen = '18' if (!defined $minlen);
$nameswitch = 'Yes' if (defined !$nameswitch);
$abundance = '5' if (!defined $abundance);
$tempfiles = 'No' if ($tempfiles eq "");

# mandatory
die "$usage" unless ($reads && $blastdb);

if ($tempfiles eq "No") {
	$remove_tempfiles = 1;
}
else {
	$remove_tempfiles = 0;
}

##############################################
#
#  Prerequisites
#
##############################################

# check for presence of prerequisites
my %paths;
$paths{cutadapt} = which('cutadapt');
$paths{makeblastdb} = which('makeblastdb');
$paths{blastn} = which('blastn');
$paths{fastq_quality_trimmer} = which('fastq_quality_trimmer');
$paths{fastq_to_fasta} = which('fastq_to_fasta');
$paths{fastx_collapser} = which('fastx_collapser');
$paths{shortreads} = which('shortreads.pl');
$paths{clusterFastaSameSeed} = which('clusterFastaSameSeed.pl');
$paths{blast2mirlist2} = which('blast2mirlist2.pl');
if ($paths{cutadapt} eq "") {
	print "The following command was not found: cutadapt. Please install cutadapt 1.4.1 or higher and include it in the system PATH if adapter trimming is desired.";
}
if ($paths{makeblastdb} eq "") {
	die "The following command was not found: makeblastdb. Please install NCBI BLAST+ 2.2.28 or higher and include it in the system PATH.";
}
if ($paths{blastn} eq "") {
	die "The following command was not found: blastn. Please install NCBI BLAST+ 2.2.28 or higher and include it in the system PATH.";
}
if ($paths{fastq_quality_trimmer} eq "") {
	die "The following command was not found: fastq_quality_trimmer. Please install FASTX-Toolkit 0.0.13 or higher and include it in the system PATH.";
}
if ($paths{fastq_to_fasta} eq "") {
	die "The following command was not found: fastq_to_fasta. Please install FASTX-Toolkit 0.0.13 or higher and include it in the system PATH.";
}
if ($paths{fastx_collapser} eq "") {
	die "The following command was not found: fastx_collapser. Please install FASTX-Toolkit 0.0.13 or higher and include it in the system PATH.";
}
if ($paths{shortreads} eq "") {
	die "The following command was not found: shortreads.pl. Please install MIRPIPE 1.0 or higher and include the directory in the system PATH.";
}
if ($paths{clusterFastaSameSeed} eq "") {
	die "The following command was not found: clusterFastaSameSeed.pl. Please install MIRPIPE 1.0 or higher and include the directory in the system PATH.";
}
if ($paths{blast2mirlist2} eq "") {
	die "The following command was not found: blast2mirlist2.pl. Please install MIRPIPE 1.0 or higher and include the directory in the system PATH.";
}



if ($blastdb ne "None") {
	print "--- Creating custom blast database\n";
	system("makeblastdb -in $blastdb -dbtype nucl -out custom");
	$mirbase = "custom";
}

##############################################
#
#  Logfile
#
##############################################

print  "Using as reads:\t$reads\nUsing as reference:\t$mirbase\nWill print logfile:\t$print_log\nWill include temporary files:\t$gzip_tempfiles\n";
print "Additional parameters:\nmax mismatches:\t$mismatches\nmin counts:\t$mincount\nmin abundance:\t$abundance\n";
print "Working on: $ENV{'PWD'}\n";

###############################################
# 
#  Run Parameters
#
###############################################

#my $dirname=File::Spec->rel2abs( __FILE__ );	#get path bearing scripts
#$dirname=dirname($dirname)."/";			#ex: /opt/mirpipe/

my $trimmed='trimmed.fastq';
my $trimmed_fasta='trimmed.fasta';
my $short='trimmed.short.fasta';
my $collapsed='trimmed.short.collapsed.fasta';
my $clustered='trimmed.short.collapsed.seedcluster';
my $cutadapted='clean.fastq';

###############################################
# 
#  Pipeline
#
###############################################

my $DATA = `head -n 1 $reads`;

for($DATA) {
	if(/^>/) {
		print "\n\n--- Read format fasta: NO quality and low length trimming\n"; 
		$trimmed_fasta=$reads;
	}
	elsif(/^@/) { 
		print "\n\n--- Read format fastq: do quality and low length trimming\n"; 
		if ($adpt ne "No") {
			system('cutadapt -m 1 -e 0.1 -f fastq -a '.$adpt.' -o '.$cutadapted.' '.$reads);
			$reads = $cutadapted;
			if (-s $reads == 0) {
				print "No reads left, will quit now.\n";
				 if ($? == -1) {
					print "failed to execute: $!\n";
				}
				exit(0);
			}
		}
		system('fastq_quality_trimmer -v -Q33 -i '.$reads.' -o '.$trimmed.' -t '.$qval.' -l '.$minlen);
		system('fastq_to_fasta -Q33 -i '.$trimmed.' -o '.$trimmed_fasta);
#		system('fastq2fasta --to-fasta '.$trimmed.' > '.$trimmed_fasta);
	}
	else {
		die("Reads not in fasta or fastq format. Please choose valid input file.");
	}
}

print "\n\n--- Remove reads longer than 28 bp\n"; 
system('shortreads.pl '.$trimmed_fasta.' '.$short.' '.$maxlen);

if (-s $short == 0) {
	print "No reads <= $maxlen left, will quit now.";
	exit(0);
}

print "\n\n--- Collapse identical reads\n"; 
system('fastx_collapser -Q33 -v -i '.$short.' -o '.$collapsed);

print "\n\n--- Cluster isomiRs and remove clusters with a read count < $mincount\n"; 
system("clusterFastaSameSeed.pl -f $collapsed -o $clustered -m $mincount");

print "\n\n--- BLASTN versus reference (ex: miRBase v20 mature sequences)\n"; 
system('blastn -query '.$clustered.'.fasta -db '.$mirbase.' -out blastout.tsv -num_alignments 15 -word_size 7 -evalue 10 -dust no -strand plus -outfmt "7 qseqid sseqid pident length qlen mismatch gaps evalue bitscore" -num_threads 8');
my $reads_mapped=`grep "# Query:" blastout.tsv | wc | tr -s \' \' | cut -d \' \' -f2`;
print "Collapsed reads with at least one alignment vs. the reference DB: $reads_mapped";


if($nameswitch eq "Yes") {
	print "\n\n--- Count reads per reference seq and cluster miRNA names on the family level (mismatches <= $mismatches, abundance >= $abundance)\n";
	system("blast2mirlist2.pl blastout.tsv $clustered.fasta $mismatches $abundance 0");
}
else {
	print "\n\n--- Count reads per reference seq without clustering miRNA names on the family level (mismatches <= $mismatches, abundance >= $abundance)\n";
	system("blast2mirlist2.pl blastout.tsv $clustered.fasta $mismatches $abundance 1");
}


if($remove_tempfiles) {
	#print "\n\n--- Gzip tempfiles\n";
	#system("tar cfvz $tempfiles blastout.tsv blastout.tsv.mirnas blastout.tsv.cluster trimmed.short.collapsed.seedcluster.fasta");
	system("rm -f blastout.tsv trimmed.short.collapsed.seedcluster.fasta trimmed.fasta trimmed.fastq trimmed.short.collapsed.seedcluster trimmed.short.fasta trimmed.short.fasta.long custom.nhr custom.nin custom.nsq trimmed.short.collapsed.fasta");
}

#rename final output
system("mv blastout.tsv.mirnas mirpipe_mirna.tsv");
system("mv blastout.tsv.cluster mirpipe_cluster.tsv");
system("mv blastout.tsv.fasta mirpipe_cluster.fasta");
