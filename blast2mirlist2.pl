#!/usr/bin/env perl -w
# version	:1.2
# author	:Jens Preussner<jens.preussner@mpi-bn.mpg.de>
# 
# date		:2014-10-07

use strict;
use Graph::Undirected;
use List::Util qw(max sum);

# ------------------------------------------------------------------
# GLOBAL PARAMETERS.
# ------------------------------------------------------------------
my $delimitor=", ";		#delimitor to place between clustered mirs of equal score (ambiguous) if no unique hit could be found
my $nodup = 0;			#set to 1: do not cluster ambiguous mirs
my $verbose = 1; #Print status messages at the end of each code block
my $debug = 0; #Print status messages from within blocks (large output)

# ------------------------------------------------------------------
# USER SPECIFIED PARAMETERS.
# ------------------------------------------------------------------
my $blastFile = $ARGV[0];		#blast results
my $fastaFile = $ARGV[1];		#reads
my $maxMismatches = $ARGV[2];	#maximum mismatches
my $minAbundance = $ARGV[3];		#Remove read sequences from a cluster that are less than x% as abundant as the most abundant sequence
my $nofamily =  (defined($ARGV[4]) ? $ARGV[4] : 0);		#if == 1:  do not collapse mirs to family level (mmu-miR-1-1-3p -> miR-1)

#print "Remove read sequences from a cluster that are less than $abundance % as abundant as the most abundant sequence in that cluster\n" if $verbose;
$minAbundance = $minAbundance * 0.01 if $minAbundance > 0;

# ------------------------------------------------------------------
# USEFUL FUNCTIONS.
# ------------------------------------------------------------------
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

# ------------------------------------------------------------------
# READ IN FASTA FILE.
# As a result, all sequences are stored in the %reads hash
# ------------------------------------------------------------------
my %reads = ();
open(_F, $fastaFile) || die ("Could not read file $fastaFile\n");
while(<_F>) {
	if($_ =~ m/^\>/){
		$_ =~ s/\>//;
		chomp($_);
		my ($name,$count) = split("-",$_);
		print "Read $name with $count counts\n" if $debug;
		$_ = <_F>;
		chomp($_);
		$reads{$name}{"seq"}=$_;
		$reads{$name}{"count"} = $count;
	}
}
close (_F);
print "Read ". keys(%reads) . " sequences in total.\n" if $verbose;

# ------------------------------------------------------------------
# READ IN BLAST RESULTS.
# As a result, an undirected graph is build containing the information
# which microRNAs share a read (are hit equally good by the same read).
# Additionally, per microRNA, the reads are stored in a hash of arrays.
# As a third result, the read are annotated with their best matches.
# ------------------------------------------------------------------
my $graph = Graph::Undirected->new;
my %mirHits = ();
my %graphEdgesEvidence = ();
my %graphReadsEvidence = ();

# Temporary variables, only used in this block
my $blockQuery=""; # stores the query of the actual blast hit
my %blockMinMismatches=(); # counts how often each subject is hit by the query
my $foundAHit = 0; # A flag that indicates if a certain query has a matching microRNA.
my $firstRun = 1;
my $bestScore = $maxMismatches;

open(_F, $blastFile) || die ("could not read file $blastFile\n");
while(<_F>) {
	if($_ =~ m/^\# BLAST/) { #tabular blast new record query block OR final line of file, initialize all parameters
		if ( $firstRun ) {
			$firstRun = 0;
			next;
		}
		if( $foundAHit ) {
			my($name,$count) = split("-",$blockQuery);
			# ------------------------------------------------------
			# A hit is found. We sort the microRNAs for homology by
			# ascending number of mismatches versus the query.
			# If two microRNAs have the same score, the will be linked
			# in the undirected graph.
			# ------------------------------------------------------
			my @sorted = (sort { $blockMinMismatches{$a} <=> $blockMinMismatches{$b} } keys %blockMinMismatches);
			my $i = 0;
			my $bestMir = "";
			foreach (@sorted){
				if ($i == 0) {
					$graph->add_vertex($_);
					$bestMir = $_;
					print "Added to graph: best microRNA for read $name is $_.\n" if $debug;
					# ------------------------------------------------------
					# Here we will save the query to the microRNA
					# ------------------------------------------------------
					push(@{ $mirHits{$_} }, $name) if (exists($mirHits{$_}) and !($name ~~ @{ $mirHits{$_} }));
					$mirHits{$_} = [ $name ] if (!exists($mirHits{$_}));
					push(@{ $reads{$name}{"annotation"} }, $bestMir) if (exists($reads{$name}{"annotation"}) and !($bestMir ~~ @{ $reads{$name}{"annotation"} }));
					$reads{$name}{"annotation"} = [$bestMir] if (!exists($reads{$name}{"annotation"}));
					print "Annotation: Added $_ to read $name.\n" if $debug;
				}
				elsif ($blockMinMismatches{$_} <= $blockMinMismatches{$bestMir}) {
					$graph->add_vertex($_);
					$graph->add_edge($_,$bestMir);
					# ------------------------------------------------------
					# We will keep track of the edges, to have the possibility
					# to remove them later during graph correction.
					# ------------------------------------------------------
					my $edge = $_.":".$bestMir;
					my $edge_return = $bestMir.":".$_;
					push(@{ $graphEdgesEvidence{$name} }, $edge) if (exists($graphEdgesEvidence{$name}));
					$graphEdgesEvidence{$name} = [ $edge ] if (!exists($graphEdgesEvidence{$name}));
					push(@{ $graphEdgesEvidence{$name} }, $edge_return) if (exists($graphEdgesEvidence{$name}));
					$graphEdgesEvidence{$name} = [ $edge_return ] if (!exists($graphEdgesEvidence{$name}));
					push(@{ $graphReadsEvidence{$edge} }, $name) if (exists($graphReadsEvidence{$edge}));
					$graphReadsEvidence{$edge} = [ $name ] if (!exists($graphReadsEvidence{$edge}));
					push(@{ $graphReadsEvidence{$edge_return} }, $name) if (exists($graphReadsEvidence{$edge_return}));
					$graphReadsEvidence{$edge_return} = [ $name ] if (!exists($graphReadsEvidence{$edge_return}));
					print "Added to graph: equal best microRNA for read $name is $_.\n" if $debug;
					# ------------------------------------------------------
					# Here we will save the query to the microRNA
					# ------------------------------------------------------
					push(@{ $mirHits{$_} }, $name) if (exists($mirHits{$_}) and !($name ~~ @{ $mirHits{$_} }));
					$mirHits{$_} = [ $name ] if (!exists($mirHits{$_}));
					push(@{ $reads{$name}{"annotation"} }, $_) if !($_ ~~ @{ $reads{$name}{"annotation"} });
					print "Annotation: Added $_ to read $name.\n" if $debug;
				}
				$i = 1;
			}
		}
		$bestScore = $maxMismatches;
		$blockQuery="";
		%blockMinMismatches=();
		$foundAHit = 0;
	}
	elsif($_ =~ m/^\d/ ) {
		# ------------------------------------------------------
		# The blast result line is split into single fields.
		# If the calculated score is better than the previous one,
		# we save the score of the current block.
		# Additionally we add 1 to the count of the microRNA
		# ------------------------------------------------------ 
		my($queryID, $subjectID, $identity, $alignmentLength, $queryLength, $mismatches, $gaps, $evalue, $bitScore) = split("\t",$_);
		my $score = abs($queryLength-$alignmentLength) + $mismatches + $gaps;
		if( ($score <= $bestScore) ) { #if less mismatches than previous best hit
			if ($nofamily != 1) { #collapse mir name to family level
				$subjectID =~ s/^\w{3,}\-//;					#hit mirna: remove first 3 characters: dre-miR-203b-1-3p -> miR-203b-1-3p
				$subjectID =~ s/-(5|3)p(\.[0-9])?$//;				#hit mirna: remove final -3p/-5p: 
				my @n=split("-",$subjectID);					#count number of "-" left
				$subjectID =~ s/-(\d+)$// if ($#n == 2); #if == 2: remove optional precursor number, remove final "-" followed by any number
			}
			$blockQuery = $queryID;
			if(exists($blockMinMismatches{$subjectID})) {
				$blockMinMismatches{$subjectID} = $score if $score < $blockMinMismatches{$subjectID};
			}
			else {
				$blockMinMismatches{$subjectID} = $score;
			}
			$bestScore = $score;
			$foundAHit = 1;
		}
	}
}
close (_F);

print keys(%mirHits) . " microRNAs were hit by all reads in total.\n" if $verbose;

# ------------------------------------------------------------------
# VARIABLE SUMMARY UP TO NOW
# Hash %reads : Contains information on all reads.
# $reads{$key}{"seq"} : Nucleotide sequence of read $key
# $reads{$key}{"count"} : Read count
# $reads{$key}{"annotation"} : (array) All microRNAs that matched equally well to the read $key
# Graph $graph : An undirected graph with information on which microRNAs share a read.
# Hash %mirHits : Contains information on which microRNA is connected to which reads
# $mirHits{$key} : (array) All reads that hit the microRNA $key.
# Hash %graphEdgesEvidence: Contains information on which read is responsible for which edge.
# $graphEdgesEvidence{$key} : (array) All edges that are caused by read $key.
# Hash %graphReadsEvidence: Contains information on which edge is caused by which reads.
# $graphReadsEvidence{$key} : (array) All reads that cause the edge $key.
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# GRAPH ANALYSIS
# We define a threshold for each connected component in the graph.
# Based on this threshold, we delete read-microRNA connections
# if the read count is lower than this threshold.
# ------------------------------------------------------------------
my @cc = $graph->connected_components();

foreach my $component (@cc) {
	my @componentCounts = ();
	foreach my $vertex (@$component) {
		foreach my $k (@{ $mirHits{$vertex} }) {
			push(@componentCounts, $reads{$k}{"count"})
		}
	}
	my $threshold = (max @componentCounts) * $minAbundance;
	print "Threshold for $component is $threshold.\n" if $debug;
	# ------------------------------------------------------
	# After the threshold is defined, we sort out all reads
	# with counts lower than the threshold.
	# We might end up with microRNAs without read support,
	# i.e. vertices without links to reads, so those are
	# deleted immediately.
	# Additionally, we might end up with edges that have to
	# be deleted, because their supporting read was deleted.
	# ------------------------------------------------------
	foreach my $vertex (@$component) {
		my @deletedReads = grep { $reads{$_}{"count"} <= $threshold } @{ $mirHits{$vertex} };
		@{ $mirHits{$vertex} } = grep { $reads{$_}{"count"} > $threshold } @{ $mirHits{$vertex} };
		# remove vertices without read support immediately
		$graph->delete_vertex($vertex) if (@{ $mirHits{$vertex} } == 0);
		# ------------------------------------------------------
		# We will removed edges, that are not supported any more
		# by looking at all reads with counts below the threshold.
		# If edges were introduced by those reads, the read evidence
		# is deleted. If an edge remains without further read evidence,
		# it is ultimately removed from the graph, creating two new
		# connected components.
		# ------------------------------------------------------
		print "----------------------------\nWorking on $vertex:\n" if $debug;
		foreach my $r (@deletedReads) {
			print "Read $r will be deleted in ...\n" if $debug;
			if(exists($graphEdgesEvidence{$r})) {
				my @edges = uniq(@{ $graphEdgesEvidence{$r} });
				foreach my $e (@edges) {
					my @fromTo = split(":",$e); 
					my $e_return = $fromTo[1].":".$fromTo[0];
					print "edge $e and $e_return: " if $debug;
					# Delete the read evidence from the edge.
					@{ $graphReadsEvidence{$e} } = grep {$_ ne $r} @{ $graphReadsEvidence{$e} };
					@{ $graphReadsEvidence{$e_return} } = grep {$_ ne $e_return} @{ $graphReadsEvidence{$e_return} };
					print @{ $graphReadsEvidence{$e} }." and ".@{ $graphReadsEvidence{$e_return} }." read evidences left..\n" if $debug;
					if (@{ $graphReadsEvidence{$e} } == 0 and @{ $graphReadsEvidence{$e_return} } == 0) {
						$graph->delete_edge($fromTo[0],$fromTo[1]);
						$graph->delete_edge($fromTo[1],$fromTo[0]);
					}
				}
			}
		}
	}
}

# ------------------------------------------------------
# OUTPUT PREPARATION
# We need to recalculate the connected components, since
# we may have delete some vertices in the previous step.
# ------------------------------------------------------
my @fastaOutput = ();
my @clusterOutput = ();
my @mirlistOutput = ();
@cc = $graph->connected_components();
foreach my $component (@cc) {
	my $componentID = "";
	my $componentMembers = join(",",@$component);
	my @componentReads = ();
	foreach my $vertex (@$component) {
		my $mostabundantSequence = "";
		my $mostabundantCount = 0;
		my $percentAmbiguity = 0;
		# ------------------------------------------------------
		# PREPARE OUTPUT 1: FASTA FILE
		# All reads that are still part of a component, their count
		# and their annotation
		# ------------------------------------------------------
		push(@fastaOutput, @{ $mirHits{$vertex} });
		# ------------------------------------------------------
		# PREPARE OUTPUT 2: CLUSTER LIST
		# The component ID, the reads sequence, the reads count and its annotation. 
		# ------------------------------------------------------
		$componentID = $graph->connected_component_by_vertex($vertex);
		foreach my $k (@{ $mirHits{$vertex} }) {
			push(@componentReads, $k) unless ($k ~~ @componentReads);
			if($reads{$k}{"count"} > $mostabundantCount) {
				$mostabundantCount = $reads{$k}{"count"};
				$mostabundantSequence = $reads{$k}{"seq"};
			}
		}
		# ------------------------------------------------------
		# PREPARE OUTPUT 3: MICRORNA LIST
		# The microRNA, its total count, percentage of ambigous
		# reads, the component ID, the most abundant sequence,
		# counts of the most abundant sequence and a list of all
		# cluster members.
		# ------------------------------------------------------
		# For ambigous reads, we select all counts of those reads that have one or more commata in their annotation
		my @ambigousReads = map { @{ $reads{$_}{"annotation"} } > 1 ? $reads{$_}{"count"} : 0 } @{ $mirHits{$vertex} };
		my $ambigousSum = sum @ambigousReads;
		my @vertexCounts = map {$reads{$_}{"count"}} @{ $mirHits{$vertex} };
		my $vertexSum = sum @vertexCounts;
		push(@mirlistOutput,[$vertex, $vertexSum, sprintf("%.2f",$ambigousSum/$vertexSum), $componentID, $mostabundantSequence, $mostabundantCount, $componentMembers])
	}
	foreach my $r (@componentReads) {
		push(@clusterOutput,[$componentID, $reads{$r}{"seq"}, $reads{$r}{"count"}, $componentMembers, join(",",@{ $reads{$r}{"annotation"} })]);
	}
}

# ------------------------------------------------------
# OUTPUT 1: FASTA FILE
# All reads that are still part of a component, their count
# and their annotation
# ------------------------------------------------------

@fastaOutput = uniq(@fastaOutput);
print "After cleaning, $#fastaOutput reads are left for output.\n" if $verbose;

open(OUTF, ">".$blastFile.".fasta") || die ("Could not read file $blastFile.fasta\n");
foreach my $f (@fastaOutput) {
	print OUTF ">".join(",",@{$reads{$f}{"annotation"}})." count=".$reads{$f}{"count"}."\n".$reads{$f}{"seq"}."\n";
}
close(OUTF);

# ------------------------------------------------------
# OUTPUT 2: CLUSTER LIST
# The component ID, the reads sequence, the reads count and its annotation. 
# ------------------------------------------------------
@clusterOutput = uniq(@clusterOutput);

my @clusterOutput_sorted = sort {
	$a->[0] <=> $b->[0] || 
	$b->[2] <=> $a->[2] 
} @clusterOutput;

print "After cleaning, $#clusterOutput_sorted microRNA clusters are left.\n" if $verbose;

open(OUTC, ">".$blastFile.".cluster") || die ("Could not read file $blastFile.cluster\n");
foreach (@clusterOutput_sorted) {
	my $line = join("\t",@{$_});
	print OUTC $line."\n";
}
close (OUTC);

# ------------------------------------------------------
# OUTPUT 3: MICRORNA LIST
# ------------------------------------------------------
my @mirlistOutput_sorted = sort {
	$a->[3] <=> $b->[3] || 
	$b->[1] <=> $a->[1] 
} @mirlistOutput;

open(OUTM, ">".$blastFile.".mirnas") || die ("Could not read file $blastFile.mirnas\n");
foreach (@mirlistOutput_sorted) {
	my $line = join("\t",@{$_});
	print OUTM $line."\n";
}
close(OUTM);

print "Done.\n" if $verbose;
