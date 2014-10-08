mirpipe
=======

Source code for the MIRPIPE miRNA detection and quantification pipeline.

Introduction
============

A host of applications exist in the field of RNASeq miRNA analyses, which attempt to quantify transcripts based on reference genomes or reference miRNAs, using different degrees of granularity including subsequences and precursors of miRNAs, deviating handling of technical and biological variation, offering a range of interfaces to the user (web, Linux command line, graphical user interfaces), and bearing variable demands considering user expertise in informatics and biological interpretation.

MIRPIPE offers the following unique combination of features:
- Robust identification of miRNAs able to deal with comparisons spanning multiple species necessary for cases without a reference genome or previously identified miRNAs in the examined organism
- Proper handling of biological variation (isomiRs, cluster of related miRNAs)
- Proper handling of technical variation (i.e. sequencing errors) introduced by 454 or IonTorrent sequencers
- Output of one count value per identified miRNA or miRNA family able to directly serve as input for differential expression analyses
- Interfaces for both use cases: quick online analyses as well as inclusion in Linux workflow scripts
- Freeware and Open Source

MIRPIPE is specifically tuned to be robust versus technical variation (i.e. errors) introduced by sequencing technologies like 454 or IonTorrent, as well as biological variation due to sequence differences originating from the evolutionary distance between query and reference species. These effects are addressed on multiple levels (isomiR handling, minimum read copy number, clustering which removes comparatively low abundance reads, centering on the miRNA family). MIRPIPE is open source and runs fully automatic offering both, a public web interface and a Linux console version that can be integrated into existing workflows.

System requirements and installation
====================================

MIRPIPE works on current PC hardware without the need for server-grade amounts of RAM or CPU.

## Software requirements

- Linux
- Perl 5.10
- Perl Graph module (specifically Graph::Undirected)
- FASTX-Toolkit 0.0.13
- BLASTN 2.2.28+
- Cutadapt 1.4.1 (only if adapter removal is used)

All of these tools should be included in the PATH variable to permit system-wide access for the user executing MIRPIPE.

## Installation

- Unzip scripts: tar -zxvf mirpipe.tar.gz
- Add directory to $path variable: export PATH=$PATH:/current/mirpipe/dir

## Testing the installation

Switch to the test subdirectory included in mirpipe and run ./test.sh. If you see a message indicating success, your installation is ready to go. Otherwise, missing dependencies will be shown as error messages.

Usage
=====

mirpipe.pl -file <reads(.fastq/.fasta)> -ref <reference mirna(.fasta)> [options]
Running MIRPIPEs main script mirpipe.pl without any parameters will print out available options.

Contribution and contact
========================

## How to contribute?
Feel free to fork this repository, add your changes and send a pull request.

## Contributors

- Carsten Künne (carsten.kuenne@mpi-bn.mpg.de)
- Jens Preußner (jens.preussner@mpi-bn.mpg.de)
- Mario Looso (mario.looso@mpi-bn.mpg.de)
- Mario Herzog

## Citation

Kuenne C, Preussner J, Herzog M, Braun T and Looso M: MIRPIPE - quantification of microRNAs in niche model organisms. Bioinformatics (2014), doi: 10.1093/bioinformatics/btu573

