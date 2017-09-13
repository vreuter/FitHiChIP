					FitHiChIP 
					---------

Developed by
Dr. Sourya Bhattacharyya 

Supervisor: Dr. Ferhat Ay
Vijay-Ay lab
La Jolla Institute for Allergy and Immunology
San Diego, CA 92037, USA

******************************
For a detailed documentation, please check Manual.pdf within this repository.
******************************

Introduction
-------------

FitHiChIP computes statistically significant interactions between different pairs of genomic intervals of a given 
HiChIP or PLAC-seq data. Currently `cis' interactions are only supported. 

\section{Installation}

Requires a Linux environment (with bash) with Python (development version: 2.7.13), R (development version: 3.3.3) and 
perl (default version) installed. Also requires following packages: 

	1) Python packages: pypairix (0.1.7) (https://github.com/4dn-dcic/pairix), pysam (0.10.0) (https://github.com/pysam-developers/pysam)
	2) R libraries splines, fdrtool, ggplot2, optparse, GenomicRanges, Parallel
	3) Package MACS2 for peak calling (https://github.com/taoliu/MACS)

In addition, creation of paired end alignment files (in either BAM or PAIRIX formats) require the following packages: 

1) samtools (http://www.htslib.org/doc/samtools.html), 2) BWA (http://bio-bwa.sourceforge.net/), 3) Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), 4) picard tool (http://broadinstitute.github.io/picard/), 5) pairix (https://github.com/4dn-dcic/pairix), 6) bam2pairs (https://github.com/4dn-dcic/pairix)


Files and Directories
----------------------

1) FitHiChIP.sh: main executable
2) Preprocess: folder containing scripts and functions for creating paired end alignment from input fastq files (for details and corresponding commands, users may refer to the manual).
	A) Prep_PLAC.sh: Alignment files useful for PLAC seq pipeline
	B) Prep.sh: Alignment file for HiChIP data
	C) Align.sh: paired end alignment in BAM and PAIRIX format.
3) src: Few source codes.
4) TestData: sample dataset (1 pairix and 1 bam formatted)
5) scripts: Sample scripts to invoke the datasets provided.


Types of interactions
---------------------

1) ALL to ALL: Interactions between every pairs of bins. User does not need to 
provide any peak detection file.
2) Peak to Peak: Possible when user provides a separate peak detection output file (using the 
-P option). Otherwise, for a BAM alignment file, FitHiChIP computes the peaks using MACS2. 
Interactions between any pairs of fixed size peak segments are computed.
3) Peak to Non Peak: Interaction between a peak segment (fixed size 
bin overlapping with a peak region) and a non-peak segment (fixed size bin 
which does not overlap with any peaks).
4) Peak to ALL: Encapsulates both peak to peak and peak to non peak interactions.

Command line options
--------------------

For details, and example commands, please refer to the manual.

FitHiChIP.sh -I inpfile -M method -o outdir -n prefix [-P peakfile] [OPTIONS]

inpfile: input alignment file (BAM or PAIRIX)
method: 3 (default - ALL to ALL interactions) / 2 (peak to non peak) / 1 (peak to peak) / 0 (peak to peak)
outdir: Output directory for results.
prefix: Prefix string of output file names.
peakfile: peak detection output file (required for methods 0, 1, or 2).
OPTIONS:
	-b Binsize: Size of a bin in bp (default = 5000, means 5 Kb bins).
	-t Threads: Number of threads (default 1).
	-L LowDistThres: Lower distance threshold of interaction (in terms of bp). Default 20000 (20 Kb).
 	-U UppDistThres: Upper distance threshold of interaction (in terms of bp). Default 2000000 (2 Mb).
 	-f FitHiCBinMethod: 1 (default and recommended) or 0. If 1, equal occupancy (contact count) bins are employed for FitHiC. Else equal length bins are used.
 	-N NBins: Max no of bins (equal occupancy or equal length) employed in FitHiC. Default 200.
 	-B BinomDistr: 0 (default) or 1. If 1, binomial distribution model between the observed genomic distance and the contact count is also computed (Duan et. al. 2010) and compared with the FitHiC output.
 	-q QVALUE: Minimum FDR (q-value) cutoff for interaction detection [default = 0.01].
 	-v verbose: 1 or 0 (default); if 1, time log is printed along with the output.
 	-D DrawFig: 1 or 0 (default); if 1, various analysis plots regarding the performance of FitHiChIP are generated.
 	-g GSIZE: If MACS2 is used for peak detection, its genome size parameter. Default 'hs'.


Example commands
----------------

1) FitHiChIP.sh -I LongAlign.bam -P ShortPeakFile.Peaks -M 0 -o /home/sourya/FitHiChIP/ -n 'PREFIX' -b 5000 -t 4 -L 20000  -U 2000000

2) FitHiChIP.sh -I inp.bam -M 3 -o /home/sourya/FitHiChIP/ -n 'tempFitHiChIP' -b 5000  -t 8

3) FitHiChIP.sh -I inp.pairix -P sample.peak.bed -M 2 -o /home/sourya/FitHiChIP/ -n 'tempFitHiChIP' -b 5000 -t 8 -B 1

4) FitHiChIP.sh -I inp.bam -M 1 -o /home/sourya/FitHiChIP/ -n 'tempFitHiChIP' -b 5000  -t 8 -L 10000 -U 3000000


Contact
---------

For any queries, please e-mail:

Sourya Bhattacharyya (sourya@lji.org)
Ferhat Ay (ferhatay@lji.org)





