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
HiChIP or PLAC-seq data. Currently `cis' interactions are only supported. Requires a Linux environment (with bash) 
with Python (development version: 2.7.13), R (development version: 3.3.3) and 
perl (default version) installed. 


Files and Directories
----------------------

1) FitHiChIP.sh: main executable

2) sample_script.sh: An example of invoking FitHiChIP.sh with a sample configuration file.

3) configfile_*: various configuration files (samples - not to be directly used) for processing BAM, PAIRIX or HiC-pro pipeline 
generated validpairs file are provided for example.

4) Preprocess: folder containing scripts and functions for creating paired end alignment from input fastq files (for details and corresponding commands, users may refer to the manual).
	
	A) Prep_PLAC.sh: Creating alignment files (in BAM or PAIRIX formats) and corresponding ChIP seq peaks, from input 
	fastq reads (PLAC-seq data).

	B) Preprocess_HiChIP_fastq.sh: Creating alignment files (in BAM or PAIRIX formats) and corresponding ChIP seq peaks, 
	from input fastq reads (HiChIP data)
	
	C) Preprocess_HiChIP_BAM.sh: It uses HiChIP paired end alignment, applies HiChIP specific strand filtering and computes 
	ChIP seq peaks.


5) src: Directory containing source codes.

6) bin: associated executables for executing FitHiChIP pipeline on different input formats.


*** For different types of interactions, and command line options, please refer to the manual.


Contact
---------

For any queries, please e-mail:

Sourya Bhattacharyya (sourya@lji.org)

Ferhat Ay (ferhatay@lji.org)





