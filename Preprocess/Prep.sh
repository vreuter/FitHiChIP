#!/bin/bash -ex
#PBS -l nodes=1:ppn=4
#PBS -l mem=10GB
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

#===============
# Code to accompany FitHiChIP
# preprocessing of the input fastq / bam files to produce a sorted alignment file with duplicate removal
# Also supports restriction enzyme based cutting of the input paired end segments
# Note: the program uses BWA aligner
# 	Input: fastq files or already an aligned bam file
# 	Optional input: bed file storing the restriction enzyme cuts
# 	Output: sorted + duplicate removed alignment file (bam format)
# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# usage info
usage(){
cat << EOF

usage: 

A) When fastq reads are provided as the input:
./Prep.sh [-h] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-G BWA_GENOME] [-c CUT_ENZ] [-d outdir] [-p Picard_Exec] [-q MAPQ_Thr] [-i INSERT_SIZE] [-m MAX_MEM] [-t THREADS]

B) When an existing alignment is provided as the input:
./Prep.sh [-h] [-X Aligned_Read] [-n PREFIX] [-c CUT_ENZ] [-d outdir] [-p Picard_Exec] [-q MAPQ_Thr] [-i INSERT_SIZE] [-m MAX_MEM] [-t THREADS]

Example:

A) Given input fasta files and restriction enzyme file, and a quality threshold of 20:
./Prep.sh -f R1.fq.gz -r R2.fq.gz -n 'demo_prefix' -G 'bwa_ref_genome.fa' -c RE_cutter.sites.bed -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar -q 20

B) Given input fasta files but no restriction enzyme file (and a default quality threshold):
./Prep.sh -f R1.fq.gz -r R2.fq.gz -n 'demo_prefix' -G 'bwa_ref_genome.fa' -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar

C) Given input alignment file (BAM format) - here user needs to apply the restriction enzyme (if required) before providing the bam file as the input
./Prep.sh -X inp_align.bwa -n 'demo_prefix' -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar

Output:

A sorted and filtered alignment file with duplicate segments removed

Options:    

  -- required:
	-f  FASTQ1           R1 of pair-end sequencing data  [.fq|.gz|.bz2].
	-r  FASTQ2           R2 of pair-end sequencing data [.fq|.gz|.bz2].
	-X  Aligned_Cut 	 Paired alignment of the reads (if user has the alignment file). 
	-G  BWA_GENOME       BWA indexed reference genome (if the fasta reads are provided as inputs) (mandatory parameter)
	-p  picard_exec      Location of picard executable (.jar file) (mandatory parameter)
	-d  OutDir 			 Set the output directory which will store all the results [Default: current working directory]
	-n  PREFIX           Prefix of output file names [Default: empty string].	
  -- optional:
  	-c  CUT_ENZ          File with the restriction cut site (bed formatted)
  	-q  MAPQ_Thr		 Threshold of mapping quality [Default = 30]
	-i  INSERT_SIZE      Insert size cutoff for long-range pairs [Default = 10000]. It should be at least 1000.
	-t  INT              Set number of sorting, BWA mapping threads [Default = 1].
	-m  MAX_MEM          Set max memory of duplication removal [Default = 1G].
	
EOF
}

# default parameters
MIN_INSERT_SIZE=10000
OutDir=`pwd`	#'./'
# the number of threads used for execution
THREADS=1
# maximum memory to be used
MAX_MEM="1G"
# default input fastq files
FASTQ1=""
FASTQ2=""
# default input alignment file (paired end) with RE based cutting
Aligned_Cut=""
# default value of the restriction enzyme cutter file
# for HiChIP experiment, no cutter site is provided at all
# in such cases, aligned files (from the Fastq files) are straightaway applied duplicate removal technique
CUT_ENZ=""
# this is the default quality threshold employed for each alignment
MAPQ_Thr=30
# default prefix 
PREFIX=""

# # default executable of the picard tool
# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

while getopts "f:r:n:G:c:i:d:t:m:p:X:q:" opt;
do
	case "$opt" in
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		n) PREFIX=$OPTARG;;
		G) GENOME=$OPTARG;;
		c) CUT_ENZ=$OPTARG;;
		i) MIN_INSERT_SIZE=$OPTARG;;
		d) OutDir=$OPTARG;;
		t) THREADS=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		p) picard_exec=$OPTARG;;
		X) Aligned_Cut=$OPTARG;;
		q) MAPQ_Thr=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [ $MIN_INSERT_SIZE < 1000 ]; then 
	echo "error: -i '$MIN_INSERT_SIZE' could not smaller than 1000."; 
	exit 1;
fi

if [[ -z $FASTQ1 ]] || [[ -z $FASTQ2 ]]; then
	if [[ -z $Aligned_Cut ]]; then
		echo 'User neither provided an input alignment file, nor provided two fastq reads - quit !!'
		exit 1
	fi
fi

if [[ ! -z $FASTQ1 ]] && [[ ! -z $FASTQ2 ]]; then
	if [[ -z $GENOME ]]; then
		echo 'User provided two fastq reads but did not provide the reference BWA genome for alignment - quit !!'
		exit 1
	fi
fi

if [[ -z $picard_exec ]]; then
	echo 'User did not provide the absolute path of Picard executable (for duplicate removal) - quit !! '
	exit 1
fi

# output directory for storing the results
mkdir -p $OutDir
echo '**** Output directory for storing the results: '$OutDir

# this is a log file which stores all the performance measures
# and obsrvations
logfile=$OutDir/$PREFIX.long_short_int.logfile

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)

cd $script_dir
#----------------------------------

# output directory containing the results of alignment
BWADir=$OutDir'/Alignment_MAPQ'$MAPQ_Thr
mkdir -p $BWADir
echo "Directory containing the alignment: "$BWADir

#===========================================
# the first condition when the input fasta files are provided 
# perform the alignment using the BWA aligner
if [[ ! -z $FASTQ1 ]] && [[ ! -z $FASTQ2 ]]; then
	#===================================================
	# align the fasta files, to form the paired end read alignment
	# Note: Here we do not remove information regarding '/chrM/d;/random/d;/chrUn/d' since we found that sometimes the 
	# paired end read merging gets hampered
	echo '*** Inputs are the fasta files !!!'
	echo '===>>>> Performing the BWA alignment of the respective FASTA files'
	# alignment of the read 1
	# and removal of the chimeric reads
	if [ ! -s $BWADir/bwa_out_R1_Chimeric.sam ]; then
		bwa mem -t $THREADS $GENOME $FASTQ1 2>$BWADir/$PREFIX.R1.bwa.log | perl chimeric.pl - > $BWADir/bwa_out_R1_Chimeric.sam
	fi
	# alignment of the read 2
	# and removal of the chimeric reads
	if [ ! -s $BWADir/bwa_out_R2_Chimeric.sam ]; then
		bwa mem -t $THREADS $GENOME $FASTQ2 2>$BWADir/$PREFIX.R2.bwa.log | perl chimeric.pl - > $BWADir/bwa_out_R2_Chimeric.sam
	fi
	# now read both the reads to find out the read pairs
	# which are uniquely mapped together with a specified quality
	# also mitichondrial chromosomes, random contigs are deleted
	if [ ! -s $BWADir/bwa_out.PairFilter.sam ]; then
		python pair_up_filter.py -f $BWADir/bwa_out_R1_Chimeric.sam -r $BWADir/bwa_out_R2_Chimeric.sam -q $MAPQ_Thr 2>$BWADir/$PREFIX.PairFilterlog > $BWADir/bwa_out.PairFilter.sam
	fi

	#===================================================
	# filter the alignment based on RE
	if [ ! -s $BWADir/$PREFIX.paired.cis.RE.filtered.bam ]; then
		echo '===>>>> Applying the RE based cutting + filtering operation'
		if [[ ! -z $CUT_ENZ ]]; then
			# restriction enzyme based cut of the paired end filtered sam file
			samtools view -S -h -L $CUT_ENZ $BWADir/bwa_out.PairFilter.sam | python cutter_sites_filter.py $MIN_INSERT_SIZE - 2>>$logfile | samtools view -bhS - > $BWADir/$PREFIX.paired.cis.RE.filtered.bam
		else
			# there is no restriction enzyme provided 
			# possibly the experiment corresponds to HiChIP data
			# just convert to bam format
			samtools view -bhS $BWADir/bwa_out.PairFilter.sam > $BWADir/$PREFIX.paired.cis.RE.filtered.bam
		fi
	fi

	# sort the alignment file
	if [ ! -s $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam ]; then
		samtools sort $BWADir/$PREFIX.paired.cis.RE.filtered.bam $BWADir/$PREFIX.paired.cis.RE.filtered.sorted
	fi

else
	# condition 2 - if the input file is already aligned
	# possibly with the restriction enzyme based cutting as well
	if [[ ! -z $Aligned_Cut ]]; then

		# # debug - sourya
		# # for the input alignment file, check for any single read (ideally it should be read pairs - two consecutive lines)
		# samtools view $Aligned_Cut | awk '{a[$1]+=1;}END{for (i in a) {if (a[i] < 2) {print i, a[i];} } }' > $BWADir/check_single_read_input_alignment.txt

		#===================================================
		# extract the CIS read pairs
		
		# comment - sourya
		# samtools view -h $Aligned_Cut | sed '/chrM/d;/random/d;/chrUn/d' - | awk -v OFS='\t' '{ if(substr($1, 1, 1)=="@"){print}; if(substr($1, 1, 1)!="@" && ($7=="=" || "$3"=="$7")) {$9=$8-$4;print} }' - | samtools view -bhS - > $BWADir/$PREFIX.paired.cis.bam
		
		# add - sourya
		if [ ! -s $BWADir/$PREFIX.temp.sam ]; then
			samtools view -h $Aligned_Cut > $BWADir/$PREFIX.temp.sam
		fi
		if [ ! -s $BWADir/bwa_out.PairFilter.sam ]; then
			python pair_up_filter_BAM.py -f $BWADir/$PREFIX.temp.sam -q $MAPQ_Thr 2>$BWADir/$PREFIX.PairFilterlog > $BWADir/bwa_out.PairFilter.sam
		fi

		# # debug - sourya
		# # for the generated alignment file, check for any single read (ideally it should be read pairs - two consecutive lines)
		# samtools view $BWADir/$PREFIX.paired.cis.bam | awk '{a[$1]+=1;}END{for (i in a) {if (a[i] < 2) {print i, a[i];} } }' > $BWADir/check_single_read_cis_file.txt
		
		# if [ ! -f $BWADir/$PREFIX.paired.trans.bam ]; then
		# 	# TRANS read pairs in a bam file
		# 	# comment - sourya
		# 	# samtools view -h $Aligned_Cut | sed '/chrM/d;/random/d;/chrUn/d' - | awk '{ if(substr($1, 1, 1)=="@" || ($3!=$7 && $7!="=")) print }' | samtools view -bhS - > $BWADir/$PREFIX.paired.trans.bam
		# 	# add - sourya
		# 	samtools view -h $Aligned_Cut | awk '{ if(substr($1, 1, 1)=="@" || ($3!=$7 && $7!="=")) print }' | samtools view -bhS - > $BWADir/$PREFIX.paired.trans.bam
		# fi

		#===================================================
		# filter the alignment based on RE
		if [ ! -s $BWADir/$PREFIX.paired.cis.RE.filtered.bam ]; then
			echo '===>>>> Applying the RE based cutting + filtering operation'
			if [[ ! -z $CUT_ENZ ]]; then
				# restriction enzyme based cut of the paired end filtered sam file
				samtools view -S -h -L $CUT_ENZ $BWADir/bwa_out.PairFilter.sam | python cutter_sites_filter.py $MIN_INSERT_SIZE - 2>>$logfile | samtools view -bhS - > $BWADir/$PREFIX.paired.cis.RE.filtered.bam
			else
				# there is no restriction enzyme provided 
				# possibly the experiment corresponds to HiChIP data
				# just convert to bam format
				samtools view -bhS $BWADir/bwa_out.PairFilter.sam > $BWADir/$PREFIX.paired.cis.RE.filtered.bam
			fi
		fi

		# sort the alignment file
		if [ ! -s $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam ]; then
			samtools sort $BWADir/$PREFIX.paired.cis.RE.filtered.bam $BWADir/$PREFIX.paired.cis.RE.filtered.sorted
		fi

	else
		echo 'No input FASTQ / Paired end alignment files are provided -- exit !!!'
		return
	fi
fi


#===========================================
# remove the temporary files
#===========================================
if [ -f $BWADir/bwa_out_R1.sam ]; then
	rm $BWADir/bwa_out_R1.sam
fi
if [ -f $BWADir/bwa_out_R2.sam ]; then
	rm $BWADir/bwa_out_R2.sam
fi
if [ -f $BWADir/bwa_out_R1_Chimeric.sam ]; then
	rm $BWADir/bwa_out_R1_Chimeric.sam
fi
if [ -f $BWADir/bwa_out_R2_Chimeric.sam ]; then
	rm $BWADir/bwa_out_R2_Chimeric.sam
fi
if [ -f $BWADir/bwa_out.PairFilter.sam ]; then
	rm $BWADir/bwa_out.PairFilter.sam
fi
if [ -f $BWADir/$PREFIX.paired.cis.RE.filtered.bam ]; then
	rm $BWADir/$PREFIX.paired.cis.RE.filtered.bam
fi
if [ -f $BWADir/$PREFIX.temp.sam ]; then
	rm $BWADir/$PREFIX.temp.sam
fi

#===========================================
# now remove the duplicates of the produced sorted alignment file
# use the Picard tool
if [ ! -f $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam ]; then
	java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates INPUT=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam OUTPUT=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$BWADir/$PREFIX.picard_metrics.txt
fi

# index the generated bam file
samtools index $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
