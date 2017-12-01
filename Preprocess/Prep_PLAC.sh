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

# This script is particularly targeted to process PLAC seq data

# Note: the program uses BWA aligner

# 	Input: fastq files or already an aligned bam file
# 	Optional input: bed file storing the restriction enzyme cuts
# 	Output: 1) sorted + duplicate removed alignment file (bam format)
# 			2) Alignment file is further divided in short (read < 1 Kb) and long (read > 10 Kb) segments
# 			3) In addition, peaks from the short segments are derived. 
# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# usage info
usage(){
cat << EOF

A) When fastq reads are provided as the input:
./Prep_PLAC.sh [-h] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-G BWA_GENOME] [-c CUT_ENZ] [-d outdir] [-p Picard_Exec] [-q MAPQ_Thr] [-i INSERT_SIZE] [-m MAX_MEM] [-t THREADS]

B) When an existing alignment is provided as the input:
./Prep_PLAC.sh [-h] [-X Aligned_Read] [-n PREFIX] [-c CUT_ENZ] [-d outdir] [-p Picard_Exec] [-q MAPQ_Thr] [-i INSERT_SIZE] [-m MAX_MEM] [-t THREADS]

Example:

A) Given input fasta files and restriction enzyme file, and a quality threshold of 20:
./Prep_PLAC.sh -f R1.fq.gz -r R2.fq.gz -n 'demo_prefix' -G 'bwa_ref_genome.fa' -c RE_cutter.sites.bed -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar -q 20

B) Given input fasta files but no restriction enzyme file (and a default quality threshold):
./Prep_PLAC.sh -f R1.fq.gz -r R2.fq.gz -n 'demo_prefix' -G 'bwa_ref_genome.fa' -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar

C) Given input alignment file (BAM format) - here user needs to apply the restriction enzyme (if required) before providing the bam file as the input
./Prep_PLAC.sh -X inp_align.bwa -n 'demo_prefix' -m "8G" -t 4 -i 10000 -d /home/sourya/ -p /home/picard-tools/picard-tools-2.7.1/picard.jar

Output:

1) A sorted and filtered alignment file with duplicate segments removed
2) Separate alignment files of short reads (read < 1 Kb) and long reads (read > 10 Kb)
3) Peaks derived from the short reads

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
	-g  GSIZE		 	 Genome size parameter for MACS2 peak calling (default = 'hs')

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
# genome size parameter for MACS2 peak calling
GSIZE='hs'

# # default executable of the picard tool
# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

while getopts "f:r:n:G:c:i:d:t:m:p:X:q:g:" opt;
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
		g) GSIZE=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ $MIN_INSERT_SIZE -lt 1000 ]]; then 
	echo "error: -i '$MIN_INSERT_SIZE' could not smaller than 1000."; 
	exit 1;
fi

if [[ -z $FASTQ1 || -z $FASTQ2 ]]; then
	if [[ -z $Aligned_Cut ]]; then
		echo 'User neither provided an input alignment file, nor provided two fastq reads - quit !!'
		exit 1
	fi
fi

if [[ ! -z $FASTQ1 && ! -z $FASTQ2 ]]; then
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
		samtools sort -o $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam $BWADir/$PREFIX.paired.cis.RE.filtered.bam 
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
			samtools sort -o $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam $BWADir/$PREFIX.paired.cis.RE.filtered.bam 
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

#==============================================
# convert the duplicate removed BAM file into pairix format
# the output file name becomes: $BWADir/$PREFIX'.paired.cis.RE.bsorted.pairs.gz'
# and the corresponding index file becomes $BWADir/$PREFIX'.paired.cis.RE.bsorted.pairs.gz.px2'
if [[ ! -f $BWADir/$PREFIX'.paired.cis.RE.bsorted.pairs.gz' || ! -s $BWADir/$PREFIX'.paired.cis.RE.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $BWADir/$PREFIX'.paired.cis.RE.filtered.sorted.nodup.bam' $BWADir/$PREFIX'.paired.cis.RE'
fi

#===========================================
# create the short and long range segments from the generated alignment file
# **** applicable for the PLAC seq pipeline ****
#===========================================
SegDir=$OutDir'/Segments'
mkdir -p $SegDir
echo "Directory storing the short and long range segments: "$SegDir

echo "splitting cis reads into long-range ($MIN_INSERT_SIZE) and short-range (<1kb) ..." 

if [ ! -f $SegDir/$PREFIX.long.bam || ! -f $SegDir/$PREFIX.short.bam ]; then
	# long range mapped read
	samtools view -h $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam | awk -v var="$MIN_INSERT_SIZE" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) >= var) print }' | samtools view -bhS - > $SegDir/$PREFIX.long.bam
	# short range mapped read
	samtools view -h $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam | awk 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) < 1000) print }' | samtools view -bhS - > $SegDir/$PREFIX.short.bam
	# logs
	UNIQ_PAIRED_CIS_FILTERED_SORTED_LONG_NODUP=$(samtools view $SegDir/$PREFIX.long.bam | wc -l | awk '{print $1}')
	UNIQ_PAIRED_CIS_FILTERED_SORTED_SHORT_NODUP=$(samtools view $SegDir/$PREFIX.short.bam | wc -l | awk '{print $1}')

	echo "number of monoclonal cis long-range pairs $((UNIQ_PAIRED_CIS_FILTERED_SORTED_LONG_NODUP/2))" | tee -a $logfile
	echo "number of monoclonal cis short-range pairs $((UNIQ_PAIRED_CIS_FILTERED_SORTED_SHORT_NODUP/2))" | tee -a $logfile

	# index creation of the resulting bam files
	# this is a crucial step - Sourya
	# since in the interaction processing phase (process_int.sh)
	# the alignment BAM files need to be random accessed via the pysam functions fetch() and pileup()
	# which is possible if the corresponding index file is present
	samtools index $SegDir/$PREFIX.long.bam
	samtools index $SegDir/$PREFIX.short.bam
fi

#===============================================
# convert the long and short BAM files into pairix format
if [[ ! -f $SegDir/$PREFIX'.long.bsorted.pairs.gz' || ! -s $SegDir/$PREFIX'.long.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $SegDir/$PREFIX.long.bam $SegDir/$PREFIX'.long'
fi
if [[ ! -f $SegDir/$PREFIX'.short.bsorted.pairs.gz' || ! -s $SegDir/$PREFIX'.short.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $SegDir/$PREFIX.short.bam $SegDir/$PREFIX'.short'
fi
#===========================================
# now obtain the peak segments from the generated short read file
# this will be useful for mapping the peaks from the short segments to the long range reads
#===========================================
macs2dir=$OutDir'/PeaksAnchors_ShortSegment'
mkdir -p $macs2dir

# apply the short range alignment file on MACS2
# Note: we use the default MACS2 setting for duplicate peaks 
macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $SegDir/$PREFIX.short.bam

grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed

#===========================================
# here we get the peak segments from the input complete alignment file
#===========================================
macs2dir=$OutDir'/PeaksAnchors_ALL'
mkdir -p $macs2dir

# apply the complete input alignment file on MACS2
# Note: we use the default MACS2 setting for duplicate peaks 
macs2 callpeak -f AUTO -g $GSIZE --keep-dup all --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam

grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
