#!/bin/bash -ex
#PBS -l nodes=1:ppn=4
#PBS -l mem=10GB
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

#===============
# Preprocessing Code of the HiChIP alignment data
# to accompany FitHiChIP framework

# Input fastq files (two single end reads) are processed to generate the following segments:
# 1) An alignment file generated from the fastq reads
# 2) Reads of short segments (< 1 Kb default) with the filtered read pair orientation
# 3) Reads of long segments (> 10 Kb default) with the filtered read pair orientation
# 4) Peaks detected from the short read segments
# 5) Peaks detected from the complete input alignment (if corresponding option is specified)
# 
# Parameters:
# 1) Input fastq read 1
# 2) Input fastq read 2
# 3) Output directory (optional)
# 4) Boolean option to get the peaks of the complete input alignment
# 5) Picard executable (for duplicate removal)
# 6) Distance threshold of the short reads (default 1000: means < 1 Kb)
# 7) Distance threshold of the long reads (default 10000: means > 10 Kb)

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# usage info
usage(){
cat << EOF

Options:    

	-f  FASTQ1           R1 of pair-end sequencing data  [.fq|.gz|.bz2].
	-r  FASTQ2           R2 of pair-end sequencing data [.fq|.gz|.bz2].
	-G  BWA_GENOME       BWA indexed reference genome (if the fasta reads are provided as inputs) (mandatory parameter)
	-c  CUT_ENZ          File with the restriction cut site (bed formatted)
	-d  OutDir 			 Output base directory [Default: current working directory]
	-p  picard_exec      Location of picard executable (.jar file) (mandatory parameter)
	-C  CompletePeak	 Boolean option - if 1, prints the peaks for the complete input alignment file also. Default 0.
  	-q  MAPQ_Thr		 Threshold of mapping quality [Default = 30]
	-m  MAX_MEM          Set max memory of duplication removal [Default = 1G].
	-g  GSIZE		 	 Genome size parameter for MACS2 peak calling (default = 'hs')
	-n  PREFIX           Prefix of output file names [Default: empty string].
	-S 	ShortReadDistThr Threshold of Distance defining the short reads. Default 1000 (< 1 Kb)
	-L  LongReadDistThr  Threshold of Distance defining the long reads. Default 10000 (> 10 Kb)
	-t  INT              Set number of sorting, BWA mapping threads [Default = 1].
EOF
}

# default parameters
OutDir=`pwd`	#'./'
# the number of threads used for execution
THREADS=1
# maximum memory to be used
MAX_MEM="1G"
# default input fastq files
FASTQ1=""
FASTQ2=""
# default value of the restriction enzyme cutter file
CUT_ENZ=""
# this is the default quality threshold employed for each alignment
MAPQ_Thr=30
# default prefix 
PREFIX=""
# genome size parameter for MACS2 peak calling
GSIZE='hs'
# printing peaks of the complete "InpFile"
CompletePeak=0

# default values of short and long read thresholds
ShortReadDistThr=1000
LongReadDistThr=10000

# # default executable of the picard tool
# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

while getopts "f:r:G:c:n:d:m:p:q:g:C:S:L:t:" opt;
do
	case "$opt" in
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		G) GENOME=$OPTARG;;
		c) CUT_ENZ=$OPTARG;;
		n) PREFIX=$OPTARG;;
		d) OutDir=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		p) picard_exec=$OPTARG;;
		q) MAPQ_Thr=$OPTARG;;
		g) GSIZE=$OPTARG;;
   		C) CompletePeak=$OPTARG;;
		S) ShortReadDistThr=$OPTARG;;
		L) LongReadDistThr=$OPTARG;;
		t) THREADS=$OPTARG;;
		?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

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

SegDir=$OutDir'/Segments_HiChIP'
mkdir -p $SegDir
echo "Directory storing the short and long range segments: "$SegDir

#==========================
# filenames used in this script
#==========================

BWA_R1_Alignfile=$BWADir/bwa_out_R1_Chimeric.sam
BWA_R2_Alignfile=$BWADir/bwa_out_R2_Chimeric.sam
PairedFilteredReadFile=$BWADir/$PREFIX.paired.cis.RE.filtered.bam
SortPrefix=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted
rmDupFile=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam

ShortReadFile=$SegDir'/'$PREFIX'.cis.short.bam'
LongReadFile=$SegDir'/'$PREFIX'.cis.long.bam'
DumpedPairsSAMFile=$SegDir'/'$PREFIX'.cis.DumpedPairs.sam'
DumpedPairsBAMFile=$SegDir'/'$PREFIX'.cis.DumpedPairs.bam'

#==========================
# align the fasta files, to form the paired end read alignment
echo '===>>>> Performing the BWA alignment of the respective FASTA files'

# alignment of the read 1 and removal of the chimeric reads
if [ ! -s $BWA_R1_Alignfile ]; then
	bwa mem -t $THREADS $GENOME $FASTQ1 2>$BWADir/$PREFIX.R1.bwa.log | perl chimeric.pl - > $BWA_R1_Alignfile
fi

# alignment of the read 2 and removal of the chimeric reads
if [ ! -s $BWA_R2_Alignfile ]; then
	bwa mem -t $THREADS $GENOME $FASTQ2 2>$BWADir/$PREFIX.R2.bwa.log | perl chimeric.pl - > $BWA_R2_Alignfile
fi

# merge two single end read files
# subject to HiChIP specific constraints (read orientation, quality mapping)
# also perform RE based cutting and filtering if available
echo '===>>>> Merging the single end reads'
#if [ ! -s $PairedFilteredReadFile ]; then
	if [[ ! -z $CUT_ENZ ]]; then
		echo '===>>>> Applying the RE based cutting + filtering operation'
		python pair_up_filter_HiChIP.py -f $BWA_R1_Alignfile -r $BWA_R2_Alignfile -q $MAPQ_Thr -W $DumpedPairsSAMFile 2>$BWADir/$PREFIX.PairFilterlog | samtools view -S -h -L $CUT_ENZ - | samtools view -bhS - > $PairedFilteredReadFile
	else
		python pair_up_filter_HiChIP.py -f $BWA_R1_Alignfile -r $BWA_R2_Alignfile -q $MAPQ_Thr -W $DumpedPairsSAMFile 2>$BWADir/$PREFIX.PairFilterlog | samtools view -bhS - > $PairedFilteredReadFile
	fi
#fi

# sort the alignment file
#if [ ! -s $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.bam ]; then
	samtools sort -o $SortPrefix'.bam' $PairedFilteredReadFile 
#fi

#if [ ! -f $BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam ]; then
	java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates INPUT=$SortPrefix'.bam' OUTPUT=$rmDupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$BWADir/$PREFIX.picard_metrics.txt

	# index the generated bam file
	samtools index $rmDupFile
#fi

#==============================================
# convert the duplicate removed BAM file into pairix format
# the output file name becomes: $BWADir'/'$PREFIX'.paired.cis.RE.bsorted.pairs.gz'
# and the corresponding index file becomes $BWADir'/'$PREFIX'.paired.cis.RE.bsorted.pairs.gz.px2'
if [[ ! -f $BWADir'/'$PREFIX'.paired.cis.RE.bsorted.pairs.gz' || ! -s $BWADir'/'$PREFIX'.paired.cis.RE.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $rmDupFile $BWADir'/'$PREFIX'.paired.cis.RE'
fi

#===========================================
# create the short and long range segments from the generated alignment file
# applicable for the HiChIP pipeline
#===========================================

# now create the short 'cis' reads (with the filtered read orientation) (< ShortReadDistThr)
# and place it in an output alignment file
#if [ ! -f $ShortReadFile ]; then
	samtools view -h $rmDupFile | awk -v dt="$ShortReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) < dt) print }' | samtools view -bhS - > $ShortReadFile
	samtools index $ShortReadFile
#fi

# now get the long 'cis' reads (with the filtered read orientation) (> LongReadDistThr)
# and place it in an output alignment file
#if [ ! -f $LongReadFile ]; then
	samtools view -h $rmDupFile | awk -v var="$LongReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) >= var) print }' | samtools view -bhS - > $LongReadFile
	samtools index $LongReadFile
#fi

#===============================================
# convert the long and short BAM files into pairix format
if [[ ! -f $SegDir'/'$PREFIX'.cis.long.bsorted.pairs.gz' || ! -s $SegDir'/'$PREFIX'.cis.long.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $LongReadFile $SegDir'/'$PREFIX'.cis.long'
fi
if [[ ! -f $SegDir'/'$PREFIX'.cis.short.bsorted.pairs.gz' || ! -s $SegDir'/'$PREFIX'.cis.short.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $ShortReadFile $SegDir'/'$PREFIX'.cis.short'
fi

# convert the dumped reads to bam format and 
# remove the temporary sam file
#if [ -f $DumpedPairsSAMFile ]; then
	samtools view -bhS $DumpedPairsSAMFile > $DumpedPairsBAMFile
	rm $DumpedPairsSAMFile
#fi

#===========================================
# now obtain the peak segments from the generated short read file
# this will be useful for mapping the peaks from the short segments to the long range reads
#===========================================
macs2dir=$OutDir'/PeaksAnchors_ShortSegment_Fastq_Latest'
mkdir -p $macs2dir

# apply the short range alignment file on MACS2
# Note: we use the default MACS2 setting for duplicate peaks
# we also use a FDR cutoff of 0.01, as suggested in the HiChIP paper (Mumbach et. al. 2016) 
macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $ShortReadFile -q 0.01
grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed

#===========================================
# here we get the peak segments from the input complete alignment file
#===========================================
if [ $CompletePeak == 1 ]; then
	macs2dir=$OutDir'/PeaksAnchors_CompleteAlignment_Fastq_Latest'
	mkdir -p $macs2dir

	# apply the complete input alignment file on MACS2
	# Note: we use the default MACS2 setting for duplicate peaks 
	# we also use a FDR cutoff of 0.01, as suggested in the HiChIP paper (Mumbach et. al. 2016) 
	macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $rmDupFile -q 0.01
	grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed
fi

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------



