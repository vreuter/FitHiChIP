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

# Input alignment (BAM) file is processed to generate the following segments:
# 1) Reads of short segments (< 1 Kb default) with the filtered read pair orientation
# 2) Reads of long segments (> 10 Kb default) with the filtered read pair orientation
# 3) Peaks detected from the short read segments
# 4) Peaks detected from the complete input alignment (if corresponding option is specified)
# 
# Parameters:
# 1) Input alignment file
# 2) Output directory (optional)
# 3) Boolean option to get the peaks of the complete input alignment
# 4) Picard executable (for duplicate removal)
# 5) Boolean option specifying whether the input alignment file has two lines for a paired end read (1), or a single line (0)
#    - default 1
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

	-I  InpFile          Input alignment file (mandatory parameter)
	-d  OutDir 			 Output base directory [Default: current working directory]
	-p  picard_exec      Location of picard executable (.jar file) (mandatory parameter)
	-C  CompletePeak	 Boolean option - if 1, prints the peaks for the complete input alignment file also. Default 0.
	-R  RepeatLine 		 Boolean option - if 1, signifies that two lines in the bam file denote one paired end read. 
						 Else, different lines of the bam file denote different reads. Default 1.
  	-q  MAPQ_Thr		 Threshold of mapping quality [Default = 30]
	-m  MAX_MEM          Set max memory of duplication removal [Default = 1G].
	-g  GSIZE		 	 Genome size parameter for MACS2 peak calling (default = 'hs')
	-n  PREFIX           Prefix of output file names [Default: empty string].
	-S 	ShortReadDistThr Threshold of Distance defining the short reads. Default 1000 (< 1 Kb)
	-L  LongReadDistThr  Threshold of Distance defining the long reads. Default 10000 (> 10 Kb)

EOF
}

# default parameters
OutDir=`pwd`'/'
# maximum memory to be used
MAX_MEM="1G"
# default input file
InpFile=""
# this is the default quality threshold employed for each alignment
MAPQ_Thr=30
# default prefix 
PREFIX=""
# genome size parameter for MACS2 peak calling
GSIZE='hs'
# printing peaks of the complete "InpFile"
CompletePeak=0
# assuming that the paired end reads are spanned in two lines
RepeatLine=1

# default values of short and long read thresholds
ShortReadDistThr=1000
LongReadDistThr=10000

# threshold of FDR employed for MACS2
# in the HiChIP pipeline (Mumbach et. al. 2016), 0.01 is employed
Q_val_thr=1e-2

# # default executable of the picard tool
# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

while getopts "I:n:d:m:p:q:g:C:R:S:L:" opt;
do
	case "$opt" in
		I) InpFile=$OPTARG;;
		n) PREFIX=$OPTARG;;
		d) OutDir=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		p) picard_exec=$OPTARG;;
		q) MAPQ_Thr=$OPTARG;;
		g) GSIZE=$OPTARG;;
   		C) CompletePeak=$OPTARG;;
		R) RepeatLine=$OPTARG;;
		S) ShortReadDistThr=$OPTARG;;
		L) LongReadDistThr=$OPTARG;;
		?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $InpFile ]]; then
	echo 'Input BAM alignment file is not provided - quit !!'
	exit 1
fi

if [[ -z $picard_exec ]]; then
	echo 'PATH of picard executable tool is not provided - quit !!'
	exit 1
fi

# output directory for storing the results
mkdir -p $OutDir
echo '**** Output directory for storing the results: '$OutDir

SegDir=$OutDir'/Segments_HiChIP'
mkdir -p $SegDir
echo "Directory storing the short and long range segments: "$SegDir

#=======================================
# first create a sorted and duplicate removed BAM file
# with 1) cis reads 2) opposite strands 3) valid chromosomes
# 4) quality thtesholds
# if the flag AND 16 (bitwise AND) is 1, then the strand is mapped for the reverse strand
# else it is mapped in the forward strand
# also check if the length of the chromosome name is <= 5 (to delete the other chromosome contigs)
# *** also calculate the insert size (difference between two read positions)
#=======================================
rmDupFile=$SegDir'/'$PREFIX'.cis.sort.rmdup.bam'
tempSortFilePrefix=$SegDir'/'$PREFIX'.cis.Align.temp'

if [ $RepeatLine == 1 ]; then
	samtools view -h $InpFile | awk -v qt="$MAPQ_Thr" '{if(substr($1, 1, 1)!="@") {l1=$0; chr1=$3; flag1=$2; q1=$5; strand1=(and(16, flag1)); getline; l2=$0; chr2=$3; flag2=$2; q2=$5; strand2=(and(16, flag2)); {if (chr1==chr2 && length(chr1)<=5 && q1>=qt && q2>=qt && strand1 != strand2) {printf "%s\n%s\n", l1, l2}}} else {print $0}}' | awk '{if(substr($1, 1, 1)!="@") {$9=($8-$4); print $0} else {print $0}}' | perl -ne 'chomp $_; $_=~s/\s/\t/g; if($_!~/\@PG/){print "$_\n";}' | samtools view -bhS - > $tempSortFilePrefix'.bam' 
fi

# now sort the generated bam file
samtools sort -o $tempSortFilePrefix'.sort.bam' $tempSortFilePrefix'.bam' 

# remove duplicates
java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates INPUT=$tempSortFilePrefix'.sort.bam' OUTPUT=$rmDupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$SegDir/$PREFIX.cis.picard_metrics.txt

samtools index $rmDupFile

#==============================================
# convert the duplicate removed BAM file into pairix format
# the output file name becomes: $SegDir'/'$PREFIX'.cis.bsorted.pairs.gz'
# and the corresponding index file becomes $SegDir'/'$PREFIX'.cis.bsorted.pairs.gz.px2'
if [[ ! -f $SegDir'/'$PREFIX'.cis.bsorted.pairs.gz' || ! -s $SegDir'/'$PREFIX'.cis.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $rmDupFile $SegDir'/'$PREFIX'.cis'
fi

#==================================
# get the reads with length < ShortReadDistThr
# check read Distance
#==================================
ShortReadFile=$SegDir'/'$PREFIX'.cis.short.bam'
samtools view -h $rmDupFile | awk -v dt="$ShortReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) < dt) print }' | samtools view -bhS - > $ShortReadFile
samtools index $ShortReadFile

#==================================
# get the reads with length > LongReadDistThr
# check read Distance
#==================================
LongReadFile=$SegDir'/'$PREFIX'.cis.long.bam'
samtools view -h $rmDupFile | awk -v var="$LongReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) >= var) print }' | samtools view -bhS - > $LongReadFile
samtools index $LongReadFile

#===============================================
# convert the long and short BAM files into pairix format
if [[ ! -f $SegDir'/'$PREFIX'.cis.long.bsorted.pairs.gz' || ! -s $SegDir'/'$PREFIX'.cis.long.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $LongReadFile $SegDir'/'$PREFIX'.cis.long'
fi
if [[ ! -f $SegDir'/'$PREFIX'.cis.short.bsorted.pairs.gz' || ! -s $SegDir'/'$PREFIX'.cis.short.bsorted.pairs.gz.px2' ]]; then
	bam2pairs $ShortReadFile $SegDir'/'$PREFIX'.cis.short'
fi

#=======================================
# 'cis' reads which are dumped pairs (same strand orientation)
#=======================================
samtools view -h $InpFile | awk '{if(substr($1, 1, 1)!="@") {l1=$0; chr1=$3; flag1=$2; strand1=(and(16, flag1)); getline; l2=$0; chr2=$3; flag2=$2; strand2=(and(16, flag2)); {if(chr1==chr2 && strand1 == strand2) {printf "%s\n%s\n", l1, l2}}} else {print $0}}' - | samtools view -bhS - > $SegDir'/'$PREFIX'.cis.DumpedPairs.bam'

#===========================================
# now obtain the peak segments from the generated short read file
# this will be useful for mapping the peaks from the short segments to the long range reads
#===========================================
macs2dir=$OutDir'/PeaksAnchors_ShortSegment_BAM_Latest'
mkdir -p $macs2dir

# apply the short range alignment file on MACS2
# Note: we use the default MACS2 setting for duplicate peaks
# we also use a FDR cutoff of 0.01, as suggested in the HiChIP paper (Mumbach et. al. 2016) 
macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $ShortReadFile -q $Q_val_thr
grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed

#===========================================
# here we get the peak segments from the input complete alignment file
#===========================================
if [ $CompletePeak == 1 ]; then
	macs2dir=$OutDir'/PeaksAnchors_CompleteAlignment_BAM_Latest'
	mkdir -p $macs2dir

	# apply the complete input alignment file on MACS2
	# Note: we use the default MACS2 setting for duplicate peaks 
	# we also use a FDR cutoff of 0.01, as suggested in the HiChIP paper (Mumbach et. al. 2016) 
	macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $rmDupFile -q $Q_val_thr
	grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed

fi



