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

#=================================
# this program accompanies the FitHiChIP code
# to create an alignment file from two input Fastq files 
# Paired end Alignment files in both BAM and PAIRIX data formats are generated
# the alignment is created using Bowtie2
# alignment is sorted and duplicates are removed
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: 

./align.sh [-h] [-f FASTQ1] [-r FASTQ2] [-g BOWTIE2_GENOME] [-d OUTDIR] [-n PREFIX] [-d OutDir] [-P picard_exec] [-t THREADS] [-a ALIGNVALIDMAX] [-l MAXFRAGLEN] [-m MAX_MEM] [-q MAPQ_THR]
Example:
./align.sh -f R1.fq.gz -r R2.fq.gz -n 'demo_prefix' -g '/home/bowtie2_index/hg19/hg19' -a 4 -m "4G" -l 1000 -d '/home/FitHiChIP_Align' -t 4 -P '/home/picard-tools/picard-tools-2.7.1/picard.jar' -q 30

Options:    

  -- required:
	-f  FASTQ1           Read one of pair-end sequencing data  [.fq|.gz|.bz2].
	-r  FASTQ2           Read two of pair-end sequencing data [.fq|.gz|.bz2].
	-g  BOWTIE2_GENOME   Bowtie2 indexed reference genome.
	-n  PREFIX           Prefix of output files.
	-d  OutDir 			 Output directory which will store all the results
	-P  picard_exec		 Full path of the Picard tool executable (used for the duplicate removal)
  -- optional:
	-t  INT              Set number of threads used for Bowtie2 [default= 1].
	-m  MAX_MEM          Set max memory of duplication removal [default = "8G"].
	-a  ALIGNVALIDMAX	 Set the number of (max) valid alignments which will be searched (default = 4)
	-l  MAXFRAGLEN 		 Set the maximum fragment length to be used for Bowtie2 alignment (default = 2000)
	-q  MAPQ_THR		 Quality value threshold, below which the mapped reads are to be removed (Default 30)

EOF
}

#============================
# reference packages / executables
#============================
# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

#============================
# a few bowtie2 related parameters
# the multimapping flag - at most, this no of valid alignments will be searched
ALIGNVALIDMAX=4

# maximum fragment length considered for valid paired end alignments
MAXFRAGLEN=2000

# the number of threads used for execution
THREADS=1

# set the default output directory
OutDir=`pwd`

# maximum memory allotted
MAX_MEM="8G"

# default prefix 
PREFIX=""

# threshold of mapq quality
MAPQ_THR=30

# default fasta formatted input files
FASTQ1=""
FASTQ2=""

while getopts "f:r:n:g:t:m:d:a:l:q:P:" opt;
do
	case "$opt" in
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		n) PREFIX=$OPTARG;;
		g) GENOME=$OPTARG;;
		d) OutDir=$OPTARG;;
		t) THREADS=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		a) ALIGNVALIDMAX=$OPTARG;;
		l) MAXFRAGLEN=$OPTARG;;
		P) picard_exec=$OPTARG;;
		q) MAPQ_THR=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $FASTQ1 ]] || [[ -z $FASTQ2 ]]; then
	echo 'User did not provide one of the fastq reads - quit !!'
	exit 1
fi

if [[ -z $GENOME ]]; then
	echo 'User did not provide Bowtie2 reference genome - quit !! '
	exit 1
fi

if [[ -z $picard_exec ]]; then
	echo 'User did not provide the absolute path of Picard executable - quit !! '
	exit 1
fi

echo '**** OutDir: '$OutDir
mkdir -p $OutDir

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

bowtie2_BAM_prefix=$OutDir$PREFIX'.align.sort.MAPQ'$MAPQ_THR
bowtie2_logfile=$OutDir$PREFIX'.align.log'

# the --mm option is used for memory mapped I/O: fast parallel execution
# output of Bowtie is a sam file
# get the uniquely mapped reads by the flag 1804 - indicates discarding any improper mapping
# source: https://github.com/kundajelab/training_camp/wiki/2.3.-Processing-the-aligned-reads
# employ the quality threshold of 20
# also remove the mitochondrial chromosome (indicated by chrM)
if [[ ! -f $bowtie2_BAM_prefix'.bam' || ! -s $bowtie2_BAM_prefix'.bam' ]]; then
	bowtie2 -k $ALIGNVALIDMAX --mm --threads $THREADS -X $MAXFRAGLEN -x $GENOME -1 $FASTQ1 -2 $FASTQ2 2>$bowtie2_logfile | sed '/chrM/d;/random/d;/chrUn/d;/chrY/d' - | samtools view -Shb -F 1804 -q $MAPQ_THR - | samtools sort - $bowtie2_BAM_prefix
fi

# now remove any PCR duplicates using Picard tool
if [[ ! -f $bowtie2_BAM_prefix'.rmdup.bam' || ! -s $bowtie2_BAM_prefix'.rmdup.bam' ]]; then
	java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates INPUT=$bowtie2_BAM_prefix'.bam' OUTPUT=$bowtie2_BAM_prefix'.rmdup.bam' ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$bowtie2_BAM_prefix'.picard_metrics.txt'
fi

# index the generated bam file
if [[ ! -f $bowtie2_BAM_prefix'.rmdup.bam.bai' || ! -s $bowtie2_BAM_prefix'.rmdup.bam.bai' ]]; then
	samtools index $bowtie2_BAM_prefix'.rmdup.bam'
fi

# now convert the duplicate removed BAM file to the target PAIRIX format
# the output file name becomes: $bowtie2_BAM_prefix'.bsorted.pairs.gz'
# and the corresponding index file becomes $bowtie2_BAM_prefix'.bsorted.pairs.gz.px2'
if [[ ! -f $bowtie2_BAM_prefix'.bsorted.pairs.gz' || ! -s $bowtie2_BAM_prefix'.bsorted.pairs.gz' ]]; then
	bam2pairs $bowtie2_BAM_prefix'.rmdup.bam' $bowtie2_BAM_prefix
fi

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------

