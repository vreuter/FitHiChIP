#!/bin/bash

#===============
# main executable of FitHiChIP
# When the input and output data format can be in 
# 1) BAM format
# 2) PAIRIX format (https://github.com/4dn-dcic/pairix)

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# usage info
usage(){
cat << EOF

usage: ./FitHiChIP.sh [-h] [-I InpFile] [-M METHOD] [-P PEAKFILE] [-o OUTDIR] [-t THREADS] [-n PREFIX] [-b BIN_SIZE] [-q QVALUE] [-v verbose (1 / 0)] [-L LowDistThres] [-U UppDistThres] [-B 0/1] [-D 0/1] [-f 0/1] [-N NBINS_FitHiC]

*************************
Examples using BAM file as input:
*************************

1) FitHiChIP.sh -I inp.bam -M 3 -n 'FitHiChIP' -v 1 -t 8 -L 20000 -U 2000000 -f 1 -D 1 -B 1 -N 200 -b 5000 -t 8 -o '/home/sourya/FitHiChIP/'
Interactions among all possible segments (fixed size bins)
2) FitHiChIP.sh -I inp.bam -M 1 -P 'inpPeak.narrowPeak' -n 'FitHiChIP' -v 1 -t 8 -L 20000 -U 2000000 -f 1 -D 1 -B 1 -N 200 -b 5000 -t 8 -o '/home/sourya/FitHiChIP/'
Interactions between peaks (obtained from the peak input file)
3) FitHiChIP.sh -I inp.bam -M 2 -n 'FitHiChIP' -v 1 -t 8 -L 20000 -U 2000000 -f 1 -D 1 -B 1 -N 200 -b 5000 -t 8 -o '/home/sourya/FitHiChIP/'
Interactions from peak segments to the non-peak segments (Here the peak information is not separately 
provided. But the alignment file is in BAM format. So, peaks are computed from this alignment file, using the package MACS2.
These peaks are then used for finding the interactions.
4) FitHiChIP.sh -I inp.bam -M 0 -n 'FitHiChIP' -v 1 -t 8 -L 20000 -U 2000000 -f 1 -D 1 -B 1 -N 200 -b 5000 -t 8 -o '/home/sourya/FitHiChIP/'
Interactions from the peak segments to ALL fixed size bins (peaks or non-peaks) are computed.

*************************
Examples using PAIRIX file as input:
*************************
FitHiChIP.sh -I inp.pairix -M 2 -P 'inpPeak.narrowPeak' -n 'FitHiChIP' -v 1 -t 8 -L 20000 -U 2000000 -f 1 -D 1 -B 1 -N 200 -b 5000 -t 8 -o '/home/sourya/FitHiChIP/'


Options:

  -- required:
   	-I  Inp			Input alignment file (complete alignment). The file can be either in BAM format, 
	       			or in PAIRIX format. Conversion between BAM to PAIRIX can be done by 
	       			using the utility bam2pairs (https://github.com/4dn-dcic/pairix)
	-M  Met 		Can be 1, 2, or 3 (default 3), where:
					3 (ALL): Finds the statistical significant interactions among all the binned segments.
					Here user does not need to supply the peak detection output file.
					2 (Peak to Non Peak): Finds the statistical significant interactions from the peaks to the non peak segments 
					(subject to the binning employed).
					In the interaction file, Left side is a peak segment and right side is a non peak segment. 
					Possible only if either the user supplies the peak file via -P option, or 
					supplies the complete alignment file in BAM format.
					1 (Peak to Peak): Finds the statistical significant interactions between the peaks 
					(subject to the binning employed).
					In the interaction file, Both left and right sides are peak segments.  
					Possible only if either the user supplies the peak file via -P option, or 
					supplies the complete alignment file in BAM format.
					0 (Peak to ALL): Finds the statistical significant interactions from the peaks 
					(subject to the binning employed) to ALL fixed size binned segments.
					In the interaction file, Left side is a peak segment. 
					Right side can be either a peak or a non peak segment.
					Possible only if either the user supplies the peak file via -P option, or 
					supplies the complete alignment file in BAM format.
   	-P  Peak 		Input file containing the peak regions (obtained by the packages like MACS2, etc.). 
	       			The file should be in .bed format.
	       			When the method (-M option) is either 0, 1 or 2, this file is used.
	       			If the user does not supply the peak file but provides the input alignment file 
	       			in BAM format, peaks are computed using MACS2.
	       			For alignment file in PAIRIX format, user must supply this peak file separately.
	-o  OutDir 		Output directory storing the output results
	-n  PREFIX      	Prefix string of output files (Default = 'FitHiChIP').
	-b  BINSIZE     	Size of the bins [default = 5000], in bases, for interaction detection.

  -- optional:
    	-t  INT        		Set number of threads for peak calling and finding the interactions [default = 1].
    	-g  GSIZE		Genome size parameter for MACS2 (if peak detection is sought). Default = 'hs'
    	-q  QVALUE      	Minimum FDR (q-value) cutoff for interaction detection [default = 1e-2].
    	-v  verbose		Specified as 1 or 0; if 1, generate a timing profile information
   	-L  LowDistThr 		Lower distance threshold of interaction between two segments (default = 20000)
    	-U  UppDistThr		Upper distance threshold of interaction between two segments (default = 200000)
	-f  FitHiCBin		If 1 (default), equal occupancy (contact count) bins will be employed in the FitHiC code. 
	 				Else, equal length bins will be used.
	-N  NBins		In the FitHiC model, this is the max no of bins (equal occupancy or equal length) 
	 				that is to be employed. Default 200.
	-B  Binom		If 1, computes the binomial distribution model between the observed 
					genomic distance and the contact count, 
					with respect to the paper Duan et. al. (2010). Default is 0.
	-D  Draw		Specified as 1 or 0. If 1, draws the figures of various statistics / analysis. Default 0.

EOF
}

# default params
THREADS=1

# default prefix
PREFIX='FitHiChIP'

# this corresponds to an interval of 5 Kb resolution (as mentioned in the main paper)
# if the interaction between short segment peaks and the long segments are to be detected
# then this bin size is used to get the anchor region corresponding to every peak
BIN_SIZE=5000

# genome size parameter for MACS2
GSIZE='hs'

# sample p and q values
# PVALUE=1e-5
QVALUE=1e-2

# default output directory for storing the results
OutDir=`pwd`

# default option of Verbose (output of the timing profile)
Verbose=0

# default value of lower distance threhold between two interacting segments
LowDistThres=20000

# default value of upper distance threhold between two interacting segments
UppDistThres=2000000

# binning method for FitHiC technique
# 0 - equal length bin, 1 for equal occupancy bin (default)
FitHiCBinMethod=1

# The file with peak detection output
PeakFile=""

# if 1, computes the binomial distribution model of genomic distance vs contact count, 
# with respect to the paper Duan et. al.
BinomDistrModel=0

# default value of the number of bins that is to be employed for FitHiC
NBins=200

# default value of method (type of interaction detection)
Method=3

# default value of plotting analysis figures
DrawFig=0

while getopts "I:P:t:n:b:q:o:v:L:U:f:B:N:M:D:g:" opt;
do
	case "$opt" in
		I) INPFILE=$OPTARG;;
		P) PeakFILE=$OPTARG;;
		n) PREFIX=$OPTARG;;
		t) THREADS=$OPTARG;;
		b) BIN_SIZE=$OPTARG;;
		q) QVALUE=$OPTARG;;
		o) OutDir=$OPTARG;;
		v) Verbose=$OPTARG;;
		L) LowDistThres=$OPTARG;;
		U) UppDistThres=$OPTARG;;
		f) FitHiCBinMethod=$OPTARG;;
		B) BinomDistrModel=$OPTARG;;
		N) NBins=$OPTARG;;
		M) Method=$OPTARG;;
		D) DrawFig=$OPTARG;;
		g) GSIZE=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $INPFILE ]]; then
	echo 'No input alignment file is provided - exit !!'
	exit 1
fi

# input file format - 1 means bam file, 2 means pairix file
extension="${INPFILE##*.}"
if [[ "$extension" == "bam" ]]; then
	inpfilefmt=1
else
	inpfilefmt=2
fi

# check if the input peak file is not provided for the pairix file format
if [[ $Method != 3 ]]; then
	if [[ -z $PeakFILE ]]; then
		if [[ $inpfilefmt == 2 ]]; then
			echo 'User requires interaction from peaks but did not provide the peak detection output file (via -P option) and neither has provided a BAM alignment file !!' 
			exit 1
		fi
	fi
fi

#==============================
# convert the input files to their absolute path
# input file
INPFILE="$(cd "$(dirname "$INPFILE")"; pwd)/$(basename "$INPFILE")"
# peak file if provided
if [[ ! -z $PeakFILE ]]; then
	PeakFILE="$(cd "$(dirname "$PeakFILE")"; pwd)/$(basename "$PeakFILE")"
fi
# output directory
if [ ! -d $OutDir ]; then
    mkdir -p $OutDir
fi
OutDir=`cd "$OutDir"; pwd`
#==============================

#==============================
# important - sourya
# first change the current working directory to the directory containing this script
# it is useful when the script is invoked from a separate directory
#==============================
currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir

#==============================
# if the target method is 0, 1, or 2
# if user does not provide a peak detection output but provides the bam alignment file
# we compute the peaks via MACS2
#==============================
if [[ $Method != 3 ]]; then
	macs2dir=$OutDir'/PeaksAnchors'
	mkdir -p $macs2dir
	if [[ ! -z $PeakFILE ]]; then
		if [ ! -s $macs2dir/$PREFIX.anchors.bed ]; then

			# echo 'Extend the input peaks in either direction to span the BIN SIZE'
			
			# option 1 - flank shorter peaks to 5 Kb but keep longer peaks (> 5 Kb) intact
			# grep -vwE "chrM" $PeakFILE | sort -k5,5rn - | awk -v EXTEND_LEN="$((BIN_SIZE / 2))" 'function max(x,y) {return x>y?x:y}; function min(x,y) {return x>y?y:x}; {printf "%s\t%d\t%d\t%s\t%s\n", $1, max(1, min($2,(($2+$3)/2)-EXTEND_LEN)), max($3,(($2+$3)/2)+EXTEND_LEN), $4, $5 }' - > $macs2dir/$PREFIX.anchors.bed

			# option 2 - flank shorter peaks to 5 Kb and trim longer peaks to 5 Kb as well
			# thus enforce fixed sized bins (= bin size specified)
			# grep -vwE "chrM" $PeakFILE | sort -k5,5rn - | awk -v EXTEND_LEN="$((BIN_SIZE / 2))" 'function max(x,y) {return x>y?x:y}; {printf "%s\t%d\t%d\t%s\t%s\n", $1, max(1, ($2+$3)/2-EXTEND_LEN), ($2+$3)/2+EXTEND_LEN, $4, $5 }' - > $macs2dir/$PREFIX.anchors.bed	

			# option 3 - just copy the peaks as it is
			# after removing mitochondrial and random genome
			# print the first five fields
			grep -vwE "chrM" $PeakFILE | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed
		fi
	else
		# input alignment BAM file is applied on MACS2
		if [ ! -s $macs2dir/$PREFIX'_peaks.narrowPeak' ]; then
			# Note: we do not filter any peaks but rather call all the peaks
			macs2 callpeak -f AUTO -g $GSIZE --keep-dup all --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $INPFILE
		fi
		if [ ! -s $macs2dir/$PREFIX'.anchors.bed' ]; then
			# option 1
			# extend the peaks on either side by "EXTEND_LEN" and store in in a separate anchor file 
			# filter the peaks to remove chrM information - important
			# the size of anchor region is the same as the bin size
			# Note: during extension of the peaks in either side, check that the values should be positive
			# awk -v EXTEND_LEN="$((BIN_SIZE / 2))" 'function max(x,y) {return x>y?x:y}; {printf "%s\t%d\t%d\t%s\t%s\n", $1, max(($2+$3)/2-EXTEND_LEN, 1), ($2+$3)/2+EXTEND_LEN, $4, $5 }' <(grep -vwE "chrM" $macs2dir/$PREFIX\_summits.bed | sort -k5,5rn -) > $macs2dir/$PREFIX.anchors.bed

			# option 2 - copy the peaks as it is
			grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed
		fi
	fi

	# check if peak count > 0
	# otherwise exit
	npeak=`cat $macs2dir/$PREFIX.anchors.bed | wc -l`
	if [[ $npeak -eq 0 ]]; then
		echo 'Number of peaks = 0 - quit !!'
		exit 1
	fi 
fi

#==========================
# directory storing all the outputs
# 1) GenFitHiCBaseDir: stores the basic and sorted interactions
# 2) Separate directory underlying $GenFitHiCBaseDir storing the spline fit or binomial distribution
#==========================
if [ $Method == 3 ]; then
	GenFitHiCBaseDir=$OutDir'/FitHiChIP_ALL2ALL'
elif [ $Method == 2 ]; then
	GenFitHiCBaseDir=$OutDir'/FitHiChIP_Peak2NonPeak'
elif [ $Method == 1 ]; then
	GenFitHiCBaseDir=$OutDir'/FitHiChIP_Peak2Peak'
else
	GenFitHiCBaseDir=$OutDir'/FitHiChIP_Peak2ALL'
fi
GenFitHiCBaseDir=$GenFitHiCBaseDir'_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres

echo "Base directory of generalized FitHiC: "$GenFitHiCBaseDir
mkdir -p $GenFitHiCBaseDir

#==================================
# write the sample configuration to a text file
#==================================
ConfFile=$GenFitHiCBaseDir/Configuration.txt
echo -e "Configurations used: \n Complete input alignment file : $INPFILE \n Method of interactions: $Method " > $ConfFile
if [[ $Method != 3 ]]; then
	if [[ -z $PeakFILE ]]; then
		echo -e "\n Peak file is computed from the input alignment using MACS2" >> $ConfFile
	else
		echo -e "\n Peak input file: $PeakFILE" >> $ConfFile
	fi
fi
echo -e "\n Bin size: $BIN_SIZE \n Interaction distance lower threshold: $LowDistThres \n Interaction distance higher threshold: $UppDistThres \n FitHiCBinMethod: $FitHiCBinMethod \n Output directory (results): $OutDir " >> $ConfFile

#==================================
# generate the timing profile text file (if the option is specified)
#==================================
if [ $Verbose == 1 ]; then
	OutTimeFile=$GenFitHiCBaseDir'/TimingProfile.txt'
	echo " ================ Time profiling =========== " > $OutTimeFile
	SECONDS=0
fi

#======================================================
# file depicting the interaction between different genomic segments
# initial list of interactions
Interaction_Initial_File=$GenFitHiCBaseDir/$PREFIX.anchors.interactions.initial.bed

# this file stores the coverage values for individual genomic intervals 
# with respect to fixed size bins
CoverageFile=$GenFitHiCBaseDir/$PREFIX.Coverage.bed

# file containing the final set of interactions
# either between different genomic segments
# or between peaks and non peaks
Interaction_File=$GenFitHiCBaseDir/$PREFIX.anchors.interactions.bed

# file containing the list of interactions along with different features for individual fixed size bins
# used to model the statistical significance
Complete_Interaction_File_Features=$GenFitHiCBaseDir/$PREFIX.Complete_InteractionFeatures.bed

# file containing the list of interactions
# sorted by the distance between interacting segments
InteractionFile_SortedGeneDist=$GenFitHiCBaseDir/$PREFIX.Compl_IntFeat.sortedGenDist.bed

#=====================================================
# Interaction between binned segments
# according to the input method and desired interaction type

# currently two output files are generated:
# 1) file containing the interactions among different genomic segments (Interaction_Initial_File and Interaction_File)
# 2) file containing coverage of different genomic segments (CoverageFile)
#=====================================================

# Peak to Peak, Peak to ALL or ALL to ALL interactions
if [[ $Method == 1 || $Method == 0 || $Method == 3 ]]; then
	# initial file - repeated interval pairs
	if [ ! -s $Interaction_Initial_File ]; then
		if [[ $Method == 1 || $Method == 0 ]]; then
			python ./src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -U $UppDistThres -L $LowDistThres -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile
		else
			python ./src/FitHiChIP_Interactions.py -M $Method -l $INPFILE -o $Interaction_Initial_File -U $UppDistThres -L $LowDistThres -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile
		fi
	fi
	if [ $Verbose == 1 ]; then
		duration=$SECONDS
		echo " ++++ Computing the interactions - time (in seconds): $duration" >> $OutTimeFile
		SECONDS=0
	fi
	
	# check if interaction count > 0
	# otherwise exit
	nic=`cat $Interaction_Initial_File | wc -l`
	if [[ $nic -eq 0 ]]; then
		echo 'Number of interaction = 0 - quit !!'
		exit 1
	fi 

	# duplicate pair of intervals - merge interactions
	# require sorted interactions
	if [ ! -s $Interaction_File ]; then
		if [[ $Method == 1 || $Method == 3 ]]; then
			# for 1) peak to peak and 2) all to all interactions, left side coordinate 
			# is always lower than the right side

			# Option 1 - choose the sum of contact counts for duplicate entries -  *** Incorrect ***
			# sort -k1,1 -k2,2n -k5,5n $Interaction_Initial_File | awk '{if (NR>1) {if (p1==$1 && p2==$2 && p3==$3 && p4==$4 && p5==$5 && p6==$6) {p7+=$7; REPT=1} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7; REPT=0;}}; p1=$1;p2=$2;p3=$3;p4=$4;p5=$5;p6=$6; {if (REPT==0) {p7=$7;}}} END {if (REPT==0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7}}' - > $Interaction_File

			# Option 2 - choose the maximum value of contact count for duplicate entries
			sort -k1,1 -k2,2n -k5,5n $Interaction_Initial_File | uniq | awk 'function max(x,y) {return x < y ? y : x} {if (NR>1) {if (p1==$1 && p2==$2 && p3==$3 && p4==$4 && p5==$5 && p6==$6) {p7=max($7,p7); REPT=1} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7; REPT=0;}}; p1=$1;p2=$2;p3=$3;p4=$4;p5=$5;p6=$6; {if (REPT==0) {p7=$7;}}} END {if (REPT==0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7}}' - > $Interaction_File
		else
			# for 1) peak to all, or 2) peak to non peak interactions, left side is the peak region

			# Special case: Peak to all interactions - when both intervals are peaks
			# Duplicate detection for the case L1-R1 and L2-R2 where L2=R1, and R2=L1
			# Solution: Swap according to coordinates - marker field 1 (swap) or 0 (no swap) as an extra 8th field
			# then merge the duplicates and revert back the entries
			awk -v x=0 -v y=1 '{if ($2 < $5) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"x} else {print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"\t"y}}' $Interaction_Initial_File | sort -k1,1 -k2,2n -k5,5n | awk 'function max(x,y) {return x < y ? y : x} {if (NR>1) {if (p1==$1 && p2==$2 && p3==$3 && p4==$4 && p5==$5 && p6==$6) {p7=max($7,p7); p8=max($8,p8); REPT=1} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7"\t"p8; REPT=0;}}; p1=$1;p2=$2;p3=$3;p4=$4;p5=$5;p6=$6; {if (REPT==0) {p7=$7; p8=$8;}}} END {if (REPT==0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7"\t"p8}}' - | awk '{if ($8 == 1) {print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' - > $Interaction_File
		fi
	fi
	if [ $Verbose == 1 ]; then
		duration=$SECONDS
		echo " ++++ Duplicate detect + merge pairs of interactions - time (in seconds): $duration" >> $OutTimeFile
		SECONDS=0
	fi

	# check if interaction count > 0
	# otherwise exit
	nic=`cat $Interaction_File | wc -l`
	if [[ $nic -eq 0 ]]; then
		echo 'Number of interaction = 0 - quit !!'
		exit 1
	fi 
fi

# Peak to non-Peak interactions
# LHS: peak region. RHS: non peak region
# interactions are unique - no duplicate entries
if [ $Method == 2 ]; then
	if [ ! -s $Interaction_File ]; then
		python ./src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_File -U $UppDistThres -L $LowDistThres -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile
	fi
	if [ $Verbose == 1 ]; then
		duration=$SECONDS
		echo " ++++ Computing the interactions - time (in seconds): $duration" >> $OutTimeFile
		SECONDS=0
	fi
	# check if interaction count > 0
	# otherwise exit
	nic=`cat $Interaction_File | wc -l`
	if [[ $nic -eq 0 ]]; then
		echo 'Number of interaction = 0 - quit !!'
		exit 1
	fi 
fi

#=======================================
# Analysis - plot the distribution of coverage for peaks and non-peaks
#if [ $DrawFig == 1 ]; then
#	Rscript ./Analysis/PlotCoverageDistr.r $CoverageFile $GenFitHiCBaseDir'/Other_Plots/'
#	if [ $Verbose == 1 ]; then
#		duration=$SECONDS		
#		echo " Plotting coverage of peaks - time (in seconds): $duration" >> $OutTimeFile
#		SECONDS=0
#	fi
#fi

#=======================================
# Extended interaction matrix (contact count) with bias features
# Currently bias is computed by normalized read depth
#=======================================
if [ ! -s $Complete_Interaction_File_Features ]; then
	Rscript ./src/Significance_Features.r -I $Interaction_File -E $CoverageFile -O $Complete_Interaction_File_Features
fi

if [ $Verbose == 1 ]; then
	duration=$SECONDS		
	echo " ++++ Computing extended interaction file with bias values - time (in seconds): $duration" >> $OutTimeFile
	SECONDS=0
fi

# derive the contact count column
cccol=`cat $Interaction_File | tail -n 1 | awk '{print NF}' -`
echo 'Contact count col: '$cccol

# derive the number of columns in the final interaction file (with bias features)
totcol=`cat $Complete_Interaction_File_Features | tail -n 1 | awk '{print NF}' -`
echo 'Total number of columns for the complete feature interactions: '$totcol

# sort the generated interaction file according to the distance between interacting regions (ascending order)
# and for equal distance, decreasing order of contact count
if [ ! -s $InteractionFile_SortedGeneDist ]; then
	coln=`expr $totcol + 1`
	awk -v OFS='\t' 'function abs(v) {return v < 0 ? -v : v} {print $0"\t"abs($5-$2)}' $Complete_Interaction_File_Features | sort -k${coln},${coln}n -k${cccol},${cccol}nr | cut -f1-${totcol} - > $InteractionFile_SortedGeneDist
fi

if [ $Verbose == 1 ]; then
	duration=$SECONDS
	echo " ++++ Sorting the interactions according to the distance between intervals: - time (in seconds): $duration" >> $OutTimeFile
	SECONDS=0
fi

# Plot the binned contact count vs genomic distance - consider equal length interval bins
if [ $DrawFig == 1 ]; then
	Rscript ./Analysis/binquantile.r $InteractionFile_SortedGeneDist $NBins $cccol $GenFitHiCBaseDir'/Other_Plots/BinQuantile.pdf'
fi

#=================================================================
# model the genomic distance vs contact count as a binomial distribution (Duan et. al. 2010)
# find out the statistical significant interactions
#=================================================================

if [ $BinomDistrModel == 1 ]; then
	if [ $Verbose == 1 ]; then
		SECONDS=0
	fi
	
	BinomDistrDir=$GenFitHiCBaseDir'/BinomDistr'
	mkdir -p $BinomDistrDir

	BinomDistr_IntFile=$BinomDistrDir/$PREFIX.interactions_BinomDistr.bed
	BinomDistr_Int_FiltFile=$BinomDistrDir/$PREFIX.interactions_BinomDistr_FILTER.bed
	BinomDistr_Int_Filt_PeakCount_File=$BinomDistrDir/$PREFIX.interactions_BinomDistr_FILTER.Peakcount.bed
	BinomDistr_PeakCCDistr_Text=$BinomDistrDir/$PREFIX.anchorPeakCCDistr.bed

	# Modeling the statistical significance by binomial distribution (Duan et. al. 2010) - main function
	if [ ! -s $BinomDistr_IntFile ]; then
		if [ $Verbose == 1 ]; then
			Rscript ./src/Interaction_Binom.r $InteractionFile_SortedGeneDist $BinomDistr_IntFile $NBins $cccol $OutTimeFile
		else
			Rscript ./src/Interaction_Binom.r $InteractionFile_SortedGeneDist $BinomDistr_IntFile $NBins $cccol
		fi
	fi

	if [ $Verbose == 1 ]; then
		duration=$SECONDS
		echo " ++++ Binomial distribution ---- Interactions + prob distr + P value computation + sorting : - time (in seconds): $duration" >> $OutTimeFile
		SECONDS=0
	fi

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	# The Q value is placed as the last column of the derived interaction file
	if [ ! -s $BinomDistr_Int_FiltFile ]; then
		awk -v q="$QVALUE" '{if ($NF < q && $NF > 0) {print $0}}' $BinomDistr_IntFile > $BinomDistr_Int_FiltFile
	fi

	# no of significant interactions (Q value based filtering)
	nsigbinom=`cat $BinomDistr_Int_FiltFile | wc -l`

	# check the no of significant interactions associated with individual peak segments 
	if [[ $Method != 3 ]]; then
		if [ ! -s $BinomDistr_Int_Filt_PeakCount_File ]; then
			# At least 10 significant interactions are required (an empirical threshold)
			if [[ $nsigbinom -gt 10 ]]; then
				cat $BinomDistr_Int_FiltFile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $BinomDistr_Int_Filt_PeakCount_File
			else
				echo 'number of significant interactions for binomial distribution < 10 - skip the peak count distribution function'
			fi
		fi
	fi

	if [[ $DrawFig == 1 && $Method != 3 ]]; then
		# Distribution of significant contact counts for individual peaks
		res_outdir=$BinomDistrDir'/Results'
		mkdir -p $res_outdir
		if [ -f $BinomDistr_Int_Filt_PeakCount_File ]; then
			if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltBinomDuan.pdf ] || [ ! -f $BinomDistr_PeakCCDistr_Text ]; then
				Rscript ./Analysis/ContactCountDistr.r $macs2dir/$PREFIX.anchors.bed $BinomDistr_Int_Filt_PeakCount_File $res_outdir/$PREFIX.PeakCCDistr_FiltBinomDuan.pdf $BinomDistr_PeakCCDistr_Text
			fi
		fi
	fi
fi

#==================================
# Model the statistical significance with FitHiC
# Note: we have two options - 1) use the probability only based on the contact count (BiasCorr = 0)
# 2) account for the bias factors (BiasCorr = 1)
#==================================

#for BiasCorr in 0 1; do
BiasCorr=0

	# according to the binning type (equal length vs equal occupancy)
	if [ $FitHiCBinMethod == 0 ]; then
		GenFitHiCDir=$GenFitHiCBaseDir'/FitHiC_EqLenBin'
	else
		GenFitHiCDir=$GenFitHiCBaseDir'/FitHiC_EqOccBin'
	fi
	if [ $BiasCorr == 1 ]; then
		GenFitHiCDir=$GenFitHiCDir'_BiasCorr'
	fi
	mkdir -p $GenFitHiCDir

	if [ $Verbose == 1 ]; then
		SECONDS=0
	fi

	FitHiC_Pass1_outfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1.bed
	FitHiC_Pass1_Filtfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.bed
	FitHiC_Pass1_Filt_PeakCountfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.Peakcount.bed
	FitHiC_Pass1_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_LogQ.bed
	FitHiC_Pass1_PeakCCDistr_Text=$GenFitHiCDir/$PREFIX.anchorPeakCCDistr.bed

	# Modeling the statistical significance by FitHiC - main function
	if [ ! -s $FitHiC_Pass1_outfile ]; then
		if [ $Verbose == 1 ]; then
			Rscript ./src/Interaction.r $InteractionFile_SortedGeneDist $FitHiC_Pass1_outfile $FitHiCBinMethod $BiasCorr $NBins $DrawFig $cccol $OutTimeFile 
		else
			Rscript ./src/Interaction.r $InteractionFile_SortedGeneDist $FitHiC_Pass1_outfile $FitHiCBinMethod $BiasCorr $NBins $DrawFig $cccol
		fi
	fi

	if [ $Verbose == 1 ]; then
		duration=$SECONDS
		echo " ++++ Spline Interactions + prob distr + P value computation + sorting : - time (in seconds): $duration" >> $OutTimeFile
		SECONDS=0
	fi

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	if [ ! -s $FitHiC_Pass1_Filtfile ]; then
		awk -v q="$QVALUE" '{if ($NF < q && $NF > 0) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
	fi

	# no of significant interactions (FitHiC)
	nsigFitHiC=`cat $FitHiC_Pass1_Filtfile | wc -l`

	# Check the no of significant contacts associated with each short peak segment
	if [[ $Method != 3 ]]; then
		if [ ! -s $FitHiC_Pass1_Filt_PeakCountfile ]; then
			# At least 10 significant interactions are required (empirical threshold)
			if [[ $nsigFitHiC -gt 10 ]]; then
				cat $FitHiC_Pass1_Filtfile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $FitHiC_Pass1_Filt_PeakCountfile
			else
				echo 'number of significant interactions for spline distribution < 10 - skip the peak count distribution function'
			fi
		fi
	fi

	# code for the second pass is not used for the moment - sourya

	# #==============================
	# # now the filtered interaction file (with respect to spline based binomial distribution)
	# # is applied for the second pass of the spline
	# # the modeled 2nd pass spline is applied on the original interaction bed file (without any p-value filtering)
	# #==============================
	# Rscript ./src/splinepass2.r $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass1_FILTER.bed $GenFitHiCDir/$PREFIX.anchors.interactions.sorted.bed $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass2.bed

	# # filter the generated 2nd pass spline bed file with respect to the P value
	# M=`cat $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass2.bed | wc -l`

	# if [ ! -f $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass2_FILTER.bed ]; then
	# 	awk -v m=$M '{if ($9 < (1/m)) {print $0}}' $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass2.bed > $GenFitHiCDir/$PREFIX.anchors.interactions.sorted_SplinePass2_FILTER.bed
	# fi

	#-----------------------------------
	# following scripts are for result statistics
	#-----------------------------------
	if [ $DrawFig == 1 ]; then
		# the R file takes the spline fitted interaction file (without q-value based filtering)
		# and returns the contact count distribution for two different sets of interactions
		# separated by the Q value threshold of 0.01
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 1 ]]; then
			Rscript ./Analysis/result_summary.r $FitHiC_Pass1_outfile $cccol $QVALUE
		else
			echo 'Number of significant spline interaction <= 1 - no result summary'
		fi
	fi

	# the filtered interaction (with respect to the spline) file is used to create a session file
	# for applying in WashU epigenome browser
	# for that, a special file containing only the interacting chromosome intervals 
	# and the log of Q value is created
	if [ ! -s $FitHiC_Pass1_LogQ_file ]; then
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 1 ]]; then
			Rscript ./Analysis/FiltSplineLogQ.r $FitHiC_Pass1_Filtfile $FitHiC_Pass1_LogQ_file
		else
			echo 'Number of significant spline interaction <= 1 - no session log of Q value'
		fi
	fi

	if [[ $DrawFig == 1 && $Method != 3 ]]; then
		res_outdir=$GenFitHiCDir'/Results'
		mkdir -p $res_outdir
		
		if [ -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
			# Distribution of significant contact counts (FitHiC) for individual peaks
			if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf ] || [ ! -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
				Rscript ./Analysis/ContactCountDistr.r $macs2dir/$PREFIX.anchors.bed $FitHiC_Pass1_Filt_PeakCountfile $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf $FitHiC_Pass1_PeakCCDistr_Text
			fi
		fi

		# Comparison of significant contact counts for both Binomial distribution and FitHiC
		# with respect to the peak segments (anchor segments from MACS2)
		if [ $BinomDistrModel == 1 ]; then
			if [ -f $BinomDistr_PeakCCDistr_Text ] && [ -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
				if [ ! -f $GenFitHiCDir/$PREFIX.PeakCC_Binom_FitHiC_Compare.txt ]; then
					Rscript ./Analysis/PeakBasisContactDistr.r $BinomDistr_PeakCCDistr_Text $FitHiC_Pass1_PeakCCDistr_Text $GenFitHiCDir/$PREFIX.PeakCC_Binom_FitHiC_Compare.txt
				fi
			fi
		fi

		# # For individual peaks, find and plot the contact counts 
		# # for both unfiltered and filtered (FDR) interactions
		# Rscript ./Analysis/ShortPeakCount.r $FitHiC_Pass1_outfile $cccol $QVALUE
	fi

	#====================================
	# compare the significant interactions between binomial distribution and FitHiC
	#====================================
	if [ $BinomDistrModel == 1 ]; then
		if [[ $nsigbinom -gt 1 && $nsigFitHiC -gt 1 ]]; then
			res_outdir=$GenFitHiCDir'/Results'
			mkdir -p $res_outdir
			if [ $FitHiCBinMethod == 0 ]; then
				if [ ! -f $res_outdir/Interaction_BinomDistr_SplineEqLen.pdf ]; then
					Rscript ./Analysis/Binned_Interaction_Overlap.r $BinomDistr_Int_FiltFile $FitHiC_Pass1_Filtfile $res_outdir/Interaction_BinomDistr_SplineEqLen.pdf $cccol 'BinomDistr' 'SplineEqLen' 
				fi
			else
				if [ ! -f $res_outdir/Interaction_BinomDistr_SplineEqOcc.pdf ]; then
					Rscript ./Analysis/Binned_Interaction_Overlap.r $BinomDistr_Int_FiltFile $FitHiC_Pass1_Filtfile $res_outdir/Interaction_BinomDistr_SplineEqOcc.pdf $cccol 'BinomDistr' 'SplineEqOcc'
				fi
			fi
		else
			echo 'At least one of the binomial or spline model has significant interactions count <= 1 - so no comparison between them is Possible'
		fi
	fi

#done

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir


