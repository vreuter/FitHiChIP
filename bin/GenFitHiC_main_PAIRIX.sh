#!/bin/bash

#===============
# FitHiChIP executable
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

Options:

  -- required:
       	-I  Inp				Input alignment file (complete alignment). The file can be either in BAM format, 
       						or in PAIRIX format. Conversion between BAM to PAIRIX can be done by 
       						using the utility bam2pairs (https://github.com/4dn-dcic/pairix)
		-M 	Met 			Can be 0, 1, 2, or 3 (default 3), where:
							3 (ALL): Finds the statistical significant interactions among all the binned segments.
							2 (Peak to Non Peak): Finds the statistical significant interactions 
							from the peaks to the non peak segments (subject to the binning employed).
								In the interaction file, Left side is a peak segment and right side is a non peak segment. 
								User must provide a peak detected file via -P option. 
								Otherwise, user must provide the complete alignment file in BAM format, in which case 
								the peaks would be computed from the alignment using MACS2.
							1 (Peak to Peak): Finds the statistical significant interactions between the peak regions 
							(subject to the binning employed).
								In the interaction file, both left and right sides are peak segments.  
								Possible only if either the user supplies the peak file via -P option, or 
								supplies the complete alignment file in BAM format.
							0 (Peak to ALL): Finds the statistical significant interactions from the peaks 
							(subject to the binning employed) to ALL fixed size binned segments.
								In the interaction file, Left side is a peak segment. 
								Right side can be either a peak or a non peak segment.
								Possible only if either the user supplies the peak file via -P option, or 
								supplies the complete alignment file in BAM format.
       	-P  Peak 			Input file containing the peak regions (obtained by the packages like MACS2). 
       						The file should be in .bed format.
       						If the user does not supply the peak file but provides the input alignment file 
       						in BAM format, peaks are computed using MACS2.
       						For alignment file in PAIRIX format, user must supply this peak file separately.
		-o  OutDir 			Output directory storing the output results
		-n  PREFIX          Prefix string of output files (Default = 'FitHiChIP').
		-B  BINSIZE        	Size of the bins [default = 5000], in bases. Contacts between a pair of bins are analyzed.
        -L  LowDistThr 	 	Lower distance threshold of significant interaction between two bins. Default = 20000 (means 20 Kb)
        -U  UppDistThr	 	Upper distance threshold of interaction between two bins. Default = 2000000 (means 2 Mb)
	   	-l 	biaslowthr 		Lower threshold of bias correction (fraction) - default 0.2
	   	-h  biashighthr 	Higher threshold of bias correction (fraction) - default 5
	   	-b  BeginBiasFilter 0/1: If 1, filters the interactions at the beginning according to the bias settings. Default 0.
	   	-e  EndBiasFilter 	0/1: If 1, computes the probabilities with respect to the bias (multiply). Default 0.
	   	-m  MergeNearInt	0/1: If 1, merges nearby (significant) interactions. Default 0.

  -- optional:
        -t  INT              Set number of threads for peak calling and finding the interactions [default = 1].
        -g  GSIZE			 Genome size parameter for MACS2 (if peak detection is sought). Default = 'hs'
        -q  QVALUE           Minimum FDR (q-value) cutoff for interaction detection [default = 1e-2].
        -v  TimeProf	 	 Specified as 1 or 0; if 1, generates a timing profile information
	 	-N  NBins			 In the FitHiC model, this is the max no of bins (equal occupancy or equal length) 
	 						 that is to be employed. Default 200.
	 	-D  Draw			 Specified as 1 or 0. If 1, draws the figures of various statistics / 
	 						 analysis. Default 0.

EOF
}

#==========================
# default parameters 
#==========================

# # number of threads employed for interaction / contact detection
# THREADS=1

# # default prefix
# PREFIX='FitHiChIP'

# # Bin size of 5 Kb resolution 
# BIN_SIZE=5000

# # genome size parameter for MACS2
# GSIZE='hs'

# # sample p and q values
# # PVALUE=1e-5
# QVALUE=1e-2

# # default output directory for storing the results
# OutDir=`pwd`

# # default option of TimeProf (output of the timing profile)
# TimeProf=0

# # default value of lower distance threhold between two interacting segments
# LowDistThres=20000

# # default value of upper distance threhold between two interacting segments
# UppDistThres=2000000

# # binning method for FitHiC technique
# # 0 - equal length bin, 1 for equal occupancy bin (default)
# FitHiCBinMethod=1

# # The file with peak detection output
# PeakFile=""

# if 1, computes the binomial distribution model of genomic distance vs contact count, 
# with respect to the paper Duan et. al.
# it is fixed as 0
BinomDistrModel=0

# # default value of the number of bins that is to be employed for FitHiC
# NBins=200

# # default value of method (type of interaction detection)
# Method=0

# # default value of plotting analysis figures
# DrawFig=0

# # default values for bias correction
# biaslowthr=0.2
# biashighthr=5
# BeginBiasFilter=0
# EndBiasFilter=0

# # default value for merging nearby interactions
# MergeInteraction=0

while getopts "I:P:t:n:B:q:o:v:L:U:N:M:D:g:l:h:b:e:m:" opt;
do
	case "$opt" in
		I) INPFILE=$OPTARG;;
		P) PeakFILE=$OPTARG;;
		n) PREFIX=$OPTARG;;
		t) THREADS=$OPTARG;;
		B) BIN_SIZE=$OPTARG;;
		q) QVALUE=$OPTARG;;
		o) OutDir=$OPTARG;;
		v) TimeProf=$OPTARG;;
		L) LowDistThres=$OPTARG;;
		U) UppDistThres=$OPTARG;;
		N) NBins=$OPTARG;;
		M) Method=$OPTARG;;
		D) DrawFig=$OPTARG;;
		g) GSIZE=$OPTARG;;
		l) biaslowthr=$OPTARG;;
		h) biashighthr=$OPTARG;;
		b) BeginBiasFilter=$OPTARG;;
		e) EndBiasFilter=$OPTARG;;
		m) MergeInteraction=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

# if [[ -z $INPFILE ]]; then
# 	echo 'No input alignment file is provided - exit !!'
# 	exit 1
# fi

# input file format - 1 means bam file, 2 means pairix file
extension="${INPFILE##*.}"
if [[ "$extension" == "bam" ]]; then
	inpfilefmt=1
else
	inpfilefmt=2
fi

# # check if the input peak file is not provided for the pairix file format
# # even for all to all interactions, we require a peak detection input file
# # if [[ $Method != 3 ]]; then
# 	if [[ -z $PeakFILE ]]; then
# 		if [[ $inpfilefmt == 2 ]]; then
# 			echo 'The input alignment is in PAIRIX format, but no peak detection file is provided as input - quit !!' 
# 			exit 1
# 		fi
# 	fi
# # fi

#==============================
# # convert the input files to their absolute path
# # input file
# INPFILE="$(cd "$(dirname "$INPFILE")"; pwd)/$(basename "$INPFILE")"
# # peak file, if provided
# if [[ ! -z $PeakFILE ]]; then
# 	PeakFILE="$(cd "$(dirname "$PeakFILE")"; pwd)/$(basename "$PeakFILE")"
# fi
# output directory - convert in the absolute format
mkdir -p $OutDir
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
# if user does not provide a peak detection output 
# but provides the bam alignment file
# we compute the peaks via MACS2
# applicable for every types of interactions
#==============================
# if [[ $Method != 3 ]]; then
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

			# earlier used MACS2 command
			# there we have retained all the duplicate peaks (peaks covering the same region)
			# macs2 callpeak -f AUTO -g $GSIZE --keep-dup all --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $INPFILE

			# now we use MACS2 command with the default values of --keep-dup
			# since it is recommended in the HiC pipeline
			macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $INPFILE -q $QVALUE
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
# fi

#==========================
# directory storing all the outputs
# 1) OutBaseDir: Stores the complete set of interactions (including normalization features)
# without any distance specific constraint
# 2) GenFitHiCBaseDir: stores the interactions (with normalization features) according to the 
# distance thresholds specified, sorted interaction file (according to increasing distance) 
# and also stores the FitHiC output for this particular settings
# 3) Separate directory underlying $GenFitHiCBaseDir storing the spline fit or binomial distribution
#==========================
if [ $Method == 3 ]; then
	OutBaseDir=$OutDir'/FitHiChIP_ALL2ALL'
elif [ $Method == 2 ]; then
	OutBaseDir=$OutDir'/FitHiChIP_Peak2NonPeak'
elif [ $Method == 1 ]; then
	OutBaseDir=$OutDir'/FitHiChIP_Peak2Peak'
else
	OutBaseDir=$OutDir'/FitHiChIP_Peak2ALL'
fi
OutBaseDir=$OutBaseDir'_b'$BIN_SIZE
echo "Base directory of generalized FitHiC: "$OutBaseDir
mkdir -p $OutBaseDir

GenFitHiCBaseDir=$OutBaseDir'/L_'$LowDistThres'_U_'$UppDistThres
mkdir -p $GenFitHiCBaseDir

#==================================
# write the sample configuration to a text file
#==================================
ConfFile=$GenFitHiCBaseDir/Configuration.txt

echo "Input alignment file: $INPFILE " > $ConfFile
echo "Method of interactions: $Method " >> $ConfFile
# if [[ $Method != 3 ]]; then
	if [[ -z $PeakFILE ]]; then
		echo "\n Peak file (computed using MACS2) : $macs2dir/$PREFIX.anchors.bed " >> $ConfFile
	else
		echo "\n Peak input file: $PeakFILE" >> $ConfFile
	fi
# fi
echo "OutDir: $OutDir " >> $ConfFile
echo "BIN_SIZE: $BIN_SIZE " >> $ConfFile
echo "LowDistThr: $LowDistThr " >> $ConfFile
echo "UppDistThr: $UppDistThr " >> $ConfFile
echo "QVALUE: $QVALUE " >> $ConfFile
echo "FitHiCBinMethod: $FitHiCBinMethod " >> $ConfFile
echo "NBins: $NBins " >> $ConfFile
echo "PREFIX: $PREFIX " >> $ConfFile
echo "DrawFig: $DrawFig " >> $ConfFile
echo "Timeprof: $TimeProf " >> $ConfFile
echo "Bias pre-filtering: $BeginBiasFilter " >> $ConfFile
echo "Prob Adjust due to bias: $EndBiasFilter " >> $ConfFile
echo "Bias lower cutoff: $biaslowthr " >> $ConfFile
echo "Bias higher cutoff: $biashighthr " >> $ConfFile
echo "Merging nearby interactions: $MergeInteraction " >> $ConfFile

# echo "RefGENOME: $RefGENOME " >> $ConfFile
# echo "ChrSizeFile: $ChrSizeFile " >> $ConfFile
# echo "MappabilityFile: $MappabilityFile " >> $ConfFile
# echo "RefFastaFile: $RefFastaFile " >> $ConfFile
# echo "REFragFile: $REFragFile " >> $ConfFile
# echo "GCContentWindowSize: $GCContentWindowSize " >> $ConfFile
# echo "MappabilityWindowSize: $MappabilityWindowSize " >> $ConfFile

#==================================
# generate the timing profile text file (if the option is specified)
#==================================
if [ $TimeProf == 1 ]; then
	OutTimeFile=$GenFitHiCBaseDir'/TimingProfile.txt'
	if [ ! -f $OutTimeFile ]; then
		echo " ================ Time profiling =========== " > $OutTimeFile
	fi
	start=$(date +%s.%N)
fi

#======================================================
# This directory stores the normalization related features	
FeatureDir=$OutBaseDir'/NormFeatures'
mkdir -p $FeatureDir

# this file stores the coverage values for individual genomic intervals 
# also stores whether individual genomic intervals are peak or not
# with respect to fixed size bins
CoverageFile=$FeatureDir'/'$PREFIX'.Coverage.bed'

# This file stores coverage + isPeak + bias values for individual genomic intervals
# The bias is computed as a ratio of coverage value (for that particular interval)
# with the mean coverage value
CoverageBiasFile=$FeatureDir'/'$PREFIX'.coverage_Bias.bed'

# This directory stores the full set of interactions and all the temporary stuff
# This full set of interaction does not have any distance constraint
# so, any distance specific interactions can be readily obtained
InteractionDir=$OutBaseDir'/Interaction_ALL'
mkdir -p $InteractionDir

# file depicting the interaction between different genomic segments
# initial list of interactions (without any distance threshold)
Interaction_Initial_File=$InteractionDir/$PREFIX.interactions_initial.bed

# file storing the temporary processed interactions (from the initial set)
# (without any distance threshold)
Interaction_Initial_TempFile=$InteractionDir/$PREFIX.interactions_initial_temp.bed

# file containing the final set of interactions
# either between different genomic segments
# or between peaks and non peaks (without any distance threshold)
Interaction_File=$InteractionDir/$PREFIX.interactions_FULL.bed

# file containing the list of interactions along with 
# different features for individual fixed size bins
# used to model the statistical significance
# (without any distance threshold)
Complete_Interaction_File_Features=$InteractionDir/$PREFIX.Interactions_FULL_ALLFeatures.bed

# file containing the list of interactions and the normalization related features
# with respect to the specified distance threshold
Complete_Interaction_File_Features_DistThr=$GenFitHiCBaseDir/$PREFIX.Interactions.DistThr.bed

# file containing the list of interactions
# sorted by the distance between interacting segments
InteractionFile_SortedGeneDist=$GenFitHiCBaseDir/$PREFIX.IntFeat.sortedGenDist.bed

#=====================================================
# Interaction between binned segments
# according to the input method and desired interaction type

# currently two output files are generated:
# 1) file containing the interactions among different genomic segments (Interaction_Initial_File and Interaction_File)
# 2) file containing coverage of different genomic segments (CoverageFile)

# Note: the interaction file generated does not have any distance threshold constraint
#=====================================================

# initial file - repeated interval pairs

# sourya - change this function such that 
# 1 at the end: original interaction had interaction1 > interaction2
# 0 at the end: original interaction had interaction1 < interaction2
# each interaction has 8 fields: 6 for 2 intervals, 1 for contact count, and 1 for this 1/0

if [ ! -s $Interaction_Initial_File ]; then
	if [[ $Method == 1 || $Method == 0 || $Method == 2 ]]; then
		# Peak to Peak, Peak to ALL, or peak to non peak interactions
		# python ../src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -U $UppDistThres -L $LowDistThres -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile

		# old version - computing the coverage simultaneously
		# python ../src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile

		# new version - only computing the interactions
		python ../src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -t $THREADS -b $BIN_SIZE -f $inpfilefmt
	else
		# all to all interactions

		# earlier implementation required no peak detection output file
		# no peak detection file is required
		# python ../src/FitHiChIP_Interactions.py -M $Method -l $INPFILE -o $Interaction_Initial_File -U $UppDistThres -L $LowDistThres -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile

		# now we use the peak detection output

		# old version - computing the coverage simultaneously
		# python ../src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -t $THREADS -b $BIN_SIZE -f $inpfilefmt -c $CoverageFile

		# new version - only computing the interactions
		python ../src/FitHiChIP_Interactions.py -p $macs2dir/$PREFIX.anchors.bed -M $Method -l $INPFILE -o $Interaction_Initial_File -t $THREADS -b $BIN_SIZE -f $inpfilefmt
	fi
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Computing the interactions (initial set) - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

# check if interaction count > 0
# otherwise exit
nic=`cat $Interaction_Initial_File | wc -l`
if [[ $nic -eq 0 ]]; then
	echo 'Number of interaction = 0 - quit !!'
	exit 1
fi

# sum of contact count for the same pair of intervals with the same indicator (1 or 0)
# Note: last awk statement is used to swap the columns 7 and 8
if [ ! -s $Interaction_Initial_TempFile ]; then
	awk -v OFS='\t' '{a[$1" "$2" "$3" "$4" "$5" "$6" "$8]+=$7}END{for (i in a){print i,a[i]}}' $Interaction_Initial_File | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' - | sort -k1,1 -k2,2n -k5,5n - > $Interaction_Initial_TempFile
fi

# now check if there exists repetition in the chromosome interval pair
# and only differs in the 8th column (1 or 0)
# in such a case, select the maximum contact count, and completely discard the 8th column 
# The generated interaction file has 7 columns:
# 6 columns define the interval pair, while the 7th denote the contact count
if [ ! -s $Interaction_File ]; then
	awk 'function max(x,y) {return x < y ? y : x} {if (NR>1) {if (p1==$1 && p2==$2 && p3==$3 && p4==$4 && p5==$5 && p6==$6) {p7=max($7,p7); REPT=1} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7; REPT=0;}}; p1=$1;p2=$2;p3=$3;p4=$4;p5=$5;p6=$6; {if (REPT==0) {p7=$7;}}} END {if (REPT==0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} else {print p1"\t"p2"\t"p3"\t"p4"\t"p5"\t"p6"\t"p7}}' $Interaction_Initial_TempFile > $Interaction_File
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Duplicate detect + merge pairs of interactions - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi

# =======================================
# compute the coverage for individual genomic bins
# Note: we always have the peak detection file
# =======================================
if [ ! -f $CoverageFile ]; then
	python ../src/CoverageBin.py -p $macs2dir/$PREFIX.anchors.bed -i $INPFILE -b $BIN_SIZE -f $inpfilefmt -o $CoverageFile
fi

# =======================================
# Compute the bias for peaks and non-peaks separately 
# Curerntly, the bias values are written in the coverage file itself (overwriting)
# and it uses only the coverage values information
# later, more complex models can be incorporated
# =======================================
if [ ! -f $CoverageBiasFile ]; then
	Rscript ../src/BiasCalc.r --CoverageFile $CoverageFile --OutFile $CoverageBiasFile
	echo 'Appended bias information for individual genomic bins'
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing bias of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi

# =======================================
# Analysis - plot the distribution of coverage for peaks and non-peaks separately
if [ $DrawFig == 1 ]; then
	Rscript ../Analysis/PlotCoverageDistr.r $CoverageFile $GenFitHiCBaseDir'/Plots/'
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)		
		echo " Plotting coverage of peaks - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================================
# Extended interaction matrix (contact count) with bias features
# Currently bias is computed by normalized read depth
# Note: The interaction file with all the normalization features
# do not have any distance constraint - they can be used for multiple distance ranges
#=======================================
if [ ! -s $Complete_Interaction_File_Features ]; then
	Rscript ../src/Significance_Features.r -I $Interaction_File -E $CoverageBiasFile -O $Complete_Interaction_File_Features
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)		
		echo " ++++ Computing extended interaction file with bias values - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================================
# Generate interactions (with normalization features) 
# for the mentioned distance threshold values
# ****************** 
# Also, this routine checks if any interaction is from non-peak to peak, 
# i.e. left side is non peak and right side is peak
# in such a case, the content is swapped between these two segments
# 9th and 12th fields denote peak information
# ****************** 
#=======================================
if [ ! -f $Complete_Interaction_File_Features_DistThr ]; then
	awk -v l="$LowDistThres" -v u="$UppDistThres" 'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || ($1==$4 && abs($2-$5)>=l && abs($2-$5)<=u)) {print $0}}' $Complete_Interaction_File_Features | awk '{if (NR>1 && $9==0 && $12==1) {print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"\t"$11"\t"$12"\t"$13"\t"$8"\t"$9"\t"$10} else {print $0}}' - > $Complete_Interaction_File_Features_DistThr
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions within distance thresholds: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================================
# derive the contact count column
cccol=`cat $Interaction_File | tail -n 1 | awk '{print NF}' -`
echo 'Contact count col: '$cccol

# derive the number of columns in the final interaction file (with bias features)
totcol=`cat $Complete_Interaction_File_Features_DistThr | tail -n 1 | awk '{print NF}' -`
echo 'Total number of columns for the complete feature interactions: '$totcol

# sort the generated interaction file according to the 
# distance between interacting regions (ascending order)
# and for equal distance, decreasing order of contact count
if [ ! -s $InteractionFile_SortedGeneDist ]; then
	coln=`expr $totcol + 1`
	awk -v OFS='\t' 'function abs(v) {return v < 0 ? -v : v} {print $0"\t"abs($5-$2)}' $Complete_Interaction_File_Features_DistThr | sort -k${coln},${coln}n -k${cccol},${cccol}nr | cut -f1-${totcol} - > $InteractionFile_SortedGeneDist
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Sorting the interactions according to the distance between intervals: - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi

#=================================================================
# model the genomic distance vs contact count as a binomial distribution (Duan et. al. 2010)
# find out the statistical significant interactions
#=================================================================
if [ $BinomDistrModel == 1 ]; then
	
	BinomDistrDir=$GenFitHiCBaseDir'/BinomDistr'
	mkdir -p $BinomDistrDir

	BinomDistr_IntFile=$BinomDistrDir/$PREFIX.interactions_BinomDistr.bed
	BinomDistr_Int_FiltFile=$BinomDistrDir/$PREFIX.interactions_BinomDistr_FILTER.bed
	BinomDistr_Int_Filt_PeakCount_File=$BinomDistrDir/$PREFIX.interactions_BinomDistr_FILTER.Peakcount.bed
	BinomDistr_PeakCCDistr_Text=$BinomDistrDir/$PREFIX.anchorPeakCCDistr.bed

	# Modeling the statistical significance by binomial distribution (Duan et. al. 2010) - main function
	if [ ! -s $BinomDistr_IntFile ]; then
		if [ $TimeProf == 1 ]; then
			Rscript ../src/Interaction_Binom.r $InteractionFile_SortedGeneDist $BinomDistr_IntFile $NBins $cccol $OutTimeFile
		else
			Rscript ../src/Interaction_Binom.r $InteractionFile_SortedGeneDist $BinomDistr_IntFile $NBins $cccol
		fi
		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Binomial distribution ---- Interactions + prob distr + P value computation + sorting : - time (in seconds): $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi
	fi

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	# The Q value is placed as the last column of the derived interaction file
	if [ ! -s $BinomDistr_Int_FiltFile ]; then
		awk -v q="$QVALUE" '{if ((NR==1) || ($NF < q && $NF > 0)) {print $0}}' $BinomDistr_IntFile > $BinomDistr_Int_FiltFile
	fi

	# no of significant interactions (Q value based filtering)
	nsigbinom=`cat $BinomDistr_Int_FiltFile | wc -l`

	# check the no of significant interactions associated with individual peak segments 
	if [[ $Method != 3 ]]; then
		if [ ! -s $BinomDistr_Int_Filt_PeakCount_File ]; then
			# At least 10 significant interactions are required (an empirical threshold)
			if [[ $nsigbinom -gt 11 ]]; then
				awk 'NR>1' $BinomDistr_Int_FiltFile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $BinomDistr_Int_Filt_PeakCount_File
			else
				echo 'number of significant interactions for binomial distribution < 10 - skip the peak count distribution function'
			fi
		fi
	fi

	if [[ $DrawFig == 1 && $Method != 3 ]]; then
		# Distribution of significant contact counts for individual peaks
		# applicable except ALL to ALL interactions
		res_outdir=$BinomDistrDir'/Results'
		mkdir -p $res_outdir
		if [ -f $BinomDistr_Int_Filt_PeakCount_File ]; then
			if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltBinomDuan.pdf ] || [ ! -f $BinomDistr_PeakCCDistr_Text ]; then
				Rscript ../Analysis/ContactCountDistr.r $macs2dir/$PREFIX.anchors.bed $BinomDistr_Int_Filt_PeakCount_File $res_outdir/$PREFIX.PeakCCDistr_FiltBinomDuan.pdf $BinomDistr_PeakCCDistr_Text
			fi
		fi
	fi
fi

#==================================
# Model the statistical significance with FitHiC
# Note: we have two options - 1) use the probability only based on the contact count (BiasCorr = 0)
# 2) account for the bias factors (BiasCorr = 1)
#==================================
# according to the binning type (equal length vs equal occupancy)
if [ $FitHiCBinMethod == 0 ]; then
	GenFitHiCDir=$GenFitHiCBaseDir'/FitHiC_EqLenBin'
else
	GenFitHiCDir=$GenFitHiCBaseDir'/FitHiC_EqOccBin'
fi
if [[ $BeginBiasFilter == 0 && $EndBiasFilter == 0 ]]; then
	BiasCorr=0
else		
	GenFitHiCDir=$GenFitHiCDir'_BiasCorr_'$biaslowthr'_'$biashighthr'_b'$BeginBiasFilter'_e'$EndBiasFilter
	BiasCorr=1
fi
mkdir -p $GenFitHiCDir

# FitHiC generated output files
FitHiC_Pass1_outfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1.bed
FitHiC_Pass1_Filtfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.bed
FitHiC_Bin_Filt_CC_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.BinSpecific.bed
FitHiC_Pass1_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_LogQ.bed
FitHiC_Pass1_PeakCCDistr_Text=$GenFitHiCDir/$PREFIX.PeakCCDistr.bed
FitHiC_Pass1_Filt_MergedIntfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED.bed
FitHiC_Pass1_Filt_MergedInt_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED_LogQ.bed

#====================================
# Modeling the statistical significance by FitHiC - main function
#====================================
if [ ! -f $FitHiC_Pass1_outfile ]; then
	if [ $FitHiCBinMethod == 0 ]; then
		# equal length binning (last option)
		Rscript ../src/Interaction.r --InpFile $InteractionFile_SortedGeneDist --OutFile $FitHiC_Pass1_outfile --Norm $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --EqLenBin
	else
		# equal occupancy binning - default settings
		Rscript ../src/Interaction.r --InpFile $InteractionFile_SortedGeneDist --OutFile $FitHiC_Pass1_outfile --Norm $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter
	fi
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Applying FitHiC: - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi
echo 'Applied FitHiC'

#====================================
# Filter the interaction file with respect to significance (Q value < $QVALUE)
#====================================
if [ ! -s $FitHiC_Pass1_Filtfile ]; then
	awk -v q="$QVALUE" '{if ((NR==1) || ($NF < q && $NF > 0)) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
fi

echo 'Extracted significant interactions from FitHiC'

# no of significant interactions (FitHiC)
nsigFitHiC=`cat $FitHiC_Pass1_Filtfile | wc -l`

if [ $DrawFig == 1 ]; then
	# the R file takes the spline fitted interaction file (without q-value based filtering)
	# and returns the contact count distribution for two different sets of interactions
	# separated by the Q value threshold of 0.01
	# check for non empty interactions file
	if [[ $nsigFitHiC -gt 1 ]]; then
		Rscript ../Analysis/result_summary.r $FitHiC_Pass1_outfile $cccol $QVALUE
	else
		echo 'Number of significant spline interaction <= 1 - no result summary'
	fi
fi

#====================================
# the filtered interaction (with respect to the spline) file is used to create a session file
# for applying in WashU epigenome browser
# for that, a special file containing only the interacting chromosome intervals 
# and the log of Q value is created
#====================================
if [ ! -s $FitHiC_Pass1_LogQ_file ]; then
	# check for non empty interactions file
	if [[ $nsigFitHiC -gt 2 ]]; then
		awk '{if (NR > 1) {print $1","$2","$3"\t"$4","$5","$6"\t"(-log($NF)/log(10))}}' $FitHiC_Pass1_Filtfile > $FitHiC_Pass1_LogQ_file
	else
		echo 'There is no significant interaction - so no WashU specific session file is created !!'
	fi
fi

#===========================================
# Now check individual bins (with respect to the given bin size)
# and find the number of (significant) contacts each bin is associated
#===========================================
if [[ $Method != 3 ]]; then
	if [ ! -s $FitHiC_Bin_Filt_CC_file ]; then
		# At least 10 significant interactions are required (empirical threshold)
		if [[ $nsigFitHiC -gt 11 ]]; then
			awk 'NR>1' $FitHiC_Pass1_Filtfile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $FitHiC_Bin_Filt_CC_file
		else
			echo 'number of significant interactions for spline distribution < 10 - skip the bin specific significant contact count distribution function'
		fi
	fi
fi

# Map the significant contact count information with respect to the original peaks
# and plot the distribution of peaks vs significant contacts
if [[ $DrawFig == 1 && $Method != 3 ]]; then
	res_outdir=$GenFitHiCDir'/Results'
	mkdir -p $res_outdir
	if [ -f $FitHiC_Bin_Filt_CC_file ]; then
		# Distribution of significant contact counts (FitHiC) for individual peaks
		# applicable except ALL to ALL interactions
		if [ ! -f $res_outdir/PeakCCDistr_FiltSpline.pdf ] || [ ! -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
			Rscript ../Analysis/ContactCountDistr.r $macs2dir/$PREFIX.anchors.bed $FitHiC_Bin_Filt_CC_file $res_outdir/PeakCCDistr_FiltSpline.pdf $FitHiC_Pass1_PeakCCDistr_Text
		fi
	fi

# 	# Comparison of significant contact counts for both Binomial distribution and FitHiC
# 	# with respect to the peak segments (anchor segments from MACS2)
# 	if [ $BinomDistrModel == 1 ]; then
# 		if [ -f $BinomDistr_PeakCCDistr_Text ] && [ -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
# 			if [ ! -f $GenFitHiCDir/$PREFIX.PeakCC_Binom_FitHiC_Compare.txt ]; then
# 				Rscript ../Analysis/PeakBasisContactDistr.r $BinomDistr_PeakCCDistr_Text $FitHiC_Pass1_PeakCCDistr_Text $GenFitHiCDir/$PREFIX.PeakCC_Binom_FitHiC_Compare.txt
# 			fi
# 		fi
# 	fi

# 	# # For individual peaks, find and plot the contact counts 
# 	# # for both unfiltered and filtered (FDR) interactions
# 	# Rscript ../Analysis/ShortPeakCount.r $FitHiC_Pass1_outfile $cccol $QVALUE

fi

# #====================================
# # compare the significant interactions between binomial distribution and FitHiC
# #====================================
# if [ $BinomDistrModel == 1 ]; then
# 	if [[ $nsigbinom -gt 1 && $nsigFitHiC -gt 1 ]]; then
# 		res_outdir=$GenFitHiCDir'/Results'
# 		mkdir -p $res_outdir
# 		if [ $FitHiCBinMethod == 0 ]; then
# 			if [ ! -f $res_outdir/Interaction_BinomDistr_SplineEqLen.pdf ]; then
# 				Rscript ../Analysis/Binned_Interaction_Overlap.r $BinomDistr_Int_FiltFile $FitHiC_Pass1_Filtfile $res_outdir/Interaction_BinomDistr_SplineEqLen.pdf $cccol 'BinomDistr' 'SplineEqLen' 
# 			fi
# 		else
# 			if [ ! -f $res_outdir/Interaction_BinomDistr_SplineEqOcc.pdf ]; then
# 				Rscript ../Analysis/Binned_Interaction_Overlap.r $BinomDistr_Int_FiltFile $FitHiC_Pass1_Filtfile $res_outdir/Interaction_BinomDistr_SplineEqOcc.pdf $cccol 'BinomDistr' 'SplineEqOcc'
# 			fi
# 		fi
# 	else
# 		echo 'At least one of the binomial or spline model has significant interactions count <= 1 - so no comparison between them is Possible'
# 	fi
# fi

if [ $TimeProf == 1 ]; then
	duration=$(echo "$(date +%s.%N) - $start" | bc)
	echo " ++++ Time (in seconds) for post processing FitHiC results: $duration" >> $OutTimeFile
	start=$(date +%s.%N)
fi

#=============================
# if merging nearby interactions are enabled
# then we merge the nearby interactions from the earlier generated significant interactions
# and also create a washu browser generated compatible file
#=============================
if [ $MergeInteraction == 1 ]; then
	if [ ! -f $FitHiC_Pass1_Filt_MergedIntfile ]; then
		python ../src/CombineNearbyInteraction.py --InpFile $FitHiC_Pass1_Filtfile --OutFile $FitHiC_Pass1_Filt_MergedIntfile --headerInp 1 --binsize $BIN_SIZE
	fi
	if [ ! -s $FitHiC_Pass1_Filt_MergedInt_LogQ_file ]; then
		nint=`cat $FitHiC_Pass1_Filt_MergedIntfile | wc -l`
		if [[ $nint -gt 2 ]]; then
			# 9th field stores the Q value
			awk '{if (NR > 1) {print $1","$2","$3"\t"$4","$5","$6"\t"(-log($9)/log(10))}}' $FitHiC_Pass1_Filt_MergedIntfile > $FitHiC_Pass1_Filt_MergedInt_LogQ_file
		else
			echo 'There is no significant interaction - so no WashU specific session file is created !!'
		fi
	fi
	echo 'Merged nearby significant interactions - created washu browser compatible file for these merged interactions!!!'
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for merging nearby interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir

