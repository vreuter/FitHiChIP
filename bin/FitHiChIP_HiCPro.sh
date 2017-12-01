#!/bin/bash

#===============
# FitHiChIP executable
# when the input is a valid pairs (txt or .txt.gz) file from the HiC-pro pipeline 

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# usage info
usage(){
cat << EOF

Options:

  -- required:
   	-I  InpValidPairsFile		Input HiC-pro valid pairs (.txt / txt.gz) file
   	-T  InpBinIntervalFile		Input bin interval file obtained from the valid pairs file
   	-M  InpMatrixFile			Input contact matrix file (obtained from HiC-pro utility)
   	-P  PeakFILE				Bed file containing Peak information 
   	-o  OutDir 					Base directory storing output results.
   	-c  ChrSizeFile 			File containing the size of reference chromosome
	-B  BINSIZE        			Size of the bins [default = 5000], in bases, for detecting contacts
    -L  LowDistThres 	 	 	Lower distance threshold of contacts between two segments (default = 20000)
    -U  UppDistThres	 	 	Upper distance threshold of contacts between two segments (default = 2000000)
    -H  HiCProMatrixBuild		Executable of HiC-pro pipeline matrix building utility
    -n  PREFIX          		Prefix string of output files (Default = 'FitHiChIP')
   	-l 	biaslowthr 				Lower threshold of bias correction (fraction) - default 0.2
   	-h  biashighthr 			Higher threshold of bias correction (fraction) - default 5
   	-b  BeginBiasFilter 		0/1: If 1, filters the interactions at the beginning according to the bias settings. Default 0.
   	-e  EndBiasFilter 			0/1: If 1, computes the probabilities with respect to the bias (multiply). Default 0.
   	-m  MergeNearInt			0/1: If 1, merges nearby (significant) interactions. Default 0.
	
  -- optional:
    -q  QVALUE           		Minimum FDR (q-value) cutoff for interaction detection [default = 1e-2].
    -v  TimeProf	 		 	Specified as 1 or 0; if 1, generates a timing profile information
 	-D  DrawFig			 		Specified as 1 or 0. If 1, draws the figures of various statistics / analysis. Default 0.
 	-N  NBins			 		In FitHiC, this is the max no of bins (equal occupancy) to be employed. Default 200.

EOF
}

# # input files
# InpValidPairsFile=""
# InpBinIntervalFile=""
# InpMatrixFile=""
# PeakFILE=""

# PREFIX='FitHiChIP'

# # size of the chromosome that is to be provided
# ChrSizeFile=""

# # 5 Kb resolution
# BIN_SIZE=5000

# QVALUE=1e-2	# 0.01

# # default value of output directory is the present working directory
# OutDir=`pwd`'/'

# # lower distance threshold for two cis interactions
# LowDistThres=20000	# 20 Kb

# # upper distance threshold for two cis interactions
# UppDistThres=2000000 # 2 Mb

# # number of bins employed for FitHiC
# NBins=200

# # default value of plotting analysis figures
# DrawFig=0

# # binning method for FitHiC technique
# # 1 for equal occupancy bin (default)
# FitHiCBinMethod=1

# # option to note down the timing information
# TimeProf=0

# #=========================
# # bias correction related parameters
# #=========================

# # default lower cutoff of bias value
# biaslowthr=0.2

# # default higher cutoff of bias value
# biashighthr=5

# # boolean variable for pre-filtering the interactions according to the bias value
# BeginBiasFilter=0

# # boolean variable for probability adjustment of the interactions according to the bias value
# EndBiasFilter=0

# # Merging interactions which are near to each other
# MergeInteraction=1

while getopts "I:T:M:P:o:c:B:L:U:q:N:H:n:D:v:l:h:b:e:m:" opt;
do
	case "$opt" in
		I) InpValidPairsFile=$OPTARG;;
		T) InpBinIntervalFile=$OPTARG;;
		M) InpMatrixFile=$OPTARG;;
		P) PeakFILE=$OPTARG;;
		o) OutDir=$OPTARG;;
		c) ChrSizeFile=$OPTARG;;
		B) BIN_SIZE=$OPTARG;;
		L) LowDistThres=$OPTARG;;
		U) UppDistThres=$OPTARG;;
		q) QVALUE=$OPTARG;;
		N) NBins=$OPTARG;;
		H) HiCProMatrixBuild=$OPTARG;;
		n) PREFIX=$OPTARG;;
		D) DrawFig=$OPTARG;;
		v) TimeProf=$OPTARG;;
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

# create the output directory
mkdir -p $OutDir

#============================
# print the parameters and values
#============================
ConfFile=$OutDir/Parameters.txt

echo "InpValidPairsFile: $InpValidPairsFile " > $ConfFile
echo "InpBinIntervalFile: $InpBinIntervalFile " >> $ConfFile
echo "InpMatrixFile: $InpMatrixFile " >> $ConfFile
echo "PeakFILE: $PeakFILE " >> $ConfFile
echo "OutDir: $OutDir " >> $ConfFile
echo "ChrSizeFile: $ChrSizeFile " >> $ConfFile
echo "BIN_SIZE: $BIN_SIZE " >> $ConfFile
echo "LowDistThr: $LowDistThres " >> $ConfFile
echo "UppDistThr: $UppDistThres " >> $ConfFile
echo "QVALUE: $QVALUE " >> $ConfFile
echo "NBins: $NBins " >> $ConfFile
echo "HiCProMatrixBuild executable: $HiCProMatrixBuild " >> $ConfFile
echo "PREFIX: $PREFIX " >> $ConfFile
echo "DrawFig: $DrawFig " >> $ConfFile
echo "Timeprof: $TimeProf " >> $ConfFile
echo "Bias pre-filtering: $BeginBiasFilter " >> $ConfFile
echo "Prob Adjust due to bias: $EndBiasFilter " >> $ConfFile
echo "Bias lower cutoff: $biaslowthr " >> $ConfFile
echo "Bias higher cutoff: $biashighthr " >> $ConfFile
echo "Merging nearby interactions: $MergeInteraction " >> $ConfFile

#=======================================
# generate a file which will contain the timing profile
if [ $TimeProf == 1 ]; then
	OutTimeFile=$OutDir'/TimingProfile.txt'
	echo " ================ Time profiling =========== " > $OutTimeFile
	start=$(date +%s.%N)
fi

#==============================
# important - sourya
# first change the current working directory to the directory containing this script
# it is useful when the script is invoked from a separate directory
#==============================
currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir

# generate the matrix of Hi-C interactions (ALL)
# using HiC-pro pipeline
HiCProMatrixDir=$OutDir'/HiCPro_Matrix_BinSize'$BIN_SIZE
mkdir -p $HiCProMatrixDir

#=================
# if the matrices are not provided and the validpairs text file is provided
# then compute the interaction matrices using the HiC-pro utility
#=================
if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
	echo '*** Computing HiC-pro matrices from the input valid pairs file'

	# This directory and prefix is used to denote the generated matrices
	# using the HiC pro routine
	OutPrefix=$HiCProMatrixDir'/MatrixHiCPro'

	if [ ! -f $OutPrefix'_abs.bed' ]; then
		# check the extension of input valid pairs file
		# and extract accordingly
		if [[ $InpValidPairsFile == *.gz ]]; then
			zcat $InpValidPairsFile | $HiCProMatrixBuild --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper'  
		else
			cat $InpValidPairsFile | $HiCProMatrixBuild --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
		fi
	fi

	# now assign the matrix names to the designated variables
	InpBinIntervalFile=$OutPrefix'_abs.bed'
	InpMatrixFile=$OutPrefix'.matrix'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing the interaction matrix using HiC-Pro build_matrix utility: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================
# Now generate the list of interactions from the HiC-pro matrix data
# ALL to ALL interactions
# Both cis and trans interactions are considered (with respect to the given bin size)
# No distance threshold based filtering is used
# Interaction format:
# chr1	start1	end1	chr2	start2	end2	cc
#=======================
Interaction_Initial_File=$HiCProMatrixDir/$PREFIX.interactions.initial.bed
if [ ! -f $Interaction_Initial_File ]; then
	Rscript ../src/InteractionHicPro.r $InpBinIntervalFile $InpMatrixFile $Interaction_Initial_File
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================
# generate filtered cis interactions 
# with respect to distance thresholds
# ALL to ALL interactions
#=======================
# create a directory for individual distance thresholds
InteractionThrDir=$HiCProMatrixDir'/L_'$LowDistThres'_U'$UppDistThres
mkdir -p $InteractionThrDir
Interaction_File=$InteractionThrDir/$PREFIX.cis.interactions.DistThr.bed

if [ ! -f $Interaction_File ]; then
	awk -v l="$LowDistThres" -v u="$UppDistThres" 'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || ($1==$4 && abs($2-$5)>=l && abs($2-$5)<=u)) {print $0}}' $Interaction_Initial_File > $Interaction_File

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions within distance thresholds: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#============================
# this directory stores the features and associated data
#============================
FeatureDir=$OutDir'/NormFeatures'
mkdir -p $FeatureDir

#=================
# From the input valid paired end read file, and the given input bin size parameter
# compute the coverage of individual genome segments (bins)
# the output is a list of chromosome and bins, their coverage, and a boolean indicator 
# whether the segment overlaps with a peak
#=================
CoverageFile=$FeatureDir'/'$PREFIX'.coverage.bed'
if [ ! -f $CoverageFile ]; then
	python ../src/CoverageBin.py -i $InpValidPairsFile -p $PeakFILE -b $BIN_SIZE -o $CoverageFile -c $ChrSizeFile
	echo 'Computed initial coverage information for individual genomic bins'
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting coverage of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

# =======================================
# Compute the bias for peaks and non-peaks separately 
# Curerntly, the bias values are written in the coverage file itself (overwriting)
# and it uses only the coverage values information
# later, more complex models can be incorporated
# =======================================
CoverageBiasFile=$FeatureDir'/'$PREFIX'.coverage_Bias.bed'
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
	Rscript ../Analysis/PlotCoverageDistr.r $CoverageFile $FeatureDir'/Plots/'
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)		
		echo " Plotting coverage of peaks - time (in seconds): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=================
# templates of filenames used throughout the execution
#=================
InteractionFileName='Interactions.bed'
InteractionSortedDistFileName='Interactions.sortedGenDist.bed'

#============================
# create all to all interaction file
# with the features like read depth, mappability, GC content, 
# and number of RE sites
#============================
DirALLtoALL=$OutDir'/FitHiChIP_ALL2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
mkdir -p $DirALLtoALL
IntFileALLtoALL=$DirALLtoALL'/'$InteractionFileName

if [ ! -f $IntFileALLtoALL ]; then
	Rscript ../src/Significance_Features.r -I $Interaction_File -E $CoverageBiasFile -O $IntFileALLtoALL
fi

# derive the contact count column
cccol=`cat $Interaction_File | tail -n 1 | awk '{print NF}' -`
echo 'Contact count col: '$cccol

# derive the number of columns in the interaction file with normalization 
# related features
totcol=`cat $IntFileALLtoALL | tail -n 1 | awk '{print NF}' -`
echo 'Total number of columns for the complete feature interactions: '$totcol

if [ $TimeProf == 1 ]; then
	duration=$(echo "$(date +%s.%N) - $start" | bc)
	echo " ++++ Time (in seconds) for computing pairwise interactions among all segments: $duration" >> $OutTimeFile
	start=$(date +%s.%N)
fi

#===================
# Create the interaction files for other types of interactions
# peak to peak, peak to non peak, and peak to all
#===================
DirPeaktoPeak=$OutDir'/FitHiChIP_Peak2Peak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
DirPeaktoNonPeak=$OutDir'/FitHiChIP_Peak2NonPeak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
DirPeaktoALL=$OutDir'/FitHiChIP_Peak2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
mkdir -p $DirPeaktoPeak
mkdir -p $DirPeaktoNonPeak
mkdir -p $DirPeaktoALL

IntFilePeaktoPeak=$DirPeaktoPeak'/'$InteractionFileName
IntFilePeaktoNonPeak=$DirPeaktoNonPeak'/'$InteractionFileName
IntFilePeaktoALL=$DirPeaktoALL'/'$InteractionFileName

if [ ! -f $IntFilePeaktoPeak ]; then
	# peak to peak interactions 
	# 9th and 15th fields are 1
	awk '(NR==1) || ($9==1 && $12==1)' $IntFileALLtoALL > $IntFilePeaktoPeak
fi

if [ ! -f $IntFilePeaktoNonPeak ]; then
	# peak to non peak interactions
	# 9th field is 1, but 15th field is 0
	awk '(NR==1) || ($9==1 && $12==0)' $IntFileALLtoALL > $IntFilePeaktoNonPeak
fi

if [ ! -f $IntFilePeaktoALL ]; then
	# peak to all interactions
	# just check if 9th field is 1
	awk '(NR==1) || ($9==1)' $IntFileALLtoALL > $IntFilePeaktoALL
fi

if [ $TimeProf == 1 ]; then
	duration=$(echo "$(date +%s.%N) - $start" | bc)
	echo " ++++ Time (in seconds) for assigning different types of interactions: $duration" >> $OutTimeFile
	start=$(date +%s.%N)
fi

#==============================
# navigate through individual types of interactions (corresponding folders)
# and apply FitHiC for individual interaction types
#==============================
# $DirALLtoALL is commented for the moment - sourya
for dirname in $DirPeaktoPeak $DirPeaktoNonPeak $DirPeaktoALL $DirALLtoALL; do
	
	echo 'Processing the directory: '$dirname

	#==============
	# first create interaction files with sorted genomic distance
	#==============
	CurrIntFile=$dirname'/'$InteractionFileName
	CurrIntFileSortDist=$dirname'/'$InteractionSortedDistFileName

	coln=`expr $totcol + 1`
	if [ ! -s $CurrIntFileSortDist ]; then
		awk -v OFS='\t' 'function abs(v) {return v < 0 ? -v : v} {print $0"\t"abs($5-$2)}' $CurrIntFile | sort -k${coln},${coln}n -k${cccol},${cccol}nr | cut -f1-${totcol} - > $CurrIntFileSortDist	
	fi

	echo 'Created sorted genomic distance based interaction file'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for sorting the interactions (according to genomic distance): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	#==============
	# now apply FitHiC on the sorted gene distance based interaction matrix
	#==============
	GenFitHiCDir=$dirname'/FitHiC_EqOccBin'
	if [[ $BeginBiasFilter == 0 && $EndBiasFilter == 0 ]]; then
		BiasCorr=0
	else		
		GenFitHiCDir=$GenFitHiCDir'_BiasCorr_'$biaslowthr'_'$biashighthr'_b'$BeginBiasFilter'_e'$EndBiasFilter
		BiasCorr=1
	fi
	mkdir -p $GenFitHiCDir

	FitHiC_Pass1_outfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1.bed
	FitHiC_Pass1_Filtfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.bed
	FitHiC_Pass1_Filt_PeakCountfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.Peakcount.bed
	FitHiC_Pass1_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_WashU.bed
	FitHiC_Pass1_PeakCCDistr_Text=$GenFitHiCDir/$PREFIX.anchorPeakCCDistr.bed
	FitHiC_Pass1_Filt_MergedIntfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED.bed
	FitHiC_Pass1_Filt_MergedInt_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED_WashU.bed

	# Modeling the statistical significance by FitHiC - main function
	if [ ! -f $FitHiC_Pass1_outfile ]; then
		Rscript ../src/Interaction.r --InpFile $CurrIntFileSortDist --OutFile $FitHiC_Pass1_outfile --Norm $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter
	fi

	echo 'Applied FitHiC'

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	# also print the header line
	if [ ! -f $FitHiC_Pass1_Filtfile ]; then
		awk -v q="$QVALUE" '{if ((NR==1) || ($NF < q && $NF > 0)) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
	fi

	echo 'Extracted significant interactions from FitHiC'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for applying FitHiC (significant interactions): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi		

	# no of significant interactions (FitHiC)
	nsigFitHiC=`cat $FitHiC_Pass1_Filtfile | wc -l`

	# Check the no of significant contacts associated with each short peak segment
	if [[ $dirname != $DirALLtoALL ]]; then
		if [ ! -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
			# At least 10 significant interactions are required (empirical threshold)
			# considering the 1st line as header
			if [[ $nsigFitHiC -gt 11 ]]; then
				# skip the 1st header line
				awk 'NR>1' $FitHiC_Pass1_Filtfile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $FitHiC_Pass1_Filt_PeakCountfile
			else
				echo 'number of significant interactions for spline distribution < 10 - skip the peak count distribution function'
			fi
		fi
	fi

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

	# the filtered interaction (with respect to the spline) file is used to create a session file
	# for applying in WashU epigenome browser
	# for that, a special file containing only the interacting chromosome intervals 
	# and the log of Q value is created
	if [ ! -s $FitHiC_Pass1_LogQ_file ]; then
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 2 ]]; then
			awk '{if (NR > 1) {print $1","$2","$3"\t"$4","$5","$6"\t"(-log($NF)/log(10))}}' $FitHiC_Pass1_Filtfile > $FitHiC_Pass1_LogQ_file
		else
			echo 'There is no significant interaction - so no WashU specific session file is created !!'
		fi
	fi

	# if [[ $DrawFig == 1 && $dirname != $DirALLtoALL ]]; then
	# 	res_outdir=$GenFitHiCDir'/Results'
	# 	mkdir -p $res_outdir
	# 	if [ -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
	# 		# Distribution of significant contact counts (FitHiC) for individual peaks
	# 		if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf ] || [ ! -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
	# 			Rscript ../Analysis/ContactCountDistr.r $PeakFILE $FitHiC_Pass1_Filt_PeakCountfile $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf $FitHiC_Pass1_PeakCCDistr_Text
	# 		fi
	# 	fi
	# fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for post processing FitHiC results: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	# if merging nearby interactions are enabled
	# then we merge the nearby interactions from the earlier generated significant interactions
	# and also create a washu browser generated compatible file
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
	fi

done 	# end of directory traversal loop

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir
