#!/bin/bash

#===============
# FitHiChIP - parsing the command line arguments 
# from a configuration file
# and calling different functions / modules 
# for processing the input alignments and finding the statistical significant interactions

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

usage(){
cat << EOF

Options:
   	-C  ConfigFile		Name of the configuration file storing the parameters of FitHiChIP.
EOF
}

while getopts "C:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

#======================
# default values of the parameters
#======================

# input files
INPFILE=""
InpBinIntervalFile=""
InpMatrixFile=""
PeakFILE=""

PREFIX='FitHiChIP'

# number of threads employed for interaction / contact detection
THREADS=1

# genome size parameter for MACS2
GSIZE='hs'

# size of the chromosome that is to be provided
ChrSizeFile=""

# 5 Kb resolution
BIN_SIZE=5000

QVALUE=1e-2	# 0.01

# default value of output directory is the present working directory
OutDir=`pwd`'/'

# lower distance threshold for two cis interactions
LowDistThres=20000	# 20 Kb

# upper distance threshold for two cis interactions
UppDistThres=2000000 # 2 Mb

# number of bins employed for FitHiC
NBins=200

# default value of plotting analysis figures
DrawFig=0

# option to note down the timing information
TimeProf=0

# default value of interaction type (peak to all)
InteractionType=0

#=========================
# bias correction related parameters
#=========================

# default lower cutoff of bias value
biaslowthr=0.2

# default higher cutoff of bias value
biashighthr=5

# boolean variable for pre-filtering the interactions according to the bias value
BeginBiasFilter=0

# boolean variable for probability adjustment of the interactions according to the bias value
EndBiasFilter=0

# Merging interactions which are near to each other
MergeInteraction=1

#==============================
# read the configuration file and store various parameters
#==============================

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			#echo "Content of $param is $paramval"
			if [ $param == "InpFile" ]; then
				INPFILE=$paramval
			fi
			if [ $param == "PeakFile" ]; then
				PeakFILE=$paramval
			fi
			if [ $param == "OutDir" ]; then
				if [[ ! -z $paramval ]]; then
					OutDir=$paramval
				fi
			fi
			if [ $param == "BINSIZE" ]; then
				if [[ ! -z $paramval ]]; then
					BIN_SIZE=$paramval
				fi
			fi
			if [ $param == "LowDistThr" ]; then
				if [[ ! -z $paramval ]]; then
					LowDistThres=$paramval
				fi
			fi
			if [ $param == "UppDistThr" ]; then
				if [[ ! -z $paramval ]]; then
					UppDistThres=$paramval
				fi
			fi
			if [ $param == "QVALUE" ]; then
				if [[ ! -z $paramval ]]; then
					QVALUE=$paramval
				fi
			fi
			if [ $param == "NBins" ]; then
				if [[ ! -z $paramval ]]; then
					NBins=$paramval
				fi
			fi
			if [ $param == "BeginBiasFilter" ]; then
				if [[ ! -z $paramval ]]; then
					BeginBiasFilter=$paramval
				fi
			fi
			if [ $param == "EndBiasFilter" ]; then
				if [[ ! -z $paramval ]]; then
					EndBiasFilter=$paramval
				fi
			fi
			if [ $param == "biaslowthr" ]; then
				if [[ ! -z $paramval ]]; then
					biaslowthr=$paramval
				fi
			fi			
			if [ $param == "biashighthr" ]; then
				if [[ ! -z $paramval ]]; then
					biashighthr=$paramval
				fi
			fi
			if [ $param == "MergeInt" ]; then
				if [[ ! -z $paramval ]]; then
					MergeInteraction=$paramval
				fi
			fi
			if [ $param == "PREFIX" ]; then
				if [[ ! -z $paramval ]]; then
					PREFIX=$paramval
				fi
			fi
			if [ $param == "Draw" ]; then
				if [[ ! -z $paramval ]]; then
					DrawFig=$paramval
				fi
			fi
			if [ $param == "TimeProf" ]; then
				if [[ ! -z $paramval ]]; then
					TimeProf=$paramval
				fi
			fi
			if [ $param == "IntType" ]; then
				if [[ ! -z $paramval ]]; then
					InteractionType=$paramval
				fi
			fi
			if [ $param == "Threads" ]; then
				if [[ ! -z $paramval ]]; then
					Threads=$paramval
				fi
			fi
			if [ $param == "GSIZE" ]; then
				if [[ ! -z $paramval ]]; then
					GSIZE=$paramval
				fi
			fi
			if [ $param == "Interval" ]; then
				InpBinIntervalFile=$paramval
			fi
			if [ $param == "Matrix" ]; then
				InpMatrixFile=$paramval
			fi
			if [ $param == "ChrSizeFile" ]; then
				ChrSizeFile=$paramval
			fi
			if [ $param == "HiCProBasedir" ]; then
				HiCProBasedir=$paramval
			fi
		fi
	fi
done < $ConfigFile

#===================
# verify the input parameters
#===================

if [[ -z $INPFILE ]]; then
	echo 'No input alignment file (BAM / PAIRIX / HiC-Pro valid pairs) is provided - exit !!'
	exit 1
fi

if [[ -z $PeakFILE ]]; then
	echo 'User should provide a reference peak detection file - exit !!'
	exit 1
fi

#==============================
# here check if the configuration file has relative path names as the input
# in such a case, convert the relative path names (with respect to the location of the configuration file itself)
# in the absolute file

# directory of the configuration file
ConfigFileDir=$(dirname "${ConfigFile}")

# first go to the configuration file directory
cd $ConfigFileDir

if [[ ! -z $INPFILE ]]; then
	if [[ "${INPFILE:0:1}" != / && "${INPFILE:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		INPFILE="$(cd $(dirname $INPFILE); pwd)/$(basename $INPFILE)"
		echo 'Absolute converted path: INPFILE: '$INPFILE
	fi
fi

if [[ ! -z $PeakFILE ]]; then
	if [[ "${PeakFILE:0:1}" != / && "${PeakFILE:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		PeakFILE="$(cd $(dirname $PeakFILE); pwd)/$(basename $PeakFILE)"
		echo 'Absolute converted path: PeakFILE: '$PeakFILE
	fi
fi

if [[ ! -z $OutDir ]]; then
	if [[ "${OutDir:0:1}" != / && "${OutDir:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		OutDir="$(cd $(dirname $OutDir); pwd)/$(basename $OutDir)"
		echo 'Absolute converted path: OutDir: '$OutDir
	fi
fi

# input file format - 1 means bam file, 2 means pairix file, 3 means HiC-pro generated valid pairs file
# BAM file has an extension .bam
# PAIRIX file has an extension of .gz, and a file ".gz.px2" also exists within the same directory
extension="${INPFILE##*.}"
if [[ "$extension" == "bam" ]]; then
	inpfilefmt=1
elif [[ "$extension" == "gz" ]]; then
	if [[ -f $INPFILE'.px2' ]]; then
		inpfilefmt=2
	else 
		inpfilefmt=3
	fi
else
	inpfilefmt=3
fi

echo 'Input file format: '$inpfilefmt

#============================================
# following files are applicable if HiC-pro generated 
# valid pairs file is provided as input
#============================================
if [ $inpfilefmt == 3 ]; then

	if [[ ! -z $ChrSizeFile ]]; then
		if [[ "${ChrSizeFile:0:1}" != / && "${ChrSizeFile:0:2}" != ~[/a-z] ]]; then
			# relative path - convert to absolute path
			ChrSizeFile="$(cd $(dirname $ChrSizeFile); pwd)/$(basename $ChrSizeFile)"
			echo 'Absolute converted path: ChrSizeFile: '$ChrSizeFile
		fi
	else
		echo 'Input file is a HiC-pro valid pairs file - but Chromosome size file is not specified - exit !!'
		exit 1
	fi

	if [[ ! -z $HiCProBasedir ]]; then
		if [[ "${HiCProBasedir:0:1}" != / && "${HiCProBasedir:0:2}" != ~[/a-z] ]]; then
			# relative path - convert to absolute path
			HiCProBasedir="$(cd $(dirname $HiCProBasedir); pwd)/$(basename $HiCProBasedir)"
			echo 'Absolute converted path: HiCProBasedir: '$HiCProBasedir
		fi
	else
		if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
			echo 'Input file is a HiC-pro valid pairs file - Contact matrix files are not provided and the base directory of HiC-pro installation path is also not provided - exit !!'
			exit 1
		fi
	fi

	if [[ ! -z $InpBinIntervalFile ]]; then
		if [[ "${InpBinIntervalFile:0:1}" != / && "${InpBinIntervalFile:0:2}" != ~[/a-z] ]]; then
			# relative path - convert to absolute path
			InpBinIntervalFile="$(cd $(dirname $InpBinIntervalFile); pwd)/$(basename $InpBinIntervalFile)"
			echo 'Absolute converted path: InpBinIntervalFile: '$InpBinIntervalFile
		fi
	fi

	if [[ ! -z $InpMatrixFile ]]; then
		if [[ "${InpMatrixFile:0:1}" != / && "${InpMatrixFile:0:2}" != ~[/a-z] ]]; then
			# relative path - convert to absolute path
			InpMatrixFile="$(cd $(dirname $InpMatrixFile); pwd)/$(basename $InpMatrixFile)"
			echo 'Absolute converted path: InpMatrixFile: '$InpMatrixFile
		fi
	fi

	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
		# executable of matrix building is to be obtained from the HiC-pro base directory
		# provided as the input
		# MatrixBuildExec=$HiCProBasedir'/scripts/build_matrix'
		MatrixBuildExecSet=( $(find $HiCProBasedir -type f -name 'build_matrix') )
		len=${#MatrixBuildExecSet[@]}
		echo 'len: '$len
		if [[ $len == 0 ]]; then
			echo 'Did not find HiC-pro package installation and the utility for matrix generation - quit !!'
			exit 1
		fi
		idx=`expr $len - 1`
		echo 'idx: '$idx
		MatrixBuildExec=${MatrixBuildExecSet[idx]}
		echo -e '\n *** MatrixBuildExec: '$MatrixBuildExec
	fi

fi

# from the configuration file directory, 
# revert to the old directory
cd -

#==============================
# important - sourya
# first change the current working directory to the directory containing this script
# it is useful when the script is invoked from a separate directory
#==============================
currworkdir=`pwd`
echo 'currworkdir: '$currworkdir
currscriptdir=`dirname $0`
echo 'currscriptdir: '$currscriptdir
cd $currscriptdir
echo 'pwd before calling executables: '`pwd`

#========================
# if the file format is 1 or 2 then call the BAM / PAIRIX processing routine
# else called the HiC-pro valid pairs processing routine
#=========================

if [[ $inpfilefmt == 3 ]]; then
	currexec=`pwd`'/bin/FitHiChIP_HiCPro.sh'
	cmd=$currexec' -I '$INPFILE
	if [[ ! -z $InpBinIntervalFile ]]; then
		cmd=$cmd' -T '$InpBinIntervalFile
	fi
	if [[ ! -z $InpMatrixFile ]]; then
		cmd=$cmd' -M '$InpMatrixFile
	fi
	cmd=$cmd' -P '$PeakFILE' -o '$OutDir' -c '$ChrSizeFile' -B '$BIN_SIZE' -L '$LowDistThres' -U '$UppDistThres' -H '$MatrixBuildExec' -n '$PREFIX' -l '$biaslowthr' -h '$biashighthr' -b '$BeginBiasFilter' -e '$EndBiasFilter' -m '$MergeInteraction' -q '$QVALUE' -v '$TimeProf' -D '$DrawFig' -N '$NBins
	# execute the command
	echo 'Invoke main executable - command: '$cmd
	eval $cmd 
else
	currexec=`pwd`'/bin/GenFitHiC_main_PAIRIX.sh'
	cmd=$currexec' -I '$INPFILE' -M '$InteractionType' -P '$PeakFILE' -o '$OutDir' -n '$PREFIX' -B '$BIN_SIZE' -L '$LowDistThres' -U '$UppDistThres' -l '$biaslowthr' -h '$biashighthr' -b '$BeginBiasFilter' -e '$EndBiasFilter' -m '$MergeInteraction' -t '$Threads' -g '$GSIZE' -q '$QVALUE' -v '$TimeProf' -D '$DrawFig' -N '$NBins
	# execute the command
	echo 'Invoke main executable - command: '$cmd
	eval $cmd 
fi

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir
