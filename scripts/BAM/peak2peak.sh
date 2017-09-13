#!/bin/bash

#===============
# sample script for executing FitHiChIP
# given a BAM alignment

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# should be a sorted and indexed bam file
AlignFile='../../TestData/BAM/inp.bam'
#PeakFile='../../TestData/BAM/inp.Peaks'
OutDir='./Results_FitHiChIP_BAM/'
codeexec='../../FitHiChIP.sh'
PREFIX='sample_input_bam'


# peak to peak interaction (-M 1)
# equal occupancy spline + binomial model (-f 1 and -B 1)
# distance thresholds: 10 Kb and 3 Mb
# no of threads = 8 (-t)
#$codeexec -I $AlignFile -P $PeakFile -o $OutDir -n $PREFIX -t 8 -v 1 -L 10000 -U 3000000 -f 1 -B 1 -M 1 -D 1
$codeexec -I $AlignFile -o $OutDir -n $PREFIX -t 8 -v 1 -L 10000 -U 3000000 -f 1 -B 1 -M 1 -D 1


