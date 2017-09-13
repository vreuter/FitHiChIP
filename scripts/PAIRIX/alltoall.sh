#!/bin/bash

#===============
# sample script for executing FitHiChIP
# given a BAM alignment

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

# should be a sorted and indexed pairix file
AlignFile='../../TestData/PAIRIX/inp.bsorted.pairs.gz'
OutDir='./Results_FitHiChIP_PAIRIX/'
codeexec='../../FitHiChIP.sh'
PREFIX='sample_input_pairix'

# all to all interaction (-M 3)
# equal occupancy spline + binomial model (-f 1 and -B 1)
# distance thresholds: 10 Kb and 3 Mb
# no of threads = 8 (-t)
$codeexec -I $AlignFile -o $OutDir -n $PREFIX -t 8 -v 1 -L 10000 -U 3000000 -f 1 -B 1 -M 3 -D 1



