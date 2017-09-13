#!/usr/bin/env Rscript

#===========================================================
# R script for comparing the significant contact count distribution with respect to a given set of peaks (anchors)
# between two different statistical significance models 
# 1) binomial distribution (Duan et al 2010), 2) FitHiC (Ferhat AY et al 2014)

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript PeakBasisContactDistr.r $BinomDistrPeakCCFile $FitHiCPeakCCFile $OutFile
#===========================================================

args <- commandArgs(TRUE)

BinomDistrPeakCCFile <- args[1]
FitHiCPeakCCFile <- args[2]
OutFile <- args[3]

BinomDistrPeakCC <- read.table(BinomDistrPeakCCFile, header=F)
colnames(BinomDistrPeakCC) <- c('chr', 'start', 'end', 'Binom_Contact')
FitHiCPeakCC <- read.table(FitHiCPeakCCFile, header=F)
colnames(FitHiCPeakCC) <- c('chr', 'start', 'end', 'FitHiC_Contact')

OutData <- merge(x=BinomDistrPeakCC, y=FitHiCPeakCC, by.x=colnames(BinomDistrPeakCC)[1:3], by.y=colnames(FitHiCPeakCC)[1:3])
colnames(OutData) <- c('chr', 'start', 'end', 'Binom_Contact', 'FitHiC_Contact')

write.table(OutData, OutFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

