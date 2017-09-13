#!/usr/bin/env Rscript

#===========================================================
# R script for plotting the distribution of coverage values for different segments

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript PlotCoverageDistr.r $inpfile
#===========================================================

# function for plotting distribution
PlotDistr <- function(inpvec, plotfile, titlestr) {
	pdf(plotfile, width=14, height=10)
	plot(inpvec, main=titlestr, xlab="Peaks", ylab="Bias (ratio w.r.t mean)", col="red")
	dev.off()
}


args <- commandArgs(TRUE)

inpfile <- args[1]
inpdir <- dirname(inpfile)

OutDir <- args[2]
system(paste('mkdir -p', OutDir))

CoverageFeat <- read.table(inpfile, header=T)	
# columns: chromosome interval, read depth, is_peak
colnames(CoverageFeat) <- c("chr1","s1","e1","depth","isPeak")

PeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==1))
nonPeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==0))

if (length(PeakIDXNonZeroCoverage) > 0) {
	# coverage values (non zero) of the peak segments
	CoveragePeakVec <- CoverageFeat[PeakIDXNonZeroCoverage, 4]
	# divide by mean
	meanCoverage <- mean(CoveragePeakVec)
	CoveragePeakVec <- CoveragePeakVec / meanCoverage

	OutPlotFile <- paste0(OutDir, '/PeakCoverageDistr.pdf') 
	PlotDistr(CoveragePeakVec, OutPlotFile, "Bias distribution for peaks")
}

if (length(nonPeakIDXNonZeroCoverage) > 0) {
	# coverage values (non zero) of the non peak segments
	CoverageNonPeakVec <- CoverageFeat[nonPeakIDXNonZeroCoverage, 4]
	# divide by mean
	meanCoverage <- mean(CoverageNonPeakVec)
	CoverageNonPeakVec <- CoverageNonPeakVec / meanCoverage

	OutPlotFile <- paste0(OutDir, '/NonPeakCoverageDistr.pdf') 
	PlotDistr(CoverageNonPeakVec, OutPlotFile, "Bias distribution for non-peaks")
}






