#!/usr/bin/env Rscript

#===========================================================
# R script for plotting the distribution of coverage values for different segments
# peak and non peak segments are plotted separately

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript PlotCoverageDistr.r $inpfile
#===========================================================

# function for plotting distribution
PlotDistr <- function(inpvec, plotfile, titlestr, xlabelstr, ylabelstr) {
	pdf(plotfile, width=14, height=10)
	# comment - sourya
	# plot(inpvec, main=titlestr, xlab=xlabelstr, ylab="Bias (ratio w.r.t mean)", col="red")
	# plot a histogram where X axis denotes the bias and Y axis denotes the corresponding frequency
	# step size of the histogram is 0.1 (bias value)
	maxelem <- max(inpvec)
	stepsize <- 0.1
	no_of_breaks <- as.integer(maxelem / stepsize)
	hist(inpvec, main=titlestr, xlab=xlabelstr, ylab=ylabelstr, border="red", col="blue", breaks=no_of_breaks)
	dev.off()
}

# process the arguments
args <- commandArgs(TRUE)

# coverage file
inpfile <- args[1]
inpdir <- dirname(inpfile)

# output directory for plotting
OutDir <- args[2]
system(paste('mkdir -p', OutDir))

# read the genome coverage features
# Note: the coverage file has header information
# columns: chromosome interval, read depth, is_peak
CoverageFeat <- read.table(inpfile, header=T)	
# colnames(CoverageFeat) <- c("chr1","s1","e1","depth","isPeak")

# note the indices having peak information and non zero coverage
PeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==1))
# note the indices of non peak and non zero coverage
nonPeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==0))

# process the bias associated with peak segments
if (length(PeakIDXNonZeroCoverage) > 0) {
	# coverage values (non zero) of the peak segments
	CoveragePeakVec <- CoverageFeat[PeakIDXNonZeroCoverage, 4]
	# divide bias values by the mean
	meanCoverage <- mean(CoveragePeakVec)
	CoveragePeakVec <- CoveragePeakVec / meanCoverage
	OutPlotFile <- paste0(OutDir, '/PeakCoverageDistr.pdf') 
	PlotDistr(CoveragePeakVec, OutPlotFile, "Bias distribution for peaks", "Bias (ratio w.r.t mean)", "Peak frequency")
}

# process the bias associated with non-peak segments
if (length(nonPeakIDXNonZeroCoverage) > 0) {
	# coverage values (non zero) of the non peak segments
	CoverageNonPeakVec <- CoverageFeat[nonPeakIDXNonZeroCoverage, 4]
	# divide bias values by the mean
	meanCoverage <- mean(CoverageNonPeakVec)
	CoverageNonPeakVec <- CoverageNonPeakVec / meanCoverage
	OutPlotFile <- paste0(OutDir, '/NonPeakCoverageDistr.pdf') 
	PlotDistr(CoverageNonPeakVec, OutPlotFile, "Bias distribution for non-peaks", "Bias (ratio w.r.t mean)", "Non-Peak frequency")
}

