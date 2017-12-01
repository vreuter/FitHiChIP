#!/usr/bin/env Rscript

#===========================================================
# R script for finding the significant contact count distribution for individual peaks

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript ContactCountDistr.r $PeakFile $PeakContactFile $PlotFile $OutText
# PeakFile: Input peak file (3 columns)
# PeakContactFile: Bin specific significant contact count information (4 columns)
# PlotFile: Output file storing the plots
# OutText: Output text file
#===========================================================
suppressMessages(library(GenomicRanges))

args <- commandArgs(TRUE)

PeakFile <- args[1]
PeakContactFile <- args[2]
PlotFile <- args[3]
OutText <- args[4]

# peaks in the input peak file (PD = PeakData)
PD <- read.table(PeakFile, header=FALSE)

# peaks from the interaction file (PI = peaks from interactions)
PI <- read.table(PeakContactFile, header=TRUE)

#===========================
# comparison of the reference peaks with the first three fields of bin specific contact count file
#===========================
# Output: matrix of 2 columns
# 1st column: line no in first data (inp)
# second column: line no in second data (out)
Overlap_Peaks <- as.data.frame(findOverlaps(GRanges(PD[,1], IRanges(PD[,2], PD[,3])),GRanges(PI[,1], IRanges(PI[,2], PI[,3]))))

# This vector stores the significant contact count (filtered according to Q value) 
# for individual peaks
count <- c()

for (i in (1:nrow(PD))) {
	# check if this peak is associated with any interactions
	# this is done by finding the row in "Overlap_Peaks" whose first element is i
	idx <- which(Overlap_Peaks[,1] == i)
	if (length(idx) > 0) {
		# this vector stores the row numbers of the overlapping interactions with the given peak segment
		IntRows <- Overlap_Peaks[idx,2]
		# the contact count information is in the column number provided as an input
		count[i] <- sum(PI[IntRows, 4])
	} else {
		count[i] <- 0
	}
}

# dump the original peaks vs significant contact count in a text file
write.table(cbind(PD[,1:3], count), OutText, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

# first sort the count vector and then draw the frequency distribution
pdf(PlotFile, width=14, height=10)
plot(sort(count, decreasing=TRUE), type = "o", cex=0.5, col="red", xlab="Peak counter", ylab="Contact count (filtered)")
title(sub('\\.pdf$', '', basename(PlotFile)))
dev.off()


