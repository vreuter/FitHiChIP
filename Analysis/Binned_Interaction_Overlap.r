#!/usr/bin/env Rscript

#===========================================================
# R script for checking the mutual exclusive interactions among two different set of interactions 
# (HiC, HiChip, PLAC seq or similar interaction data between genomic fragments)
# Plots such mutual exclusive and the common interactions of both models 
# with respect to fixed size bins

# Note: Interaction files are filtered with respect to a common metric 
# (we employed Q-value based filtering in our application)

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript Binned_Interaction_Overlap.r $intfile1 $intfile2 $OutFile $name1 $name2 $col_no
# parameters: 
# 1) intfile1: interaction file generated from the first technique
# 2) intfile2: interaction file generated from the second technique
# 3) OutFile: Output file plotting the results of interaction overlap
# 4) Column_No: Column index which stores the contact count information for these interaction matrices 
# Optional parameters: 
# 5) name1: name of the first data series
# 6) name2: name of the second data series

# an example command:
# Rscript Binned_Interaction_Overlap.r NTKOCTCF_hg38.anchors.interactions_BinomDistr_FILTER.bed NTKOCTCF_hg38.anchors.interactions_SplinePass1_FILTER.bed BinomDistr_SplinePass1_FILTER_DistCC.pdf 7 'BinomDistr' 'SplineEqOcc'
#===========================================================

# package to be loaded
suppressMessages(library(GenomicRanges))

args <- commandArgs(TRUE)

# input files containing interaction matrices for two different methods
InteractionFile1 <- args[1]
InteractionFile2 <- args[2]

# output file which will store the results
OutPlotFile <- args[3]
OutPlotDir <- dirname(OutPlotFile)

# column no having the absolute contact count
ContactCol <- as.integer(args[4])

# name of the plot / graph 
Int_1_Name <- args[5]
Int_2_Name <- args[6]

colorvec <- c("red", "blue", "green", "yellow")

# load the two interactions (no header line)
Interaction1 <- read.table(InteractionFile1, header=FALSE)
Interaction2 <- read.table(InteractionFile2, header=FALSE)

#===========================
# comparison of the first three fields of two interactions
# Note: boundary locations are excluded before comparing the matching between interval ranges
# by the offset +1 and -1 in the respective conditions
#===========================
# Output: matrix of 2 columns
# 1st column: line no in first data (inp)
# second column: line no in second data (out)
ov1 <- as.data.frame(findOverlaps(GRanges(Interaction1[,1], IRanges(Interaction1[,2]+1, Interaction1[,3]-1)),GRanges(Interaction2[,1], IRanges(Interaction2[,2]+1, Interaction2[,3]-1))))

#===========================
# comparison of the next three fields between these two set of interactions
#===========================
# Output: matrix of 2 columns
# 1st column: line no in first data (inp)
# second column: line no in second data (out)
ov2 <- as.data.frame(findOverlaps(GRanges(Interaction1[,4], IRanges(Interaction1[,5]+1, Interaction1[,6]-1)),GRanges(Interaction2[,4], IRanges(Interaction2[,5]+1, Interaction2[,6]-1))))

#===========================
# Intersection of ov1 and ov2
# Common pairs of intervals
#===========================
# Output: matrix of 2 columns
# 1st column: line no in first data (inp)
# second column: line no in second data (out)
overlap_uniq <- ov1[unique(which(paste(ov1[,1], ov1[,2], sep=".") %in% paste(ov2[,1], ov2[,2], sep="."))),]

if (0) {
	# write the overlap in a separate text file
	ov1.file <- paste0(OutPlotDir, '/Ov1.txt')
	write.table(ov1, ov1.file, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

	# write the overlap in a separate text file
	ov2.file <- paste0(OutPlotDir, '/Ov2.txt')
	write.table(ov2, ov2.file, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

	# write the overlap in a separate text file
	overlap.uniq.file <- paste0(OutPlotDir, '/UniqOverlapInfo.txt')
	write.table(overlap_uniq, overlap.uniq.file, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
}

# load only the genomic distance and contact count information from both interaction data
Int1_Dist_CC <- cbind(abs(Interaction1$V2-Interaction1$V5), Interaction1[,ContactCol])
Int2_Dist_CC <- cbind(abs(Interaction2$V2-Interaction2$V5), Interaction2[,ContactCol])

# For the first interaction set (the genomic distance and contact count), create two structures
# one contains the common interaction pattern (with respect to the other interaction set)
# and the other contains the unique interactions 
# The common interaction lines are stored in the first column of the overlap structure
y <- seq(1:nrow(Int1_Dist_CC)) 
Int1_Dist_CC_Common <- subset(Int1_Dist_CC, (y %in% overlap_uniq[,1]))
Int1_Dist_CC_Uniq <- subset(Int1_Dist_CC, !(y %in% overlap_uniq[,1]))

# filter the NA's
Int1_Dist_CC_Common <- na.omit(Int1_Dist_CC_Common)
Int1_Dist_CC_Uniq <- na.omit(Int1_Dist_CC_Uniq)

# Similarly compute for the second interaction
# The common interaction lines are stored in the second column of the overlap structure
y <- seq(1:nrow(Int2_Dist_CC)) 
Int2_Dist_CC_Common <- subset(Int2_Dist_CC, (y %in% overlap_uniq[,2]))
Int2_Dist_CC_Uniq <- subset(Int2_Dist_CC, !(y %in% overlap_uniq[,2]))

# filter the NA's
Int2_Dist_CC_Common <- na.omit(Int2_Dist_CC_Common)
Int2_Dist_CC_Uniq <- na.omit(Int2_Dist_CC_Uniq)

if ((nrow(Int1_Dist_CC_Common) == 0) || (nrow(Int2_Dist_CC_Common) == 0)) {
	cat(sprintf("\n There is no common interaction between the binomial and the fithic - exit \n"))
	return
}

# sort all of these structures with respect to ascending order of genomic distance
Int1_Dist_CC_Common_Sort <- Int1_Dist_CC_Common[ order(Int1_Dist_CC_Common[,1]),]
Int2_Dist_CC_Common_Sort <- Int2_Dist_CC_Common[ order(Int2_Dist_CC_Common[,1]),]

# cat(sprintf("\n nrow Int2_Dist_CC_Uniq: %s  ncol Int2_Dist_CC_Uniq: %s ", nrow(Int2_Dist_CC_Uniq), ncol(Int2_Dist_CC_Uniq)))

# cat(sprintf("\n nrow Int1_Dist_CC_Uniq: %s  ncol Int1_Dist_CC_Uniq: %s ", nrow(Int1_Dist_CC_Uniq), ncol(Int1_Dist_CC_Uniq)))

# cat(sprintf("\n Int1_Dist_CC_Uniq: %s ", Int1_Dist_CC_Uniq))
# cat(sprintf("\n Int2_Dist_CC_Uniq: %s ", Int2_Dist_CC_Uniq))

# cat(sprintf("\n START GHOST \n"))

if (nrow(Int1_Dist_CC_Uniq) > 0) {
	if (nrow(Int1_Dist_CC_Uniq) == 1) {
		Int1_Dist_CC_Uniq_Sort <- Int1_Dist_CC_Uniq
	} else {
		Int1_Dist_CC_Uniq_Sort <- Int1_Dist_CC_Uniq[ order(Int1_Dist_CC_Uniq[,1]),]
	}
	
	# cat(sprintf("\n GHOST \n"))
	# cat(sprintf("\n nrow Int1_Dist_CC_Uniq_Sort: %s  ncol Int1_Dist_CC_Uniq_Sort: %s ", nrow(Int1_Dist_CC_Uniq_Sort), ncol(Int1_Dist_CC_Uniq_Sort)))
	# cat(sprintf("\n Int1_Dist_CC_Uniq_Sort: %s ", Int1_Dist_CC_Uniq_Sort))
	# cat(sprintf("\n GHOST \n"))
} 

if (nrow(Int2_Dist_CC_Uniq) > 0) {
	if (nrow(Int2_Dist_CC_Uniq) == 0) {
		Int2_Dist_CC_Uniq_Sort <- Int2_Dist_CC_Uniq
	} else {
		Int2_Dist_CC_Uniq_Sort <- Int2_Dist_CC_Uniq[ order(Int2_Dist_CC_Uniq[,1]),]
	}
	
	# cat(sprintf("\n nrow Int2_Dist_CC_Uniq_Sort: %s  ncol Int2_Dist_CC_Uniq_Sort: %s ", nrow(Int2_Dist_CC_Uniq_Sort), ncol(Int2_Dist_CC_Uniq_Sort)))
	# cat(sprintf("\n Int2_Dist_CC_Uniq_Sort: %s ", Int2_Dist_CC_Uniq_Sort))
}

#===========================
# the following step adds only those vectors (and legend entries)
# which are non empty
#===========================

# create the legend vectors
if (length(args) == 5) {
	legendvec <- c(paste0(Int_1_Name, "_Common"), paste0(Int_2_Name, "_Common"))
} else {
	legendvec <- c("Int1_Common", "Int2_Common")
}

# plot the lines
pdf(OutPlotFile, width=14, height=10)
plot(Int1_Dist_CC_Common_Sort[,1], Int1_Dist_CC_Common_Sort[,2], cex=0.1, lwd=0.3, col=colorvec[1], xlab="Genomic distance", ylab="Contact count")
lines(Int2_Dist_CC_Common_Sort[,1], Int2_Dist_CC_Common_Sort[,2], col=colorvec[2], cex=0.3, lwd=0.3, type="p", pch=1) 

# if the unique interaction vectors are non-empty
# they are plotted

if (nrow(Int1_Dist_CC_Uniq) > 0) {
	if (length(args) == 5) {
		legendvec <- c(legendvec, paste0(Int_1_Name, "_Uniq"))
	} else {
		legendvec <- c(legendvec, "Int1_Uniq")
	}
	lines(Int1_Dist_CC_Uniq_Sort[,1], Int1_Dist_CC_Uniq_Sort[,2], col=colorvec[3], cex=0.3, lwd=0.3, type="p", pch=1) 
}

if (nrow(Int2_Dist_CC_Uniq) > 0) {
	if (length(args) == 5) {
		legendvec <- c(legendvec, paste0(Int_2_Name, "_Uniq"))
	} else {
		legendvec <- c(legendvec, "Int2_Uniq")
	}
	lines(Int2_Dist_CC_Uniq_Sort[,1], Int2_Dist_CC_Uniq_Sort[,2], col=colorvec[4], cex=0.3, lwd=0.3, type="p", pch=1)
}

title("Genomic distance vs contact count for two different distribution models - unique and common interactions")
legend("topright",legend=legendvec, col=colorvec, lty=1, lwd=2, cex=0.8)
dev.off()


