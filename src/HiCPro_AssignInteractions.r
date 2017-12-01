#!/usr/bin/env Rscript

#===========================================================
# R script for assigning different types of interactions
# from a given collection of ALL to ALL interactions from a HiC-pro pipeline output
# using validpairs.txt and associated contact matrix
# the idea is to assign the contacts from peak to peak, peak to non peak, and peak to all interactions

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

# library to be loaded
suppressMessages(library(GenomicRanges))

args <- commandArgs(TRUE)

IntFileALLtoALL <- args[1]
PeakFILE <- args[2]
IntFilePeaktoPeak <- args[3]
IntFilePeaktoNonPeak <- args[4]
IntFilePeaktoALL <- args[5]

# temp argument
OutDir <- args[6]

# read the cis interactions (among all pairs of bins, subject to the distance thresholds)
InteractionAllData <- read.table(IntFileALLtoALL, header=T)

# number of fields in the peak file
nfieldPeak <- max(count.fields(PeakFILE, sep = "\t"))

# read the peak information (only the first three columns)
# first three columns denote chromosome, start, and end values
PeakData <- read.table(PeakFILE, colClasses = c("character", rep("integer", 2), rep("NULL", (nfieldPeak-3))), header=F)

#===========================
# checking overlap between the peak regions (in the peak data) and the first three columns of the "InteractionAllData"
# Note: boundary locations are excluded before comparing this overlap
# by the offset +1 and -1 in the respective conditions

# Output: matrix of 2 columns
# 1st column: row no in peak data (with respect to the peak file)
# second column: row no in all interaction data (with respect to the first three columns)
#===========================
ov1 <- as.data.frame(findOverlaps(GRanges(PeakData[,1], IRanges(PeakData[,2]+1, PeakData[,3]-1)),GRanges(InteractionAllData[,1], IRanges(InteractionAllData[,2]+1, InteractionAllData[,3]-1))))

# debug
# write.table(ov1, paste0(OutDir,'/ov1.txt'), row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

#===========================
# checking overlap between the peak regions (in the peak data) and the columns 4-6 
# (second interacting segment) of the "InteractionAllData"
# Note: boundary locations are excluded before comparing this overlap
# by the offset +1 and -1 in the respective conditions

# Output: matrix of 2 columns
# 1st column: row no in peak data (with respect to the peak file)
# second column: row no in all interaction data (with respect to the columns 4-6)
#===========================
ov2 <- as.data.frame(findOverlaps(GRanges(PeakData[,1], IRanges(PeakData[,2]+1, PeakData[,3]-1)),GRanges(InteractionAllData[,4], IRanges(InteractionAllData[,5]+1, InteractionAllData[,6]-1))))

# debug
# write.table(ov2, paste0(OutDir,'/ov2.txt'), row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

#===========================
# Peak to Peak Interactions involve intersection between the second columns of ov1 and ov2
# Reason: suppose, ov1 has the entry "x z" 
# where x is the row no in the peak file and z is the row no of the interaction file (with respect to first chr interval)
# now suppose ov2 has the entry "y z"
# where y is the peak file entry which matches with the second chr interval of the line no z in the interaction file
# so here z is the line no (interaction file) having peak to peak interaction
#===========================
Peak2PeakLineVec <- intersect(ov1[,2], ov2[,2])

# debug
# write.table(Peak2PeakLineVec, paste0(OutDir,'/Peak2PeakLineMatrix.txt'), row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

# write the peak to peak interactions
write.table(InteractionAllData[Peak2PeakLineVec, ], IntFilePeaktoPeak, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

#===========================
# Peak to Non Peak Interactions involve non-intersection between the second columns of ov1 and ov2
# Reason: suppose, ov1 has the entry "x z" 
# where x is the row no in the peak file and z is the row no of the interaction file (with respect to first chr interval)
# now suppose ov2 has the entry "y z"
# where y is the peak file entry which matches with the second chr interval of the line no z in the interaction file
# so the line no z should be excluded
# specifically we shall include those lines z (from the interaction file)
# such that "x z" exists but "y z" does not exist
#===========================
Peak2NonPeakLineVec <- setdiff(ov1[,2], ov2[,2])

# debug
# write.table(Peak2NonPeakLineVec, paste0(OutDir,'/Peak2NonPeakLineMatrix.txt'), row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

# write the peak to non peak interactions
write.table(InteractionAllData[Peak2NonPeakLineVec, ], IntFilePeaktoNonPeak, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

#===========================
# Peak to ALL Interactions means that 
# from the above comments we should include all line numbers z
# such that "x z" exists
# so just check the ov1 matrix
#===========================
Peak2All.rows <- unique(ov1[,2])

# write the peak to non peak interactions
write.table(InteractionAllData[Peak2All.rows, ], IntFilePeaktoALL, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)






