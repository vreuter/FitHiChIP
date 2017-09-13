#!/usr/bin/env Rscript

#===========================================================
# R script for appending the interaction file (with contact counts)
# with the features like read depth (for both intervals), peak information, and the normalized contact count

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

# disabling the scientific notations
options(scipen = 999)

library("optparse")

# Sourya - Note - the string after '--'' and the metavar field should be identical
option_list = list(
  	make_option(c("-I", "--IntFile"), type="character", default=NULL, help="File having interaction among segments", metavar="IntFile"),
	make_option(c("-O", "--OutFile"), type="character", default=NULL, help="Output file name storing interactions + features", metavar="OutFile"),
	make_option(c("-E", "--ExtraFeatureFile"), type="character", default=NULL, help="Read Depth containing file", metavar="ExtraFeatureFile")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$IntFile)) {
	print_help(opt_parser)
	stop("Interaction file (pairs of genomic intervals with contact count) is not provided - check the option -I \n", call.=FALSE)
} else {
	cat(sprintf("Input interaction file: %s \n", opt$IntFile))
}

if (is.null(opt$OutFile)) {
	print_help(opt_parser)
	stop("Output file (for storing the interactions + features) is not provided - check the option -O \n", call.=FALSE)
} else {
	cat(sprintf("Output file: %s \n", opt$OutFile))
}

if (is.null(opt$ExtraFeatureFile)) {
	print_help(opt_parser)
	stop("Additional Feature file of individual genomic intervals is not provided - check the option -E \n", call.=FALSE)
} else {
	cat(sprintf("Additional Feature file of individual genomic intervals: %s \n", opt$ExtraFeatureFile))
}

# output directory
OutDir <- dirname(opt$OutFile)
cat(sprintf("\n Output directory: %s \n", OutDir))

# load the interaction matrix (pairs of intervals and their contacts)
# contains a pair of chromosome intervals (may be peaks also) and the contact count
Interaction_Mat <- read.table(opt$IntFile, header=F)
colnames(Interaction_Mat) <- c("chr1","s1","e1","chr2","s2","e2","cc")

# File containing additional features (such as read depth) of different genomic intervals
# Note: this file has header line
AllFeatures <- read.table(opt$ExtraFeatureFile, header=T)	
# columns: chromosome interval, read depth, is_peak
colnames(AllFeatures) <- c("chr1","s1","e1","depth","isPeak")

#================================================
# add - sourya
# In the matrix "AllFeatures", read depth for individual genomic segments are integers
# we convert these depth values by dividing with the average read depth value for all of these segments
#================================================
ReadDepthVec <- AllFeatures[,4]
ReadDepthVec <- ReadDepthVec[ReadDepthVec != 0]
AvgDepth <- mean(ReadDepthVec)
AllFeatures[,4] <- AllFeatures[,4] / AvgDepth

#================================================
# now merge the interactions with the normalized read depth (called the bias) and isPeakInformation
#================================================

# merge with respect to the 1st three fields of either data (chromosome interval)
df1 <- merge(x=Interaction_Mat, y=AllFeatures, by.x=colnames(Interaction_Mat)[1:3], by.y=colnames(AllFeatures)[1:3])
colnames(df1) <- c(colnames(Interaction_Mat), "bias1", "isPeak1")

# merge with respect to the next three fields (2nd chromosome interval)
# in the merged output, these columns will be printed as the first three columns
Final_Intrc <- merge(x=df1, y=AllFeatures, by.x=colnames(Interaction_Mat)[4:6], by.y=colnames(AllFeatures)[1:3])

# merge operation changes the order of columns
# here the first three columns are exchanged with the columns 4, 5, 6
# so we have to retrieve the original order
Final_Intrc <- Final_Intrc[,c(4:6,1:3,7:ncol(Final_Intrc))]

# now update the column names of this data frame
# to correctly reflect the columns
# it will be printed in the final output file
colnames(Final_Intrc) <- c(colnames(df1), "bias2", "isPeak2")

# #===================================
# # perform a sample normalization of the data
# # here we divide the contact count by the product of read depth values
# #===================================
# NormCC <- Final_Intrc$cc / (Final_Intrc$d1 * Final_Intrc$d2)

# # replace the inf, NA and NAN value with 0
# NormCC[!is.finite(NormCC)] <- 0 

# # append the normalized contact count to the already generated interaction matrix
# # to create a new data frame
# Final_Intrc_with_NormCC <- data.frame(cbind(Final_Intrc, NormCC))

# # adjust the column names of this new data frame
# colnames(Final_Intrc_with_NormCC) <- c(colnames(Final_Intrc), "NormCC")

# write the complete interaction matrix and features in the specified output file
write.table(Final_Intrc, opt$OutFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 
