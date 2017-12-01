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
	make_option(c("-E", "--AllFeatFile"), type="character", default=NULL, help="File storing normalization features for genomic bins", metavar="AllFeatFile")
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

if (is.null(opt$AllFeatFile)) {
	print_help(opt_parser)
	stop("Additional Feature file of individual genomic intervals is not provided - check the option -E \n", call.=FALSE)
} else {
	cat(sprintf("Additional Feature file of individual genomic intervals: %s \n", opt$AllFeatFile))
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
AllFeatures <- read.table(opt$AllFeatFile, header=T)	
# columns: chromosome interval, read depth, is_peak
colnames(AllFeatures) <- c("chr1","s1","e1","Coverage","isPeak", "Bias")

#================================================
# merge the interactions with the normalization related features
# such that individual interacting pairs have also their normalization related features listed
#================================================
# merge with respect to the 1st three fields of either data (chromosome interval)
df1 <- merge(x=Interaction_Mat, y=AllFeatures, by.x=colnames(Interaction_Mat)[1:3], by.y=colnames(AllFeatures)[1:3])
colnames(df1) <- c(colnames(Interaction_Mat), "Coverage1", "isPeak1", "Bias1")

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
colnames(Final_Intrc) <- c(colnames(df1), "Coverage2", "isPeak2", "Bias2")

# write the complete interaction matrix and features in the specified output file
write.table(Final_Intrc, opt$OutFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 
