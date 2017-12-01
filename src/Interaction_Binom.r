#!/usr/bin/env Rscript

#===========================================================
# R script for analyzing the statistical significance of interactions
# this is an implementation of the paper Duan et. al.
# which models the interactions as a binomial distribution

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript interaction.r $interaction.file $out_interaction_file $nbins $Contact_Col $timing_file
#===========================================================

# process input arguments
args <- commandArgs(TRUE)

# 1st argument is the sorted interaction file (with respect to the genomic distance between the fragments)
interaction.file <- args[1]
inpdir <- dirname(interaction.file)

# second argument is the output file containing the interactions with their relative P and Q values
# computed using the binomial distribution modeling
outfile <- args[2]

# number of bins employed for binomial distribution
if (length(args) > 2) {
	nbins <- as.integer(args[3])
} else {
	nbins <- 200
}

# column no having the absolute contact count
if (length(args) > 3) {
	Contact_Col <- as.integer(args[4])
} else {
	Contact_Col <- 7
}

# we check if time profiling is enabled or not
if (length(args) > 4) {
	timeprof <- 1
	timefile <- args[5]
} else {
	timeprof <- 0
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# load the interaction data matrix
# Note: in the new interaction data (with normalized contact count)
# the data has a header information
interaction.data <- read.table(interaction.file, header=T)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time to load the interaction data file in the R structure: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

# return if the number of interaction is 0
if (nrow(interaction.data) == 0) {
	cat(sprintf("The interaction file is empty - no binomial distribution model can be computed - return !!"))
	return
}

# no of interactions pairs
numPairs <- length(interaction.data[,2])

# column storing the prior probability employed for individual interactions
Prior_Prob_CC <- c()

# a separate column of the same column length as of 'interaction.data' is used
# it will store the probability of the observed contact count for a single locus pair
# binomial distribution is employed (Duan et. al. 2010)
Binom_Prob_CC <- c()

# a separate column of the same column length as of 'interaction.data' is used
# it will store the P-value of the observed contact count for a single locus pair
# binomial distribution is employed (Duan et. al. 2010)
Binom_P_Val_CC <- c()

#=====================================================
# divide the genomic distance into b quantiles (bins) of equal occupancy
# currently we assume b = 200 bins
# Note: each bin would have equal no of entries
# but the corresponding bin interval will be variable


# error condition - sourya
# if the number of interactions is less than the argument 'nbins'
# then re-adjust the values of nbins and then assign the 
# no of entries in each of the bins
if (numPairs < nbins) {
	nbins <- numPairs
	nentry <- 1
} else {
	# no of entries (occupancies) for each bin (no of pairs)
	nentry <- floor(numPairs / nbins)	
}

# cat(sprintf("\n numPairs: %s  nbins: %s  nentry: %s \n ", numPairs, nbins, nentry))

#======================================================
# compute various probability distributionns and 
# corresponding p values
# for individual pairs of loci and their contact counts
#--------------------------------------------

if (timeprof == 1) {
	starttime <- Sys.time()
}

# Paper: Duan et. al. 2010
# binomial distribution and the p-value for it
# probability of a contact between a pair of loci, for this particular bin

for (i in (1:nbins)) {
	# start and end index of this bin (row no of the loaded data table)
	si <- (i - 1) * nentry + 1
	if (i < nbins) {
		ei <- si + nentry - 1		
	} else {
		ei <- numPairs	# last read
	}
	# no of distinct pair of loci for this bin
	no_distinct_loci <- (ei-si+1)
	# no of contacts for this bin
  	NumContact <- sum(interaction.data[si:ei, Contact_Col])
	p <- 1.0 / no_distinct_loci		
	curr_dbinom <- dbinom(interaction.data[si:ei, Contact_Col], size=NumContact, prob=p)
	curr_pbinom <- pbinom(interaction.data[si:ei, Contact_Col], size=NumContact, prob=p, lower.tail=FALSE)
	# binomial distribution based probability
 	Binom_Prob_CC <- c(Binom_Prob_CC, curr_dbinom)
 	# P value based on binomial distribution
 	Binom_P_Val_CC <- c(Binom_P_Val_CC, (curr_dbinom + curr_pbinom))
 	# prior probability (p) values employed for modeling the binomial distribution
 	Prior_Prob_CC[si:ei] <- p

 	# debug - sourya
 	# cat(sprintf("\n Prob compute - i: %s si: %s ei: %s  p: %s \n ", i, si, ei, p))
}

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Duan paper : binomial distribution: time to compute p value: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

# debug
# cat(sprintf("\n P val is computed - Before computing the Q val \n "))

if (timeprof == 1) {
	starttime <- Sys.time()
}

# from the generated P values, obtain the Q value using BH correction
Binom_QVal <- p.adjust(Binom_P_Val_CC, method = "BH")

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Duan paper : binomial distribution: time to compute Q value (from P value): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# debug
# cat(sprintf("\n After computing the Q val \n "))

# accumulate all results - also add header information
FinalData <- cbind(interaction.data, Prior_Prob_CC, Binom_Prob_CC, Binom_P_Val_CC, Binom_QVal)
colnames(FinalData) <- c(colnames(interaction.data), "p", "dbinom", "P-Value", "Q-Value")

# append the Spline distribution probability and corresponding P value as separate columns
# and write in a separate text file
temp_outfile <- paste0(inpdir, '/', 'temp_out.bed')
write.table(FinalData, temp_outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 

# now sort the file contents and write that in the final specified output file
system(paste('sort -k1,1 -k2,2n -k5,5n ', temp_outfile, '>', outfile))
# delete the temporary output file
system(paste('rm', temp_outfile))

# # sort the data according to the chromosome name, start positions of both the interacting segments
# FinalData <- FinalData[order( FinalData[,1], FinalData[,2], FinalData[,5] ), ]

# # write the data with the header information
# write.table(FinalData, outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Duan paper : binomial distribution: time to write interaction data (with P, Q values) to a file: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}


