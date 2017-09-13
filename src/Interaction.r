#!/usr/bin/env Rscript

#===========================================================
# R script for analyzing the statistical significance of interactions
# this is an implementation of the paper Ay et. al. 2014 (FitHiC)

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript interaction.r $interaction.file $out_interaction_file $bin_type $drawfig $timing_file
# parameter description:
# 1) $interaction.file: Input interaction file (sorted by genomic distance)
# 2) $out_interaction_file: Output interaction file, with the probability, P value and Q value information
# 3) $bin_type: binary value depicting the type of binning employed
# 		0: Equal number of pairs of loci per bin
# 		1: Equal occupancy (contact counts) bins (default) 
# 4) Drawfig: binary value. If 1, plots the spline fit curve.
# 5) $timing_file: If specified, notes the timing information for individual steps of the algorithm
#===========================================================

library(splines)
library(fdrtool)
library(parallel)	# library for parallel processing

# #===========================
# # parallel computation of binomial distribution on input vector
# parallel.BinomDistr <- function(inpvec) {
# 	idx <- inpvec[1]
# 	cc <- inpvec[2]
# 	tc <- inpvec[3]
# 	currprob <- inpvec[4]

# 	curr_dbinom <- dbinom(cc, size=tc, prob=currprob)
# 	curr_pbinom <- pbinom(cc, size=tc, prob=currprob, lower.tail=FALSE)
# 	return c(idx, curr_dbinom, (curr_dbinom + curr_pbinom))
# }

# #===========================

# number of processors within the system
ncore <- detectCores()
cat(sprintf("\n Number of cores in the system: %s ", ncore))

# process input arguments
args <- commandArgs(TRUE)

# 1st argument is the sorted interaction file
interaction.file <- args[1]

# second argument is the output file containing the interactions with their relative P and Q values
# computed using the Spline distribution modeling
outfile <- args[2]

# we check the binning type
binning_type <- as.integer(args[3])

# boolean value indicating the use of bias correction
BiasCorr <- as.integer(args[4])

# number of bins
nbins <- as.integer(args[5])

# we check if plotting the spline curve is allowed or not
drawfig <- as.integer(args[6])

# column no having the absolute contact count
Contact_Col <- as.integer(args[7])

# we check if time profiling is enabled or not
if (length(args) > 7) {
	timeprof <- 1
	timefile <- args[8]
} else {
	timeprof <- 0
}

# directory containing the input interaction file
inpdir <- dirname(interaction.file)

# directory to contain the spline fitted output interaction file
OutIntDir <- dirname(outfile)

# this directory will store the spline graph according to the model
outdir <- paste0(OutIntDir,'/Results')
system(paste('mkdir -p', outdir))

# this file stores the spline fitted model
if (binning_type == 0) {
	plotfile <- paste0(outdir,'/','EqLenBin_SplinePass1.pdf')
} else {
	plotfile <- paste0(outdir,'/','EqOccBin_SplinePass1.pdf')	
}

 
if (timeprof == 1) {
	starttime <- Sys.time()
}

# load the interaction data matrix
# Note: the data has a header information
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
	return
}

# absolute genomic distance for an interaction instance
gene.dist <- abs(interaction.data[,2] - interaction.data[,5])

# no of interactions pairs
numPairs <- length(gene.dist)
cat(sprintf("\n *** Total number of interactions: %s ", numPairs))

# total number of contacts for all the interactions
TotContact <- sum(interaction.data[,Contact_Col])

# lists the initial (prior) probabilities obtained from the spline fit
Prob_Val <- c()

# a separate column of the same column length as of 'interaction.data' is used
# it will store the probability of the observed contact count for a single locus pair
# binomial distribution + spline based estimation is employed (Ferhat Ay et. al. 2014)
Spline_Binom_Prob_CC <- c()

# a separate column of the same column length as of 'interaction.data' is used
# it will store the P-value of the observed contact count for a single locus pair
# binomial distribution + spline based estimation is employed (Ferhat Ay et. al. 2014)
Spline_Binom_P_Val_CC <- c()

# for each bin, count the no of distinct locus pairs
no_distinct_loci <- c()
# for each bin, no of observed contacts
NumContact <- c()
# for each bin, stores the average contact count per locus pair
avg_contact <- c()
# prior contact probability for a specific locus pair falling within a bin
prior_contact_prob <- c()
# average interaction distance for all locus pairs falling within a bin
avg_int_dist <- c()

#=====================================================
# divide the genomic distance into b quantiles (number of bins)

# Note: for equal length bins (binning_type = 0), each bin would have equal no of entries (pairs of loci)
# on the other hand, for equal occupancy bins (binning_type = 1), each bin would have equal no of contact count
# In both cases, however, corresponding bin intervals will be variable

if (timeprof == 1) {
	starttime <- Sys.time()
}

if (binning_type == 0) {
	# no of entries (occupancies) for each bin 
	# each bin have same no of locus pairs
	nentry <- floor(numPairs / nbins)	
} else {
	# no of contacts for each bin
	# each bin have same no of contacts
	ncontactbin <- floor(TotContact / nbins)
}

#==============================================

if (binning_type == 0) {
	# equal length bin case (equal locus entries)
	for (bin_idx in (1:nbins)) {
		si <- (bin_idx - 1) * nentry + 1
		if (bin_idx < nbins) {
			ei <- si + nentry - 1		
		} else {
			ei <- numPairs	# last read
		}
		no_distinct_loci[bin_idx] <- (ei-si+1)
		NumContact[bin_idx] <- sum(interaction.data[si:ei,Contact_Col])
		avg_contact[bin_idx] <- NumContact[bin_idx] / no_distinct_loci[bin_idx]
		prior_contact_prob[bin_idx] <- (avg_contact[bin_idx] / TotContact)
		avg_int_dist[bin_idx] <- mean(gene.dist[si:ei])
	}
} else {
	# equal occupancy (contact count) bin - preferred

	# keeps track of total no of contact counts for a particular bin
	cumContactCount <- 0
	nelem <- 0
	ei <- 0
	bin_idx <- 0

	# keeps track of the no of different genomic distance values encountered within this bin
	# and also the no of elements (locus pairs) having that particular genomic distance
	ValGeneDist <- c()
	nelemSameGeneDist <- c()

	while (ei < numPairs) {
		si <- ei + 1
		curr_gene_dist <- abs(interaction.data[si,2] - interaction.data[si,5])
		idx_list <- which(gene.dist == curr_gene_dist)
		nelem <- nelem + length(idx_list)
		ei <- ei + length(idx_list)
		cumContactCount <- cumContactCount + sum(interaction.data[si:ei, Contact_Col])
		
		ValGeneDist <- c(ValGeneDist, curr_gene_dist)
		nelemSameGeneDist <- c(nelemSameGeneDist, length(idx_list))

		if (cumContactCount >= ncontactbin) {
			# current cumulative contact exceeds the expected average contact per bin
			# and also, all the contacts with this specified genomic distance is covered
			# so, we fix this bin 
			bin_idx <- bin_idx + 1
			no_distinct_loci[bin_idx] <- nelem
			NumContact[bin_idx] <- cumContactCount
			avg_contact[bin_idx] <- cumContactCount / nelem
			prior_contact_prob[bin_idx] <- (avg_contact[bin_idx] / TotContact)

			avg_int_dist[bin_idx] <- 0
			for (j in (1:length(ValGeneDist))) {
				avg_int_dist[bin_idx] <- avg_int_dist[bin_idx] + ((nelemSameGeneDist[j] * 1.0) / nelem) * ValGeneDist[j]
			}

			# reset the couners
			cumContactCount <- 0
		 	ValGeneDist <- c()
		 	nelemSameGeneDist <- c()
		 	nelem <- 0
		}
	}

	# now dump the data in a text file
	OutBinfile <- paste0(OutIntDir, '/', 'Bin_Info.log')
	write.table(cbind(no_distinct_loci, avg_int_dist, NumContact, avg_contact, prior_contact_prob), OutBinfile, row.names = FALSE, col.names = c("NumLoci","AvgDist","NumContact","AvgContact","PriorProb"), sep = "\t", quote=FALSE, append=FALSE) 

}

# now prepare a spline fit where x axis is the 'avg_int_dist'
# and y axis is the 'prior_contact_prob'

# the first spline works according to the specified degree of freedom (16 - for testing)
fit <- smooth.spline(avg_int_dist, prior_contact_prob, df=16)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time for spline (df=16): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# the second spline works according to the cross validation principle
fit2 <- smooth.spline(avg_int_dist, prior_contact_prob, cv=TRUE)
	
if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time for spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# perform anti-tonic regression on the first spline
fit.mr <- monoreg(fit$x, fit$y, type="antitonic")

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time for antitonic regression on spline (df = 16): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# perform anti-tonic regression on the second spline
fit2.mr <- monoreg(fit2$x, fit2$y, type="antitonic")
	
if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time for antitonic regression on spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

# plot the data, if plotting is enabled via the parameter "drawfig"
if (drawfig == 1) {
	pdf(plotfile, width=14, height=10)
	plot(avg_int_dist, prior_contact_prob, cex=0.5, col="darkgrey", xlab="Average interaction distance", ylab="Prior contact probability")
	title("Smooth spline - antitonic regression - pass 1")
	lines(fit.mr, col="red", lwd=2)
	lines(fit2.mr, col="blue", lwd=2)
	legend("topright",legend=c("16 DF", paste(as.integer(fit2$df), "df (cv)")), col=c("red", "blue"), lty=1, lwd=2, cex=0.8)
	dev.off()		
}
	
#--------------------------------------------
# Paper: Ferhat Ay, 2014
# the probability is computed with respect to the spline plot
# which is then substituted to the binomial distribution, for P-value estimation

if (is.unsorted(gene.dist) == TRUE) {

	cat(sprintf("\n **** The interaction file is unsorted --- check ****"))

	# here the genomic distance is unsorted
	# so we employ the predict function and probability distribution computation
	# for the complete array

	# get the spline interpolation prediction from the already modeled smoothing spline object (CV spline)
	# with respect to the input genomic distance vector
	# the returned value is a component of 2 fields: (x,y)
	# where x is the input data and y is the predicted data
	
	# input will be the genomic distance
	# output is the corresponding contact probability 

	# note: the probability was computed with respect to 'numPairs' placed as the denominator
	# accordingly we modify the binomial distribution equation

	pp <- predict(fit2, gene.dist)
	for (k in (1:numPairs)) {
		# probability for this particular interaction is in the y field
		p <- pp$y[k]
		if (BiasCorr == 0) {
			# model binomial distribution with this probability
			Spline_Binom_Prob_CC[k] <- dbinom(interaction.data[k, Contact_Col], size=TotContact, prob=p)
			# compute the p value with respect to the modified distribution
			Spline_Binom_P_Val_CC[k] <- pbinom(interaction.data[k, Contact_Col], size=TotContact, prob=p, lower.tail=FALSE) + Spline_Binom_Prob_CC[k]
			# store the prior probability
			Prob_Val[k] <- p
		} else {
			# model binomial distribution with the product of probability, and the bias values of two segments
			curr_prob <- p * interaction.data[k, (Contact_Col + 1)] * interaction.data[k, (Contact_Col + 3)]
			Spline_Binom_Prob_CC[k] <- dbinom(interaction.data[k, Contact_Col], size=TotContact, prob=curr_prob)
			# compute the p value with respect to the modified distribution
			Spline_Binom_P_Val_CC[k] <- pbinom(interaction.data[k, Contact_Col], size=TotContact, prob=curr_prob, lower.tail=FALSE) + Spline_Binom_Prob_CC[k]
			# store the prior probability
			Prob_Val[k] <- p
		}
	}
} else {

	if (timeprof == 1) {
		starttime <- Sys.time()
	}

	# here the genomic distance is sorted
	# so we selectively apply the predict function to only the distinct elements of the list
	uniq.idx <- order(gene.dist)[!duplicated(gene.dist)]

	if (timeprof == 1) {
		endtime <- Sys.time()
		fp <- file(timefile, open="a")
		outstr <- paste('\n Time to find unique indices in genomic distance: ', (endtime - starttime))
		write(outstr, file=fp, append=T)
		close(fp)
	}

	if (timeprof == 1) {
		starttime <- Sys.time()
	}

	for (k in (1:length(uniq.idx))) {
		# all the locations from si to ei will have the same probability 
		# since the gene distance value is the same
		si <- uniq.idx[k]
		if (k < length(uniq.idx)) {
			ei <- uniq.idx[k+1] - 1
		} else {
			ei <- numPairs	# last read	
		}
		# compute the probability according to the genomic distance (of the start location - single element)
		# this probability value will be used for all the values within the interval (si to ei)
		pp <- predict(fit2, gene.dist[si])
		p <- pp$y

		#==========================
		# compute the probability distribution for the current interval
		#==========================
		if (BiasCorr == 0) {
			curr_dbinom <- dbinom(interaction.data[si:ei, Contact_Col], size=TotContact, prob=p)
			curr_pbinom <- pbinom(interaction.data[si:ei, Contact_Col], size=TotContact, prob=p, lower.tail=FALSE)
			# append the distribution to the final vector
		 	Spline_Binom_Prob_CC <- c(Spline_Binom_Prob_CC, curr_dbinom)
		 	Spline_Binom_P_Val_CC <- c(Spline_Binom_P_Val_CC, (curr_dbinom + curr_pbinom))	
		 	# store the prior probability as well
		 	Prob_Val[si:ei] <- p
		} else {
			# bias based probability correction

			# comment - sourya
			# for (idx in (si:ei)) {
			# 	# model binomial distribution with the product of probability, and the bias values of two segments
			# 	curr_prob <- p * interaction.data[idx, (Contact_Col + 1)] * interaction.data[idx, (Contact_Col + 3)]
			# 	curr_dbinom <- dbinom(interaction.data[idx, Contact_Col], size=TotContact, prob=curr_prob)
			# 	curr_pbinom <- pbinom(interaction.data[idx, Contact_Col], size=TotContact, prob=curr_prob, lower.tail=FALSE)
			# 	# append the distribution to the final vector
			#  	Spline_Binom_Prob_CC[idx] <- curr_dbinom
			#  	Spline_Binom_P_Val_CC[idx] <- (curr_dbinom + curr_pbinom)	
			#  	# store the prior probability as well
		 		# Prob_Val[idx] <- curr_prob
			# }

			# add - sourya

			# cat(sprintf("\n *** Before executing Parallel from si: %s to ei: %s ", si, ei))

			# now apply this matrix on the parallel version of binomial distribution computation
			# the second argument 1 means that individual rows of matrix is processed 
			# in the parallel computation
			result_binomdistr <- as.data.frame(parallel:::mclapply( si:ei , mc.cores = ncore , function(idx){
				cc <- interaction.data[idx, Contact_Col]
				pr <- p * interaction.data[idx, (Contact_Col + 1)] * interaction.data[idx, (Contact_Col + 3)]
				db <- dbinom(cc, size=TotContact, prob=pr)
				pb <- pbinom(cc, size=TotContact, prob=pr, lower.tail=FALSE)
				return(c(idx,pr,db,(db+pb)))
				} ))

			# cat(sprintf("\n *** rows: %s   columns: %s ", nrow(result_binomdistr), ncol(result_binomdistr)))

			# copy the results
			# columns : si to ei
			# rows: 4 
			# for (i in (1:ncol(result_binomdistr))) {
			# 	idx <- result_binomdistr[1,i]
			# 	Prob_Val[idx] <- result_binomdistr[2,i]
			# 	Spline_Binom_Prob_CC[idx] <- result_binomdistr[3,i]
			# 	Spline_Binom_P_Val_CC[idx] <- result_binomdistr[4,i]
			# }
			Prob_Val[si:ei] <- as.double(result_binomdistr[2,1:ncol(result_binomdistr)])
			Spline_Binom_Prob_CC[si:ei] <- as.double(result_binomdistr[3,1:ncol(result_binomdistr)])
			Spline_Binom_P_Val_CC[si:ei] <- as.double(result_binomdistr[4,1:ncol(result_binomdistr)])

			# end add - Sourya
		}
	}

	if (timeprof == 1) {
		endtime <- Sys.time()
		fp <- file(timefile, open="a")
		outstr <- paste('\n Time for Spline based distribution compute (p value) (sorted interaction file): ', (endtime - starttime))
		write(outstr, file=fp, append=T)
		close(fp)
	}
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# from the generated P values, obtain the Q value using BH correction
Spline_Binom_QVal <- p.adjust(Spline_Binom_P_Val_CC, method = "BH")	

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time to compute Q value (from P value) in spline: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# accumulate all results - also add header information
FinalData <- cbind(interaction.data, Prob_Val, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)
colnames(FinalData) <- c(colnames(interaction.data), "p", "dbinom", "P-Value", "Q-Value")

# append the Spline distribution probability and corresponding P value as separate columns
# and write in a separate text file
temp_outfile <- paste0(inpdir, '/', 'temp_out.bed')
write.table(FinalData, temp_outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 

# now sort the file contents and write that in the final specified output file
system(paste('sort -k1,1 -k2,2n -k5,5n', paste0('-k',Contact_Col,',',Contact_Col,'nr'), temp_outfile, '>', outfile))
# delete the temporary output file
system(paste('rm', temp_outfile))

# # sort the data according to the chromosome name, start positions of both the interacting segments
# # and finally with the decreasing contact count
# FinalData <- FinalData[order( FinalData[,1], FinalData[,2], FinalData[,5], -FinalData[,Contact_Col] ), ]

# # write the data with the header information
# write.table(FinalData, outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(timefile, open="a")
	outstr <- paste('\n Time to write the interaction data (with Q value) for spline based distribution: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}
