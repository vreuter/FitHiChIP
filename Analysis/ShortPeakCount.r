#!/usr/bin/env Rscript

#===========================================================
# R script for checking the no of contacts (both filtered and unfiltered) for each peaks 
# also chromosome specific contact counts are listed

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript ShortPeakCount.r $InpSplineBeforeFilterFile $Q_Thr
#===========================================================

library(ggplot2)

args <- commandArgs(TRUE)
InteractionSplineFile <- args[1]

# column no having the absolute contact count
ContactCol <- as.integer(args[2])

# global FDR threshold employed
FDR.thres <- as.double(args[3])

InpDir <- dirname(InteractionSplineFile)
OutDir <- paste0(InpDir, '/Results')
system(paste('mkdir -p', OutDir))

# List of chromosomes
ChrList <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")

# The unfiltered interaction file contains header information
IntSplineData <- read.table(InteractionSplineFile, header=T)

# for individual chromosomes, check the corresponding peak locations
for (chr_idx in (1:length(ChrList))) {

	curr.chr <- ChrList[chr_idx]

	PlotFile1 <- paste0(OutDir, '/', curr.chr, '_Contact_Count_Distr.png')
	PlotFile2 <- paste0(OutDir, '/', curr.chr, '_Distinct_Interaction_Distr.png')

	if ((file.exists(PlotFile1) == FALSE) || (file.exists(PlotFile2) == FALSE)) {

		# filter and get the cis interaction data only for the current chromosomes
		IntData.CurrChr <- IntSplineData[which(IntSplineData[,1] == curr.chr),]

		if (nrow(IntData.CurrChr) > 0) {
			#cat(sprintf("\n Scanning the information of %s \n", curr.chr))
			
			# this vector stores the midpoints of the peaks
			Traversed_Peaks <- c()
			# this vector stores the contact count considering all interactions (filtered or unfiltered) of the current peak
			CC_Unfiltered <- c()
			# this vector stores the contact count of filtered (FDR threshold) contacts of the current peak
			CC_filtered <- c()
			# this vector stores only the no of distinct interactions for the current peak
			NumInt_Unfiltered <- c()
			# this vector stores only the no of distinct interactions (FDR threshold) for the current peak
			NumInt_filtered <- c()

			# index of points
			j <- 0	

			# now traverse individual peaks of the current chromosome
			for (i in (1:nrow(IntData.CurrChr))) {
				curr_peak_midpoint <- as.integer((IntData.CurrChr[i,2] + IntData.CurrChr[i,3]) / 2)
				if ((curr_peak_midpoint %in% Traversed_Peaks) == FALSE) {
					j <- j + 1
					#cat(sprintf("\n new point: %d ", curr_peak_midpoint))
					
					# append the midpoint to the list of peaks
					Traversed_Peaks[j] <- curr_peak_midpoint
					
					# store the no of contacts (interactions) without FDR based filtering
					x <- which((IntData.CurrChr[,2] == IntData.CurrChr[i,2]) & (IntData.CurrChr[,3] == IntData.CurrChr[i,3]))
					if (length(x) > 0) {
						ncontact_unfilter <- sum(IntData.CurrChr[x,ContactCol])
					} else {
						ncontact_unfilter <- 0
					}
					CC_Unfiltered[j] <- ncontact_unfilter
					NumInt_Unfiltered[j] <- length(x)

					# store the no of contacts (interactions) with FDR based filtering
					# the last column stores the FDR value
					x <- which((IntData.CurrChr[,2] == IntData.CurrChr[i,2]) & (IntData.CurrChr[,3] == IntData.CurrChr[i,3]) & (IntData.CurrChr[,ncol(IntData.CurrChr)] < FDR.thres))
					if (length(x) > 0) {
						ncontact_filter <- sum(IntData.CurrChr[x,ContactCol])
					} else {
						ncontact_filter <- 0
					}
					CC_filtered[j] <- ncontact_filter
					NumInt_filtered[j] <- length(x)
				}
			}

			#cat(sprintf("\n No of valid points: %d ", j))

			# now plot the midpoints along with the contact counts (both filtered and unfiltered)
			# in a graph
			a1 <- data.frame(group = paste0("Contact count: Qval<", FDR.thres), x = Traversed_Peaks, y = CC_filtered)
			b1 <- data.frame(group = "Contact count: ALL", x = Traversed_Peaks, y = CC_Unfiltered)
			curr_plot1 <- ggplot(rbind(a1,b1), aes(x=x, y=y, fill=group, colour=group)) + geom_point()
			curr_plot1 + ggtitle("Individual peaks - Contact count - filtered (Qvalue) vs Unfiltered")
			ggsave(PlotFile1, width=20, height=8)

			# now plot the midpoints along with the contact counts (both filtered and unfiltered)
			# in a graph
			a2 <- data.frame(group = paste0("Distinct interactions: Qval<", FDR.thres), x = Traversed_Peaks, y = NumInt_filtered)
			b2 <- data.frame(group = "Distinct interactions: ALL", x = Traversed_Peaks, y = NumInt_Unfiltered)
			curr_plot2 <- ggplot(rbind(a2,b2), aes(x=x, y=y, fill=group, colour=group)) + geom_point()
			curr_plot2 + ggtitle("Individual peaks - No of distinct interactions - filtered (Qvalue) vs Unfiltered")
			ggsave(PlotFile2, width=20, height=8)
		}
	}
}