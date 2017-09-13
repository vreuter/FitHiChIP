#!/usr/bin/env Rscript

#===========================================================
# R script for analyzing the distribution of contacts (according to sorted genomic distance) 
# Divides the contacts into equal length (with respect to genomic distance) intervals, and then 
# counts the number of contacts in these intervals

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript binquantile.r $inpfile
#===========================================================

# process input arguments
args <- commandArgs(TRUE)

inpfile <- args[1]
inpdir <- dirname(inpfile)

# no of bins (intervals)
nbins <- as.integer(args[2])

# column of the absolute contact count
contactcol <- as.integer(args[3])

# output file for plotting
plotfile <- args[4]

outdir <- dirname(plotfile)
system(paste('mkdir -p', outdir))

# load the interaction data matrix
# Note: the interaction data will full features has header information
intdata <- read.table(inpfile, header=T)

# absolute genomic distance for an interaction instance
c1 <- abs(intdata[,2] - intdata[,5])

# range of genomic distance
rd <- max(c1) - min(c1)
# interval points for each bin
bin.interval.prcnt <- seq(0,1,by=1/nbins)*rd+min(c1)

# check the mean contact count within each intervals
# Note: some bin may have zero contact count as well
# and some bin may have a huge number of entries
cc.intrv <- c()
for (i in (1:length(bin.interval.prcnt))) {
	if (i==1) {
		ContactInfo <- intdata[,contactcol]
		vec <- ContactInfo[which(c1 <= bin.interval.prcnt[i])]
		#length(vec)
		if (length(vec) > 0) {
			cc.intrv[i] <- mean(vec)
		} else {
			cc.intrv[i] <- 0
		}	
	} else {
		ContactInfo <- intdata[,contactcol]
		vec <- ContactInfo[which((c1 > bin.interval.prcnt[i-1]) & (c1 <= bin.interval.prcnt[i]))]
		#length(vec)
		if (length(vec) > 0) {
            cc.intrv[i] <- mean(vec)
        } else {
            cc.intrv[i] <- 0
        }
	}
}	
# output pdf file which will store the plot
pdf(plotfile, width=14, height=10)
plot(cc.intrv, main="Binned contact count (mean) - equal length interval", xlab="bin no", ylab="mean contact count", col="red")
dev.off()
