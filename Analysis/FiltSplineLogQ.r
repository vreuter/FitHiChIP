#!/usr/bin/env Rscript

#===========================================================
# R script to create a file containing interaction and the log10 Q value
# corresponding to the statistical significant interactions

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================
args <- commandArgs(TRUE)
inpfile <- args[1]
outfile <- args[2]

x=read.table(inpfile, header=FALSE)
# the Q value is placed at the 13th field of the input (spline fitted filtered) interaction file
# which also includes the normalized contact count
Q=-log10(x[,ncol(x)])

# write the data in a table
# Note: the interaction file has the following specified format:
# ch1,start1,end1	chr2,start2,end2	contactcount
write.table(cbind(paste(as.character(x[,1]),x[,2],x[,3],sep=','),paste(as.character(x[,4]),x[,5],x[,6],sep=','),Q), outfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
