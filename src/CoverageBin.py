#!/usr/bin/env python

"""
Created: September 21,2017

This file takes input a HiC-pro pipeline output of validpairs file, 
Chromosome length file, and a peak detection file
it generates binned intervals from the input chromosomes
and creates coverage vector for individual chromosomes
Dependin on the input peak file, it also marks individual 
binned segments as either peak or non-peak

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import sys
import os
from optparse import OptionParser

# used to open the gzip compressed text file
import gzip 

# these two modules are used to process the pairix and bam input files
import pypairix
import pysam

#==========================================================
""" 
this class defines a particular chromosome interval
indexed by (chromosome, start, end) - dictionary key 
"""
class ChrmSegment(object):  
    def __init__(self):
        self.is_peak = False
        self.read_count = 0

    def _CheckIfPeak(self):
        return self.is_peak

    def _GetReadCount(self):
        return self.read_count

    def _IncrementRead(self, count=1):
        self.read_count = self.read_count + count

    def _SetPeak(self):
        self.is_peak = True

#===============================================
def main():
    parser = OptionParser() #(usage=usage)
    parser.add_option("-p", "--peakfile", dest="peakfile", help="Peak detection file")
    parser.add_option("-i", "--inpfile", dest="inpfile", help="Input valid pairs file (from HiC-pro pipeline)")
    parser.add_option("-b", "--binsize", dest="binsize", type="int", help="Size of bins employed. DEFAULT 5000 (indicating 5 Kb).")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file storing the coverage of individual genomic bins")
    parser.add_option("-c", "--chrsizefile", dest="chrsizefile", help="File containing chromosome size information")
    parser.add_option("-f", "--filefmt", dest="filefmt", type="int", help="Input file format (1 = BAM file, 2 = pairix file, 3 = HiC-Pro pipeline output. Default 3.")    
    parser.set_defaults(peakfile=None, inpfile=None, binsize=5000, filefmt=3, outfile=None, chrsizefile=None)
    (options, args) = parser.parse_args()

    # global variables from command line options
    global inpfile
    global peakfile
    global outfile
    global bin_size
    global chrsizefile
    global input_file_fmt

    # dictionary structure which stores the instances of the class ChrmSegment
    global SegmentDict
    SegmentDict = dict()

    # two lists which store the chromosome name and their length information
    # with respect to the input chromosome length file
    ChrNameList_Ref = []
    ChrLenList_Ref = []

    # list storing the chromosome names provided in the input validpairs file
    # ChrNameList_Current = [] 

    # input file format
    input_file_fmt = int(options.filefmt)

    if options.inpfile is not None:
        inpfile = options.inpfile
    else:
        sys.exit("Input validpairs file (output of HiC-pro pipeline) is not provided - quit !!")

    if options.peakfile is not None:
        peakfile = options.peakfile
    else:
        sys.exit("Peak detection file is not provided - quit !!")

    if options.outfile is not None:
        outfile = options.outfile
    else:
        sys.exit("Output file for storing the coverage of segments is not provided - quit !!")

    if options.chrsizefile is not None:
        chrsizefile = options.chrsizefile
    else:
        if (input_file_fmt == 1) or (input_file_fmt == 2):
            chrsizefile = None
        else:
            sys.exit("Chromosome size file is not provided even if HiC-pro pipeline based valid pairs file is the input - quit !!")

    # input bin size parameter  
    bin_size = int(options.binsize)

    #=====================================================
    # first read the reference chromosome length file
    # and assign the name and lengths of individual chromosomes
    #=====================================================
    if (input_file_fmt == 1):
        # BAM alignment file
        bam_input_file = pysam.AlignmentFile(inpfile, "rb")
        bam_input_refchr_names = list(bam_input_file.references)
        bam_input_refchr_lengths = list(bam_input_file.lengths)
        for i in range(len(bam_input_refchr_names)):
            ChrNameList_Ref.append(bam_input_refchr_names[i])
        ChrNameList_Ref.sort()
        for i in range(len(ChrNameList_Ref)):
            ChrLenList_Ref.append(bam_input_refchr_lengths[bam_input_refchr_names.index(ChrNameList_Ref[i])])
        bam_input_file.close()
    
    elif (input_file_fmt == 2):
        # PAIRIX alignment file
        align_tb = pypairix.open(inpfile)
        x = align_tb.get_chromsize()
        for i in range(len(x)):
            ChrNameList_Ref.append(x[i][0])
            ChrLenList_Ref.append(int(x[i][1]))
    
    else:
        # HiC pro pipeline output valid pairs file
        with open(chrsizefile, 'r') as fp:
            for line in fp:
                linecontents = (line.rstrip()).split()
                curr_chr = linecontents[0]
                curr_chr_size = int(linecontents[1])
                ChrNameList_Ref.append(curr_chr)
                ChrLenList_Ref.append(curr_chr_size)

    print '\n Created ChrNameList: ', ChrNameList_Ref
    print '\n Created ChrLenList: ', ChrLenList_Ref

    #=====================================================
    # now create individual genomic segments (dictionary and class instance)
    # according to the maximum chromosome length specified for individual chromosomes
    #=====================================================
    for i in range(len(ChrNameList_Ref)):
        curr_chr = ChrNameList_Ref[i]
        curr_chr_size = ChrLenList_Ref[i]
        if ((curr_chr_size % bin_size) == 0):
            interval_end = curr_chr_size
        else:
            interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
        for val in range(0, interval_end, bin_size):
            curr_key = (curr_chr, val, (val + bin_size))
            # python class instance for this key
            SegmentDict.setdefault(curr_key, ChrmSegment())

    #=====================================================
    # now read the input alignment (in either BAM or PAIRIX format) / valid pairs file
    # and compute the coverage information
    # with respect to 'cis' reads
    #=====================================================

    # scan through all the reads and then update the read count
    
    if (input_file_fmt == 1):
        # BAM file alignment
        bam_input_file = pysam.AlignmentFile(inpfile, "rb")
        for read in bam_input_file:
            # check if it is paired end read
            if (read.is_paired == True):
                # string representation of individual alignment
                read_str = (read.tostring(bam_input_file)).split()
                if (str(read_str[2]) == str(read_str[6])) or (str(read_str[6]) == "="):
                    # cis interaction
                    chrname = str(read_str[2])
                    # get the 1st read
                    # increment the read count in the corresponding dictionary entry of this bin
                    read1_pos = int(read_str[3])
                    bin_start1 = (int(read1_pos / bin_size)) * bin_size
                    curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
                    if curr_key1 in SegmentDict:
                        SegmentDict[curr_key1]._IncrementRead()
                    # get the 2nd read
                    # increment the read count in the corresponding dictionary entry of this bin
                    read2_pos = int(read_str[7])
                    bin_start2 = (int(read2_pos / bin_size)) * bin_size
                    curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
                    if curr_key2 in SegmentDict:
                        SegmentDict[curr_key2]._IncrementRead()
        # close the bam alignment file
        bam_input_file.close()
    
    elif (input_file_fmt == 2):
        # PAIRIX file alignment
        with gzip.open(inpfile,'r') as fin:    
            for line in fin:
                # ignore the comment lines
                if (line[0] != '#'):
                    contents = line.split()
                    if (str(contents[1]) == str(contents[3])):
                        # cis interaction
                        chrname = str(contents[1])
                        # get the 1st read
                        # increment the read count in the corresponding dictionary entry of this bin
                        read1_pos = int(contents[2])
                        bin_start1 = (int(read1_pos / bin_size)) * bin_size
                        curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
                        if curr_key1 in SegmentDict:
                            SegmentDict[curr_key1]._IncrementRead()
                        # get the 2nd read
                        # increment the read count in the corresponding dictionary entry of this bin
                        read2_pos = int(contents[4])
                        bin_start2 = (int(read2_pos / bin_size)) * bin_size
                        curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
                        if curr_key2 in SegmentDict:
                            SegmentDict[curr_key2]._IncrementRead()

    else:
        # HiC pro pipeline output valid pairs file
        with gzip.open(inpfile,'r') as fin:    
            for line in fin:
                contents = line.split()
                # check cis interaction
                if (str(contents[1]) == str(contents[4])):
                    chrname = str(contents[1])
                    read1_pos = int(contents[2])
                    read2_pos = int(contents[5])
                                        
                    # now process the current paired end read
                    # read position 1
                    bin_start1 = (int(read1_pos / bin_size)) * bin_size
                    # increment the read count in the corresponding dictionary entry of this bin
                    curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
                    if curr_key1 in SegmentDict:
                        SegmentDict[curr_key1]._IncrementRead()
                    # read position 2
                    bin_start2 = (int(read2_pos / bin_size)) * bin_size
                    # increment the read count in the corresponding dictionary entry of this bin
                    curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
                    if curr_key2 in SegmentDict:
                        SegmentDict[curr_key2]._IncrementRead()
    
    #=====================================================
    # scan the peak input file (if provided) and mark the corresponding 
    # chromosome intervals as 1 (having peak)
    #=====================================================
    with open(peakfile, 'r') as fp:
        for line in fp:
            linecontents = (line.rstrip()).split()
            curr_chr = linecontents[0]
            peak_start = int(linecontents[1])
            peak_end = int(linecontents[2])
            # only process those peaks whose chromosome 
            # is present in the current valid pairs file
            # this is possible when the peak information is downloaded from a reference
            if curr_chr in ChrNameList_Ref:
                interval_start = (int(peak_start / bin_size)) * bin_size
                if ((peak_end % bin_size) == 0):
                    interval_end = peak_end
                else:
                    interval_end = (int(peak_end / bin_size) + 1) * bin_size
                for val in range(interval_start, interval_end, bin_size):
                    curr_key = (curr_chr, val, (val + bin_size))
                    # mark the corresponding dictionary entry as a peak segment
                    # since it has an overlap with the MACS2 derived peak intervals
                    SegmentDict[curr_key]._SetPeak()
    
    #=====================================================
    # now write the genomic intervals and corresponding coverage + isPeak information
    #=====================================================
    fp_cov = open(outfile, 'w')
    fp_cov.write('Chr' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Coverage' + '\t' + 'IsPeak')
    for i in range(len(ChrNameList_Ref)):
        curr_chr = ChrNameList_Ref[i]
        curr_chr_size = ChrLenList_Ref[i]
        if ((curr_chr_size % bin_size) == 0):
            interval_end = curr_chr_size
        else:
            interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
        for val in range(0, interval_end, bin_size):
            curr_key = (curr_chr, val, (val + bin_size))
            fp_cov.write('\n' + str(curr_chr) + '\t' + str(val) + '\t' + str(val + bin_size) + '\t' + str(SegmentDict[curr_key]._GetReadCount()) + '\t' + str(int(SegmentDict[curr_key]._CheckIfPeak())))
    fp_cov.close()

    return

#===============================================
if __name__ == "__main__":
    main()
