#!/usr/bin/env python

"""
Created: August 21,2017

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import sys
import os

# these two modules are used to process the pairix and bam input files
import pypairix
import pysam

# for multi thread environment
import multiprocessing
import threading

# import re
# import time
# import math

from optparse import OptionParser

# used to open the gzip compressed text file
import gzip 

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
    parser.add_option("-p", "--peakfile", dest="peakfile", help="Peak detection output file")
    parser.add_option("-M", "--Method", dest="Method", type="int", help="Type of interaction to be derived (1 - peak to peak, 2 - peak to non peak, 3 - all to all). Default 3.")
    parser.add_option("-l", "--alignfile", dest="alignfile", help="Input alignment file (in BAM or PAIRIX format)")
    parser.add_option("-o", "--outfile", dest="outfile", help="output bed filename storing the interactions")
    # parser.add_option("-U", "--upperbound", dest="distUpThres", type="int", help="OPTIONAL: upper bound on the intra-chromosomal distance range. DEFAULT 200000.")
    # parser.add_option("-L", "--lowerbound", dest="distLowThres", type="int", help="OPTIONAL: lower bound on the intra-chromosomal distance range. DEFAULT 20000.")
    parser.add_option("-t", "--threads", dest="threads", help="Number of threads employed for computation. DEFAULT 1.")
    parser.add_option("-b", "--binsize", dest="binsize", type="int", help="Size of bins employed. DEFAULT 5 Kb.")
    parser.add_option("-f", "--filefmt", dest="filefmt", type="int", help="Input file format (1 = BAM file, 2 = pairix file. Default 1.")  
    # parser.add_option("-c", "--coveragefile", dest="coveragefile", help="output bed filename storing the coverage of individual genomic segments")
    # parser.set_defaults(peakfile=None, alignfile=None, outfile=None, distLowThres=20000, distUpThres=200000, threads=1, binsize=5000, Method=3, filefmt=1, coveragefile=None)
    # parser.set_defaults(peakfile=None, alignfile=None, outfile=None, threads=1, binsize=5000, Method=3, filefmt=1, coveragefile=None)
    parser.set_defaults(peakfile=None, alignfile=None, outfile=None, threads=1, binsize=5000, Method=3, filefmt=1)
    (options, args) = parser.parse_args()

    # global variables from command line options
    # global distUpThres
    # global distLowThres
    global alignfile
    global peakfile
    # global coveragefile

    global threads
    global bin_size
    global Method
    global input_file_fmt

    # dictionary structure which stores the instances of the class ChrmSegment
    global SegmentDict
    SegmentDict = dict()

    # two lists which store the chromosome name and their length information
    ChrNameList = []
    ChrLenList = []

    if options.alignfile is not None:
        alignfile = options.alignfile
    else:
        sys.exit("Alignment file is not provided - quit !!")

    # output file storing the bed formatted interactions
    if options.outfile is not None:
        outfile = options.outfile
    else:
        sys.exit("Output file (for storing the interactions) is not specified - quit !!")

    # sourya
    # if options.peakfile is not None:
    peakfile = options.peakfile

    # if options.coveragefile is not None:
    #     coveragefile = options.coveragefile

    # # lower and upper distance threshold for the interactions
    # distUpThres = int(options.distUpThres)
    # distLowThres = int(options.distLowThres)

    # other parameters' assignment
    threads = int(options.threads)
    bin_size = int(options.binsize)
    Method = int(options.Method)
    input_file_fmt = int(options.filefmt)

    if (Method != 3) and (options.peakfile is None):
        sys.exit("Peak based interaction is sought but no peak file is provided - exit!!")

    #=====================================================
    # create python dictionary containing the chromosome and bin intervals
    #=====================================================

    # first read the chromosome name and its lengths from the input alignment file
    if (input_file_fmt == 1):
        # BAM alignment file
        bam_input_file = pysam.AlignmentFile(alignfile, "rb")
        bam_input_refchr_names = list(bam_input_file.references)
        bam_input_refchr_lengths = list(bam_input_file.lengths)
        for i in range(len(bam_input_refchr_names)):
            ChrNameList.append(bam_input_refchr_names[i])
        ChrNameList.sort()
        for i in range(len(ChrNameList)):
            ChrLenList.append(bam_input_refchr_lengths[bam_input_refchr_names.index(ChrNameList[i])])
        bam_input_file.close()
    else:
        # PAIRIX alignment file
        align_tb = pypairix.open(alignfile)
        x = align_tb.get_chromsize()
        for i in range(len(x)):
            ChrNameList.append(x[i][0])
            ChrLenList.append(int(x[i][1]))

    print '\n Created ChrNameList: ', ChrNameList
    print '\n Created ChrLenList: ', ChrLenList

    # now create individual genomic segments (dictionary and class instance)
    # according to the maximum chromosome length specified for individual chromosomes
    for i in range(len(ChrNameList)):
        curr_chr = ChrNameList[i]
        curr_chr_size = ChrLenList[i]
        if ((curr_chr_size % bin_size) == 0):
            interval_end = curr_chr_size
        else:
            interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
        for val in range(0, interval_end, bin_size):
            curr_key = (curr_chr, val, (val + bin_size))
            # python class instance for this key
            SegmentDict.setdefault(curr_key, ChrmSegment())

    #=====================================================
    # scan the peak input file (if provided) and mark the corresponding 
    # chromosome intervals as 1 (having peak)
    #=====================================================
    # sourya - condition is modified
    if (peakfile is not None):   #(Method != 3):
        with open(peakfile, 'r') as fp:
            for line in fp:
                linecontents = (line.rstrip()).split()
                curr_chr = linecontents[0]
                peak_start = int(linecontents[1])
                peak_end = int(linecontents[2])
                # constraint the peak end with respect to the maximum 
                # length of this chromosome
                # provided the chromosome information exists within the current chromosome
                if curr_chr in ChrNameList:
                    peak_end = min(peak_end, ChrLenList[ChrNameList.index(curr_chr)])
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
    # Now create a list of chromosome intervals for interaction calling
    # for method = 1 or 2, peak segments (bins) will be used
    # for method = 3, all segments (bins) will be used
    # Note: currently we take input the chromosome size file (two column)
    # chromosome name and the length information
    #=====================================================
    List_Chr_Regions = []
    for i in range(len(ChrNameList)):
        curr_chr = ChrNameList[i]
        curr_chr_size = ChrLenList[i]
        if ((curr_chr_size % bin_size) == 0):
            interval_end = curr_chr_size
        else:
            interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
        for val in range(0, interval_end, bin_size):
            curr_key = (curr_chr, val, (val + bin_size))
            if (Method == 3) or (SegmentDict[curr_key]._CheckIfPeak() == True):
                subl = [curr_key[0], curr_key[1], curr_key[2]]
                List_Chr_Regions.append(subl)

    # apply a multi thread environment to process these anchor regions of short range segments
    # to find the interactions with the long range segments
    p = multiprocessing.Pool(threads)
    try:
        mapped_long_reads = p.map_async(MapShort2LongSegment, List_Chr_Regions).get(9999999)
    except KeyboardInterrupt:
        p.terminate()
        p.join()

    # open the output interaction bed file
    fout = open(outfile, 'w')

    # process the mapped long segment information (contacts)
    for mappedread_anchor in mapped_long_reads:
        # mappedread_anchor is a list of contacts for a particular anchor region
        for elem in mappedread_anchor:
            if (len(elem) > 0):
                # comment - sourya
                # #==============================================
                # # for peak to peak or for ALL to ALL interactions, coordinate based ordering is used
                # if ((Method == 1) or (Method == 3)):
                #     if (elem[1] < elem[4]):
                #         fout.write(str(elem[0]) + '\t' + str(elem[1]) + '\t' + str(elem[2]) + '\t' + str(elem[3]) + '\t' + str(elem[4]) + '\t' + str(elem[5]) + '\t' + str(elem[6]))
                #         fout.write('\n')
                #     else:
                #         fout.write(str(elem[3]) + '\t' + str(elem[4]) + '\t' + str(elem[5]) + '\t' + str(elem[0]) + '\t' + str(elem[1]) + '\t' + str(elem[2]) + '\t' + str(elem[6]))
                #         fout.write('\n')
                # else:
                #     # for peak to non peak, or peak to ALL interactions
                #     # left side contains the peak information
                #     fout.write(str(elem[0]) + '\t' + str(elem[1]) + '\t' + str(elem[2]) + '\t' + str(elem[3]) + '\t' + str(elem[4]) + '\t' + str(elem[5]) + '\t' + str(elem[6]))
                #     fout.write('\n')
                #==============================================
                # add - sourya
                # we employ coordinate based ordering
                # if the original interaction is also ordered, we insert 0 at the end
                # otherwise, we insert 1 at the end, to indicate the change in ordering
                if (elem[1] < elem[4]):
                    fout.write(str(elem[0]) + '\t' + str(elem[1]) + '\t' + str(elem[2]) + '\t' + str(elem[3]) + '\t' + str(elem[4]) + '\t' + str(elem[5]) + '\t' + str(elem[6]) + '\t' + "0")
                    fout.write('\n')
                else:
                    fout.write(str(elem[3]) + '\t' + str(elem[4]) + '\t' + str(elem[5]) + '\t' + str(elem[0]) + '\t' + str(elem[1]) + '\t' + str(elem[2]) + '\t' + str(elem[6]) + '\t' + "1")
                    fout.write('\n')


    # close the output text file
    fout.close()

    # #==================================================================
    # # if the coverage file is specified
    # # compute the coverage for individual genomic segments
    # # given the alignment file in either BAM or PAIRIX format
    # if options.coveragefile is not None:
    #     # scan through all the reads and then update the read count
    #     if (input_file_fmt == 1):
    #         # BAM file alignment
    #         bam_input_file = pysam.AlignmentFile(alignfile, "rb")
    #         for read in bam_input_file:
    #             # check if it is paired end read
    #             if (read.is_paired == True):
    #                 # string representation of individual alignment
    #                 read_str = (read.tostring(bam_input_file)).split()
    #                 if (str(read_str[2]) == str(read_str[6])) or (str(read_str[6]) == "="):
    #                     # cis interaction
    #                     chrname = str(read_str[2])
    #                     # get the 1st read
    #                     # increment the read count in the corresponding dictionary entry of this bin
    #                     read1_pos = int(read_str[3])
    #                     bin_start1 = (int(read1_pos / bin_size)) * bin_size
    #                     curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
    #                     if curr_key1 in SegmentDict:
    #                         SegmentDict[curr_key1]._IncrementRead()
    #                     # get the 2nd read
    #                     # increment the read count in the corresponding dictionary entry of this bin
    #                     read2_pos = int(read_str[7])
    #                     bin_start2 = (int(read2_pos / bin_size)) * bin_size
    #                     curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
    #                     if curr_key2 in SegmentDict:
    #                         SegmentDict[curr_key2]._IncrementRead()
    #         # close the bam alignment file
    #         bam_input_file.close()
    #     else:
    #         # PAIRIX file alignment
    #         with gzip.open(alignfile,'r') as fin:    
    #             for line in fin:
    #                 # ignore the comment lines
    #                 if (line[0] != '#'):
    #                     contents = line.split()
    #                     if (str(contents[1]) == str(contents[3])):
    #                         # cis interaction
    #                         chrname = str(contents[1])
    #                         # get the 1st read
    #                         # increment the read count in the corresponding dictionary entry of this bin
    #                         read1_pos = int(contents[2])
    #                         bin_start1 = (int(read1_pos / bin_size)) * bin_size
    #                         curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
    #                         if curr_key1 in SegmentDict:
    #                             SegmentDict[curr_key1]._IncrementRead()
    #                         # get the 2nd read
    #                         # increment the read count in the corresponding dictionary entry of this bin
    #                         read2_pos = int(contents[4])
    #                         bin_start2 = (int(read2_pos / bin_size)) * bin_size
    #                         curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
    #                         if curr_key2 in SegmentDict:
    #                             SegmentDict[curr_key2]._IncrementRead()

    #     # now write the genomic intervals and corresponding read count + peak information
    #     fp_cov = open(coveragefile, 'w')
    #     fp_cov.write('Chr' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Coverage' + '\t' + 'IsPeak')
    #     for i in range(len(ChrNameList)):
    #         curr_chr = ChrNameList[i]
    #         curr_chr_size = ChrLenList[i]
    #         if ((curr_chr_size % bin_size) == 0):
    #             interval_end = curr_chr_size
    #         else:
    #             interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
    #         for val in range(0, interval_end, bin_size):
    #             curr_key = (curr_chr, val, (val + bin_size))
    #             if curr_key in SegmentDict:
    #                 fp_cov.write('\n' + str(curr_chr) + '\t' + str(val) + '\t' + str(val + bin_size) + '\t' + str(SegmentDict[curr_key]._GetReadCount()) + '\t' + str(int(SegmentDict[curr_key]._CheckIfPeak())))
    #     fp_cov.close()


#=================================================
def MapShort2LongSegment((anchor_chr, anchor_start, anchor_end)):
    # global distLowThres
    # global distUpThres
    global alignfile
    global bin_size
    global Method
    global SegmentDict
    global input_file_fmt

    if anchor_start > anchor_end: 
        exit('Error: invalid anchor region: ' + '\t'.join([anchor_chr, anchor_start, anchor_end]))

    # dictionary storing the binned sum
    bins = dict()

    # total no of mapped reads
    totread = 0

    # start of the anchor mapped with respect to the fixed size bin interval
    anchor_start_bin = (int((anchor_start * 1.0) / bin_size)) * bin_size

    # dictionary entry corresponding to the anchor (binned) segment
    anchor_key = (anchor_chr, anchor_start_bin, (anchor_start_bin + bin_size))

    if (input_file_fmt == 1):

        # BAM alignment file (input)
        bam_input_file = pysam.AlignmentFile(alignfile, "rb")
        
        # fetch the reads with respect to the specified chromosome interval
        for _read in bam_input_file.fetch(anchor_chr, anchor_start, anchor_end):
            
            # start point of the mapped read
            cur_pos = _read.pos + _read.isize;

            # mapped read start point with respect to fixed size binning
            curr_pos_bin_start = int((cur_pos * 1.0) / bin_size) * bin_size
            
            # key corresponding to the chromosome binned interval (python dictionary)
            # of the mapped segment
            curr_key = (anchor_chr, curr_pos_bin_start, (curr_pos_bin_start + bin_size))

            # check if the mapped region falls within the chromosome region
            if curr_key in SegmentDict:
                #  Note: This condition is commented  -sourya
                # The objective is to get a generic interaction file at first
                # after that, we will filter the interactions according to the distance threshold requirement
                # check the distance threshold
                if (1): #(abs(curr_pos_bin_start - anchor_start_bin) >= distLowThres) and (abs(curr_pos_bin_start - anchor_start_bin) <= distUpThres):
                    # the mapped read is within the distance threshold
                    # methods 3 or 0 consider the other end of the read as either peak or non peak
                    # method 1 considers the other end of the read as peak
                    # method 2 considers the other end of the read as non peak
                    if (Method == 3) or (Method == 0) or ((Method == 1) and (SegmentDict[curr_key]._CheckIfPeak() == True)) or ((Method == 2) and (SegmentDict[curr_key]._CheckIfPeak() == False)):            
                        bin_id = (curr_pos_bin_start - anchor_start_bin) / bin_size 
                        if bin_id not in bins:
                            bins.setdefault(bin_id, 0)
                        bins[bin_id] += 1
                        totread += 1

    else:

        # PAIRIX alignment file (input)
        align_tb = pypairix.open(alignfile)

        #====================================
        # perform 1D query and check the distance threshold
        anchr_chr_pair = str(anchor_chr) + '|' + str(anchor_chr)
        mapped_align = align_tb.query(anchr_chr_pair, anchor_start, anchor_end)
        # # the output coordinates are basically the start of the mapped location
        # # so the distance threshold should be computed with respect to the "anchor_start"
        # for x in mapped_align:
        #     if (x[1] == x[3]) and (abs(int(x[4]) - anchor_start) >= distLowThres) and (abs(int(x[4]) - anchor_start) <= distUpThres):
        #         # Note: the bin id can be negative depending on the position of the mapped chromosome
        #         bin_id = (int(x[4]) - anchor_start) / bin_size
        #         if bin_id not in bins:
        #             bins.setdefault(bin_id, 0)
        #         bins[bin_id] += 1
        #         totread += 1
        # #====================================

        # #===========================================================
        # # the 2D query function in the pypairix package has the following syntax
        # # tb.query2D(chrom, start, end, chrom2, start2, end2)
        # # so we first decide a boundary region for searching the cis interactions
        # # in both +ve and -ve offset directions

        # # case 1 - when the given peak segment belongs to the left hand side 
        # # of the pairix file, the mapped portion will be on the right hand side
        # # and the coordinate of the mapped segment will be higher than the peak segment
        # # set the chromosome coordinates accordingly

        # start1_limit = anchor_start_bin + distLowThres
        # end1_limit = anchor_start_bin + bin_size + distUpThres

        # # map the current anchor region with the input alignment file
        # # subject to the region boundary decided
        # mapped_align = align_tb.query2D(anchor_chr, anchor_start, anchor_end, anchor_chr, start1_limit, end1_limit)

        for x in mapped_align:
            # x[1]: chr1, x[2]: chr1 position (start)
            # x[3]: chr2, x[4]: chr2 position (start)
            # check cis interactions and also the distance constraints

            # with respect to the binned interval, starting location of that bin
            mapped_bin_start = (int((int(x[4]) * 1.0) / bin_size)) * bin_size

            # key corresponding to the chromosome binned interval (python dictionary)
            # of the mapped segment
            curr_key = (anchor_chr, mapped_bin_start, (mapped_bin_start + bin_size))

            # check if the mapped region falls within the chromosome region
            if curr_key in SegmentDict:
                # Note: this condition is commented - sourya
                # The objective is to get a generic interaction file at first
                # after that, we will filter the interactions according to the distance threshold requirement
                
                # check cis interactions and also the distance constraints
                if (x[1] == x[3]):  # and (abs(mapped_bin_start - anchor_start_bin) >= distLowThres) and (abs(mapped_bin_start - anchor_start_bin) <= distUpThres):

                    # the mapped read is within the distance threshold
                    # methods 3 or 0 consider the other end of the read as either peak or non peak
                    # method 1 considers the other end of the read as peak
                    # method 2 considers the other end of the read as non peak
                    if (Method == 3) or (Method == 0) or ((Method == 1) and (SegmentDict[curr_key]._CheckIfPeak() == True)) or ((Method == 2) and (SegmentDict[curr_key]._CheckIfPeak() == False)):
                        bin_id = (mapped_bin_start - anchor_start_bin) / bin_size
                        if bin_id not in bins:
                            bins.setdefault(bin_id, 0)
                        bins[bin_id] += 1
                        totread += 1   

        # # case 2 - when the given peak segment belongs to the right hand side 
        # # of the pairix file, the mapped portion will be on the left hand side
        # # and the coordinate of the mapped segment will be lower than the peak segment
        # # set the chromosome coordinates accordingly
        # start2_limit = max(1, anchor_start_bin - distUpThres)
        # end2_limit = max(1, anchor_start_bin + bin_size - distLowThres)

        # # map the current anchor region with the input alignment file
        # # subject to the region boundary decided
        # mapped_align = align_tb.query2D(anchor_chr, start2_limit, end2_limit, anchor_chr, anchor_start, anchor_end)

        # for x in mapped_align:
        #     # x[1]: chr1, x[2]: chr1 position (start)
        #     # x[3]: chr2, x[4]: chr2 position (start)

        #     # with respect to the binned interval, starting location of that bin
        #     mapped_bin_start = (int((int(x[2]) * 1.0) / bin_size)) * bin_size

        #     # key corresponding to the chromosome binned interval (python dictionary)
        #     # of the mapped segment
        #     curr_key = (anchor_chr, mapped_bin_start, (mapped_bin_start + bin_size))

        #     # check if the mapped region falls within the chromosome region
        #     if curr_key in SegmentDict:
        #         # check cis interactions and also the distance constraints
        #         if (x[1] == x[3]) and (abs(mapped_bin_start - anchor_start_bin) >= distLowThres) and (abs(mapped_bin_start - anchor_start_bin) <= distUpThres):
        #             # the mapped read is within the distance threshold
        #             # methods 3 or 0 consider the other end of the read as either peak or non peak
        #             # method 1 considers the other end of the read as peak
        #             # method 2 considers the other end of the read as non peak

        #             # Note: the bin id can be negative depending on the position of the mapped chromosome
        #             if (Method == 3) or (Method == 0) or ((Method == 1) and (SegmentDict[curr_key]._CheckIfPeak() == True)) or ((Method == 2) and (SegmentDict[curr_key]._CheckIfPeak() == False)):
        #                 bin_id = (mapped_bin_start - anchor_start_bin) / bin_size
        #                 if bin_id not in bins:
        #                     bins.setdefault(bin_id, 0)
        #                 bins[bin_id] += 1
        #                 totread += 1 
    #=========================================================== 
    # this structure stores the results
    res = []

    for bin_id in bins:         
        # 'cis' interaction) format between any pair of segments
        mapped_start_bin = anchor_start_bin + bin_size * bin_id
        res.append([anchor_chr, anchor_start_bin, anchor_start_bin + bin_size, anchor_chr, mapped_start_bin, mapped_start_bin + bin_size, bins[bin_id]])
    
    if (input_file_fmt == 1):
        # close the bam alignment file
        bam_input_file.close()

    # return the list of interactions
    return res

#===============================================
if __name__ == "__main__":
    main()