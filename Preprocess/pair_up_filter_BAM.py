#!/usr/bin/env python
import optparse
import re
import sys
import os.path
from math import floor

def floored_percentage(val, digits):
    val *= 10 ** (digits + 2)
    return '{1:.{0}f}%'.format(digits, floor(val) / 10 ** digits)

class sam_alignment:
    """A class for sam alignment """
    def __init__(self, str):
        if len(str):
            elems = str.split()
            self.name = elems[0]
            self.flag = int(elems[1])
            self.chrom1 = elems[2]
            self.pos1 = int(elems[3])
            self.score = int(elems[4])
            self.cigar = elems[5]
            self.chrom2 = elems[6]            
            self.pos2 = int(elems[7])
            self.insert_size = int(elems[8])
            self.seq = elems[9]
            self.qual = elems[10]
            self.unknown1 = elems[11]
        else:
            exit('error: sam_alignment input is empty;')
    
    def aln_str(self):
        return '\t'.join([self.name,  str(self.flag), self.chrom1,
               str(self.pos1), str(self.score), str(self.cigar),
               self.chrom2, str(self.pos2), str(self.insert_size), 
               self.seq, self.qual, self.unknown1])

# this function checks the chromosome names
# returns TRUE / FALSE 
def CheckChromosome(align1_chr1, align1_chr2, align2_chr1, align2_chr2, Chr_List):
    if (align1_chr1 != "=") and (align1_chr1.lower() not in Chr_List):
        return False
    if (align1_chr2 != "=") and (align1_chr2.lower() not in Chr_List):
        return False
    if (align2_chr1 != "=") and (align2_chr1.lower() not in Chr_List):
        return False
    if (align2_chr2 != "=") and (align2_chr2.lower() not in Chr_List):
        return False
    return True


def main():
    """ main function """
    TOTAL_NUM_READS = 0
    PAIR_NUM_READS = 0
    TRANS_NUM_READS = 0
    CIS_NUM_READS = 0
    CIS_LONG_NUM_READS = 0

    parser = optparse.OptionParser(description='pair up HiC two mates and only keep uniquely mapped cis pairs.', usage='%prog [-h] [-f SAM1] [-r SAM2] [-q MAPQ_THR]')
    
    parser.add_option('-f', dest="SAM1", help='Paired end alignment file in .sam format.')
    parser.add_option('-q', dest="MAPQ_THR", help='Threshold of mapping quality for individual alignments')
    
    options, remainder = parser.parse_args()
    if options.SAM1:
        fin_r1_name = options.SAM1
    else:
        parser.print_help()
        exit('error: missing -f arguments;')

    # mapping quality threshold for individual alignments
    MAPQ_THR = int(options.MAPQ_THR)

    if not os.path.exists(fin_r1_name): exit('error: file %s not exist;' % fin_r1_name)

    # define a chromosome list within which all the chromosomes are considered
    Chr_List = []
    for i in range(1, 23):
        curr_chr = "chr" + str(i)
        Chr_List.append(curr_chr)
    Chr_List.append("chrx")
    Chr_List.append("chry")
    Chr_List.append("=")

    fin1 = open(fin_r1_name, "r")

    while True:
        line1 = fin1.readline()
        if (line1 == ""):
            break
        if (line1[0] == "@"):
            # header line
            print str(line1.rstrip('\n'))
        else:
            # line1 is not a header line
            # so get the next line also
            line2 = fin1.readline()

            # form the alignment structure using Pysam
            aln1 = sam_alignment(line1)
            aln2 = sam_alignment(line2)

            # sourya - this is an important condition
            # this condition checks if both the read pairs (coming from the two input files)
            # have the same name - otherwise, the code stops
            if aln1.name != aln2.name:
                if 1:
                    print '*** error in name: aln1.name ', aln1.name, '  aln2.name: ', aln2.name
                exit('error: hicmap_pair_up R1 and R2 read names not match;') 

            TOTAL_NUM_READS += 1

            # sourya - previously the alignment quality threshold was put as 10    
            # now it has been made dynamic, and generally its value is 30     
            if(aln1.score >= MAPQ_THR and aln2.score >= MAPQ_THR):
                if (CheckChromosome(aln1.chrom1, aln1.chrom2, aln2.chrom1, aln2.chrom2, Chr_List) == True):
                    PAIR_NUM_READS += 1
                    if aln1.chrom1 !=  aln2.chrom1: 
                        TRANS_NUM_READS += 1
                    else:
                        CIS_NUM_READS += 1
                        if aln1.flag == aln2.flag:
                            aln1.flag = 99
                            aln2.flag = 147
                            aln1.chrom2 = "="
                            aln2.chrom2 = "="
                            aln1.pos2 = aln2.pos1 
                            aln2.pos2 = aln1.pos1
                            aln1.insert_size =  aln2.pos1 - aln1.pos1
                            aln2.insert_size =  aln1.pos1 - aln2.pos1
                        else:
                            aln1.flag = 83
                            aln2.flag = 163
                            aln1.chrom2 = "="
                            aln2.chrom2 = "="
                            aln1.pos2 = aln2.pos1
                            aln2.pos2 = aln1.pos1
                            aln1.insert_size =  aln2.pos1 - aln1.pos1
                            aln2.insert_size =  aln1.pos1 - aln2.pos1
                        try:
                            print aln1.aln_str()
                            print aln2.aln_str()
                        except IOError:
                            try:
                                sys.stdout.close()
                            except IOError:
                                pass
                            try:
                                sys.stderr.close()
                            except IOError:
                                pass    
    
    fin1.close()
    print  >>sys.stderr, "number of totally sequenced pairs:", TOTAL_NUM_READS
    print  >>sys.stderr, "number of uniquely mapped pairs:", PAIR_NUM_READS, "( mappability =", floored_percentage(float(PAIR_NUM_READS)/TOTAL_NUM_READS, 2), ")"
    print  >>sys.stderr, "number of cis pairs:", CIS_NUM_READS, "( cis rate =", floored_percentage(float(CIS_NUM_READS)/PAIR_NUM_READS, 2), ")"
    
if __name__ == '__main__':
    main()