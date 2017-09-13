#!/usr/bin/env python
import sys
import re
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

def std_print(str):
    try:
        print str,
    except IOError:
        try:
            sys.stdout.close()
        except IOError:
            pass
        try:
            sys.stderr.close()
        except IOError:
            pass
    
    
def main():
    """ main function """
    pattern1=re.compile("@SQ\s+SN:chr[0-9a-zA-Z_]*\s+")
    pattern2=re.compile("@PG\s+ID:bwa\s+")
    NUM_LONG_PAIRS=0
    NUM_SHORT_PAIRS=0
    NUM_CIS_PAIRS=0    
    NUM_PAIRS=0

    DIST_CUTOFF=int(sys.argv[1])
    fin = sys.stdin
    read_prev = None
    while True:
        read_cur = fin.readline()
        if read_cur == "": break
        if (re.search(pattern1, read_cur) != None) or (re.search(pattern2, read_cur) != None):
            std_print(read_cur)
        else:
            if read_prev == None: 
                read_prev = read_cur
            else: 
                if read_cur.split()[0] == read_prev.split()[0]: # if they are a pair
                    NUM_PAIRS += 1
                    if read_cur.split()[2] == read_prev.split()[2]:
                        NUM_CIS_PAIRS += 1  
                        # only print long-range and short-range pairs   
                        if abs(int(read_prev.split()[8])) >= DIST_CUTOFF:
                            std_print(read_prev)
                            std_print(read_cur)
                            NUM_LONG_PAIRS += 1
                        if abs(int(read_prev.split()[8])) < 1000:
                            std_print(read_prev)
                            std_print(read_cur)
                            NUM_SHORT_PAIRS += 1
                    read_prev = read_cur
                else:
                    read_prev = read_cur
    fin.close()
    print  >>sys.stderr, "number of cis pairs close to cutter cites:", NUM_CIS_PAIRS, "( cutter site filter rate =", floored_percentage(float(NUM_CIS_PAIRS)/NUM_PAIRS, 2), ")"
    print  >>sys.stderr, "number of cis long-range pairs:", NUM_LONG_PAIRS, "( long-range rate =", floored_percentage(float(NUM_LONG_PAIRS)/NUM_CIS_PAIRS, 2), ")"
    print  >>sys.stderr, "number of cis short-range pairs:", NUM_SHORT_PAIRS, "( short-range rate =", floored_percentage(float(NUM_SHORT_PAIRS)/NUM_CIS_PAIRS, 2), ")"
    
if __name__ == '__main__':
    main()
