#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "1.0"

import sys
from pysam import AlignmentFile

if __name__ == "__main__":

    # genome name
    samfile = AlignmentFile('test.bam')
    sam_out = AlignmentFile('test_out.sam', 'w', header=samfile.header)

    # resect coord
    coord = 2694 # 1-based
    glen = 5386

    for read in samfile.fetch():
        # advance to the resect coordinate
        read.reference_start += (coord-1) # so 1 becomes 2694
        read.next_reference_start += (coord-1)

        # if we've gone over the origin
        if read.reference_start > glen:
            read.reference_start -= glen
        if read.next_reference_start > glen:
            read.next_reference_start -= glen

        # write the read
        sam_out.write(read)
    samfile.close()
    sam_out.close()
