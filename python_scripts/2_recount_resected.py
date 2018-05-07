#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "1.0"

import sys
from pysam import VariantFile

if __name__ == "__main__":

    # genome name
    vcf_name = sys.argv[1].split('.vcf')[0]
    vcf_in = VariantFile(vcf_name+'.vcf')  # auto-detect input format
    vcf_out = VariantFile(vcf_name+'_reindexed.vcf', 'w', header=vcf_in.header)

    # resect coord
    print(list((vcf_in.header.contigs)))
    coord = 2694 # 1-based
    glen = 5386

    # main loop
    for rec in vcf_in.fetch():
        rec.pos += (coord-1) # so 1 becomes 2694
        if rec.pos > glen:
            rec.pos -= glen # so 5387 becomes 1
        vcf_out.write(rec)
        #print(rec, file=vcf_out)
    vcf_in.close()
    vcf_out.close()
