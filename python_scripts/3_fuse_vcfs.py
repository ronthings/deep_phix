#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "1.0"

import sys
from pysam import VariantFile

if __name__ == "__main__":

    # genome name
    vcf_name = sys.argv[1].split('_1.vcf')[0]
    vcf_primary = VariantFile(vcf_name+'_1.vcf')  # auto-detect input format
    vcf_secondary = VariantFile(vcf_name+'_2_reindexed.vcf')  # auto-detect input format
    vcf_out = VariantFile(vcf_name+'_merged.vcf', 'w', header=vcf_primary.header)

    # identify all sites
    first_coords = [rec.pos for rec in vcf_primary.fetch()]
    second_coords = [rec.pos for rec in vcf_secondary.fetch()]
    all_coords = set(first_coords + second_coords)
    all_coords = sorted(list(all_coords))
    #print(all_coords)

    # main loop
    for site in all_coords:
        match_xlist = [rec for rec in vcf_primary.fetch() if rec.pos==site]
        match_ylist = [rec for rec in vcf_secondary.fetch() if rec.pos==site]

        if len(match_xlist) == 0: # no match
            recx_coverage = -1 # any positive number is larger than this
        else:
            assert len(match_xlist) == 1
            recx = match_xlist[0]
            recx_coverage = recx.info["DP"]

        if len(match_ylist) == 0: # no match
            recy_coverage = -1 # any positive number is larger than this
        else:
            assert len(match_ylist) == 1
            recy = match_ylist[0]
            recy_coverage = recy.info["DP"]

        if recy_coverage > recx_coverage:
            vcf_out.write(recy)
        else:
            vcf_out.write(recx)
    vcf_primary.close()
    vcf_secondary.close()
    vcf_out.close()
