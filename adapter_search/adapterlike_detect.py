#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "1.0.0"

import re, sys

if __name__ == "__main__":
    adapter = "CTGTCTCTTATA"
    motifs = []
    for i in range(1,len(adapter)):
        motifs.append(adapter[0:i+1])
    print(motifs)

    # read the FASTA file into a variable (seq)
    with open(sys.argv[1], 'rU') as fasta:
        all_lines = fasta.readlines()
        seq = ''
        for line in all_lines:
            line = line.rstrip()
            if line.startswith('>'):
                pass
            else:
                seq += line

    # analyse for repeats
    with open('adapterlike.txt', 'w') as outfile:
        # header lines (X2)
        print("motif\tstart\tend", file=outfile)
        # loop through each motif
        for unit in motifs:
            regex = re.compile("({})".format(unit), re.IGNORECASE)
            discovery = [(m.group(), m.start()+1, m.end()) for m in re.finditer(regex, seq)]
            for tup in discovery:
                print("\t".join([str(t) for t in tup]), file = outfile)

# dead darling:
# pat = re.compile(r'(A){3,}|(C){3,}|(G){3,}|(T){3,}', re.IGNORECASE)
