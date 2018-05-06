#!/usr/bin/env python

__author__ = "Ben Dickins"
__version__ = "1.0"

import sys, math

if __name__ == "__main__":

    # genome name
    genome_name = sys.argv[1].split(".fa")[0]

    # read the FASTA file into a variable (genome)
    with open(sys.argv[1], 'rU') as fasta:
        all_lines = fasta.readlines()
        genome = ""
        for line in all_lines:
            line = line.rstrip()
            if line.startswith('>'):
                chrom = line[1:]
            else:
                genome += line

    # resect the genome
    with open(genome_name + '_resected.fasta', 'w') as outfile:
        length = math.ceil(len(genome)/2)
        newseq = genome[length:] + genome[:length]

        # created new header - note the +1 because python starts at 0
        print(">RESTART_{}_RESECTED_{}".format(length+1, chrom), file=outfile)
        # now print the sequence lines
        for i in range(0,len(newseq),70):
            print(''.join(newseq[i:i+70]), file=outfile)
        # Let's report what we've done to stdout
        print("Genome resected. It now begins at position", length+1)
