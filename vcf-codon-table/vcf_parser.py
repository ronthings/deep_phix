#! /usr/bin/env python
__author__ = "Ben Dickins"
__status__ = "Prototype"
__version__ = "0.1"

# pretty print function
def prettyprint(x, y, handle):
    if len(x) == 0:
        print("INT", end="\t", file=handle)
    elif len(x) == 1:
        print(x[0], end="\t", file=handle)
    elif len(x) == 2:
        if y == 0:
            print( "{}({})".format(x[0],x[1]), end="\t", file=handle )
        else:
            print( "{}/{}".format(x[0],x[1]), end="\t", file=handle )
    else:
        print( "{}/{}".format(x[0],x[1]), end="", file=handle )
        for item in x[2:]:
            handle.write('(' + item + ')')
        handle.write('\t')

# DNA/protein position reporter - N.B. all inputs and outputs are 1-based!
def withincheck(position, gene_start, gene_end, length):
    if gene_start > gene_end: # adjustments for origin breakers
        if position <= gene_end:
            position += length
        gene_end += length
    if gene_start <= gene_end: # now the main logic
        if position >= gene_start and position <= gene_end:
            dna_pos = position - gene_start + 1
            prot_pos = -(-dna_pos//3)
            return(dna_pos, prot_pos)
        else:
            return(False)
    else:
        raise UserWarning("Gene start/end conflict (even correcting for genome length).")
# note the // operator used for prot_pos rounds down when numbers are negative
# hijacked it here with double negative for rounding up

# main loops (also accessible if called as vcf_parser.main)
def main():
    # handle command line input from user
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_fa", type=str, help="FASTA genome file")
    parser.add_argument("feat_table", type=str, help="feature table file")
    parser.add_argument("vcf_file", type=str, help="FreeBayes VCF output")
    parser.add_argument("output", type=str, help="Output file")
    args = parser.parse_args()

    # read single fasta sequence file
    with open(args.genome_fa,'rU') as file:
        genome = ''
        for line in file:
            if not line:
                pass
            elif not line.startswith('>'):
                genome += line.rstrip()
        length = len(genome)
        print("Measured genome length is", length, "bases.")

    # read coordinates from feature table
    genecoords, subsequences = {}, {}
    with open(args.feat_table,'rU') as file:
        for line in file:
            if not line:
                pass
            elif "Product" not in line:
                line = line.split('\t')
                gen = line[0].strip()
                beg, ter = int(line[1].strip()), int(line[2].strip())
                genecoords[gen] = (beg, ter)
                # coordinates so far are 1-based, so seq slices must accommodate this
                if beg <= ter:
                    subseq = genome[ beg-1:ter ] # remember python slicing is exclusive of end
                elif beg > ter:
                    subseq = genome[ beg-1:length] + genome[ 0:ter ]
                subsequences[gen] = subseq

    # genetic code #11, but not annotating alternative start codons
    # reference: http://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html#SG11
    codons = {
        'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
        'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
        'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

    # types of amino acids - N.B. this is a minimal classification
    # could switch to this EBI classification: https://en.wikipedia.org/wiki/Conservative_mutation
    groups = {
        'G':'NP', 'A':'NP', 'V':'NP', 'L':'NP', 'I':'NP',
        'F':'NP', 'M':'NP', 'P':'NP', 'W':'NP', 'S':'PO',
        'T':'PO', 'Y':'PO', 'C':'PO', 'N':'PO', 'Q':'PO',
        'D':'AC', 'E':'AC', 'H':'BA', 'K':'BA', 'R':'BA',
        '*':'ST'}

    # ask the user to identify a character that describes overlapping genes in the same ORF
    orf_char = input("Enter character found in non-ARF overlapping genes (default: *): ") or "*"

    # read vcf file checking each (alternative base at each) position against all genes
    with open(args.vcf_file, 'rU') as file:
        outfile = open(args.output,'w')
        print("Site\tProtein(s)\tAmino acid(s)\tRadicality\tCoverage\tMAC\tMAF", file=outfile)
        for line in file:
            if not line:
                pass
            elif not line.startswith("#"):
                line = line.split("\t")
                refpos = int(line[1]) # very important to note this is 1-based!!
                refnuc, altnucs = line[3], line[4].split(',')
                field_names = line[8].split(':')
                dp_idx, ad_idx = field_names.index('DP'), field_names.index('AD')
                field_num = line[9].split(':')
                dp_num, ad_num = field_num[dp_idx], field_num[ad_idx].split(',')[1:] # 0th val is ref allele

                # we may have >1 alternative allele so we must loop through these + their # of reads
                assert len(altnucs) == len(ad_num) # these should be the same length
                for base, num in zip(altnucs, ad_num):
                    prot_list, change_aa, change_type = [], [], []
                    orf_check = 0 # will be passed to prettyprint as y (positional) argument
                    change_pos = refnuc + str(refpos) + base
                    for gen, coord in genecoords.items():
                        beg, ter = coord # tuple unpacking
                        within_gene = withincheck(refpos, beg, ter, length)
                        if within_gene:
                            prot_list.append(gen)
                            if orf_char in gen:
                                orf_check = 1
                            dna_pos, prot_pos = within_gene # tuple unpacking
                            dna_pos -= 1 # this is now 0-based!!
                            codon_position = dna_pos % 3 # codon: 0, 1, 2

                            subseq = subsequences[gen]

                            if codon_position == 0:
                                refcodon = subseq[dna_pos:dna_pos+3]
                                newcodon = base + refcodon[1:3]

                            elif codon_position == 1:
                                refcodon = subseq[dna_pos-1:dna_pos+2]
                                newcodon = refcodon[0] + base + refcodon[2]

                            else:
                                refcodon = subseq[dna_pos-2:dna_pos+1]
                                newcodon = refcodon[0:2] + base

                            oldaa = codons[refcodon]
                            newaa = codons[newcodon]
                            change_aa.append(oldaa + str(prot_pos) + newaa)
                            if oldaa == newaa:
                                change_type.append("SYN")
                            elif groups[oldaa] == groups[newaa]:
                                change_type.append('CON')
                            else:
                                change_type.append('RAD')

                    print(change_pos, end="\t", file=outfile)
                    prettyprint(prot_list, orf_check, outfile)
                    prettyprint(change_aa, orf_check, outfile)
                    prettyprint(change_type, orf_check, outfile)
                    print(dp_num, num, float(num)/float(dp_num), sep="\t", end="\n", file=outfile)
        outfile.close()

if __name__ == "__main__":
    main()
