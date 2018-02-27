# VCF-Codon-Table
The vcf_parser.py script requires python3 and takes the following inputs:
* a FASTA genome file,
* an NCBI-style feature table,
* a FreeBayes output VCF file.

The script generates a filtered list of mutations in tabular form, giving the codon changes that have occurred. Please note this script is not yet production ready.

Sample usage:
```
./vcf_parser --help
./vcf_parser.py AF176034.fasta phix_coord.txt 40snps.vcf output.tsv
```
Another script (gbk_vcf_parser.py) is in development. The intention is to use Genbank inputs via [biopython](http://biopython.org/DIST/docs/tutorial/Tutorial.html) and [cyvcf2](http://brentp.github.io/cyvcf2/). Please feel free to fix this script - it isn't built yet...
