# Notes on FreeBayes and tabulation

## The key FreeBayes command
```
freebayes --fasta-reference FASTA/reference.fasta --pooled-continuous \
--min-alternate-fraction 0.01 --min-alternate-count 1 \
--min-mapping-quality 20 --min-base-quality 30 \
--no-indels --no-mnps --no-complex MAP/A_phixmapped.bam > VCF/A_phix.vcf
```

## Notes
* I notice that FreeBayes itself left aligns indels by default.
* Why do some folks use pooled discrete with ploidy=1 ?

## The key vcffilter command
```
conda install vcflib # if not already installed
vcffilter -f "SRP > 20" -f "SAP > 20" -f "EPP > 20" -f "QUAL > 30" -f "DP > 30" VCF/A_phix.vcf > VCF/A_phix_filtered.vcf
```

## Notes
* It's better to do this as a discrete step and we can use a program from vcflib
* This filter is quite stringent. For example 1460 was rejected for placement bias but it looks pretty solid on IGV.
* Paradoxically, the false negative rate may be increased by high coverage...

## Tabulation for human-readability (using a vcflib program):
```
vcf2tsv VCF/A_phix.vcf > VCF/A_phix.tsv
```

## Notes
* You can render the output file off-the-bat in Excel or in Atom if you have [this excellent package](https://atom.io/packages/tablr) installed.
* For a table with genome-context and biological effect, we can use my vcf-codon-table script.

## The key command for my vcf_parser script
```
./vcf_parser.py FASTA/phix_AF176034.fasta FASTA/phix_coord.txt VCF/A_phix.vcf VCF/A_table.tsv
```
