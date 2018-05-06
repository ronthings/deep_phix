## Procedure
# Now let's make the directory for VCF files to go in
# mkdir VCF

# Now the key FreeBayes command:
freebayes --fasta-reference FASTA/reference.fasta --pooled-continuous \
--min-alternate-fraction 0.01 --min-alternate-count 1 \
--min-mapping-quality 20 --min-base-quality 30 \
--no-indels --no-mnps --no-complex MAP/A_phixmapped.bam > VCF/A_phix.vcf

# I notice that FreeBayes itself left aligns indels by default.
# Why do some folks use pooled discrete with ploidy=1 ?

# Now it's time to filter the FreeBayes output (better to do this as a discrete step and we can use a program from vcflib):
# conda install vcflib # if not already installed
vcffilter -f "SRP > 20" -f "SAP > 20" -f "EPP > 20" -f "QUAL > 30" -f "DP > 30" VCF/A_phix.vcf > VCF/A_phix_filtered.vcf

# This filter is quite stringent. For example 1460 was rejected for placement bias but it looks pretty solid on IGV.
# Paradoxically, the false negative rate may be increased by high coverage...

# Time to tabulate our data so it is human-readable.
# Here's a simple approach (again using a vcflib program):
vcf2tsv VCF/A_phix.vcf > VCF/A_phix.tsv

# You can render the output file off-the-bat in Excel or in Atom if you have this excellent package installed:
https://atom.io/packages/tablr

# For a table with genome-context and biological effect, we can use my vcf-codon-table script from GitHub.
# Pre-requisites:
## We need to be in an environment with python3 (maybe move to your root conda environment).
## We need to put the phix_coord.txt information into the FASTA/ reference folder.

# Now we are ready:
./vcf_parser.py FASTA/phix_AF176034.fasta FASTA/phix_coord.txt VCF/A_phix.vcf VCF/A_table.tsv
