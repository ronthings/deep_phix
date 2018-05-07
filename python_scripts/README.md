# Example for command-line usage on defined sample = A

## What's there?
```
(phix) NHBLAP-MAC01:VCF bio3dickib$ ls
A_phix_1.vcf  A_phix_2.vcf  C1_phix_1.vcf C1_phix_2.vcf
```

## Re-count the second file
```
(phix) NHBLAP-MAC01:VCF bio3dickib$ ../python_scripts/2_recount_resected.py A_phix_2.vcf
['1', 'RESTART_2694_RESECTED_AF176034.1']
```

## Merge the VCFs (uses recounted file)
```
(phix) NHBLAP-MAC01:VCF bio3dickib$ ../python_scripts/3_fuse_vcfs.py A_phix_1.vcf
```
## Now take a look
```
(phix) NHBLAP-MAC01:VCF bio3dickib$ ls
A_phix_1.vcf           A_phix_2.vcf           A_phix_2_reindexed.vcf A_phix_merged.vcf      C1_phix_1.vcf          C1_phix_2.vcf
```

## Similar commands for another sample = C1
```
(phix) NHBLAP-MAC01:VCF bio3dickib$ ../python_scripts/2_recount_resected.py C1_phix_2.vcf
['1', 'RESTART_2694_RESECTED_AF176034.1']
(phix) NHBLAP-MAC01:VCF bio3dickib$ ../python_scripts/3_fuse_vcfs.py C1_phix_1.vcf
(phix) NHBLAP-MAC01:VCF bio3dickib$ ls
A_phix_1.vcf            A_phix_2.vcf            A_phix_2_reindexed.vcf  A_phix_merged.vcf       C1_phix_1.vcf           C1_phix_2.vcf           C1_phix_2_reindexed.vcf C1_phix_merged.vcf
```
