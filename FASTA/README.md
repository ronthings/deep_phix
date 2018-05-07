# Preparation of FASTA Reference Sequences

Full genome of phage starts at line 3 in _E. coli_ reference, e.g.
```
tail -n +3 ST_reference_FAKE.fasta
```

# Starting point
```
NHBLAP-MAC01:FASTA bio3dickib$ ls -1
EC_reference.fasta
ST_reference_FAKE.fasta
pUC18_L09136.fasta
phix_AF176034.fasta
pipeline.txt
resect_genome.py
```

# resect the spike-in genome
```
NHBLAP-MAC01:FASTA bio3dickib$ conda activate
(base) NHBLAP-MAC01:FASTA bio3dickib$ ./resect_genome.py pUC18_L09136.fasta
Genome resected. It now begins at position 1344
```

# resect the phiX genome
```
(base) NHBLAP-MAC01:FASTA bio3dickib$ ./resect_genome.py phix_AF176034.fasta
Genome resected. It now begins at position 2694
```

# create the resected EC reference by skipping the phiX genome and replacing with resected version
```
(base) NHBLAP-MAC01:FASTA bio3dickib$ head -2 EC_reference.fasta > EC_reference_resected.fasta
(base) NHBLAP-MAC01:FASTA bio3dickib$ cat phix_AF176034_resected.fasta >> EC_reference_resected.fasta
```
# do the same for the ST version
```
(base) NHBLAP-MAC01:FASTA bio3dickib$ cat ST_reference_FIRSTPASS.fasta > ST_reference_FIRSTPASS_resected.fasta
(base) NHBLAP-MAC01:FASTA bio3dickib$ cat phix_AF176034_resected.fasta >> ST_reference_FIRSTPASS_resected.fasta
```

# don't forget to concatenate to ST reference
```
(base) NHBLAP-MAC01:FASTA bio3dickib$ cat phix_AF176034.fasta >> ST_reference_FIRSTPASS.fasta
```

# now we have
```
(base) NHBLAP-MAC01:FASTA bio3dickib$ ls -1
EC_reference.fasta
EC_reference_resected.fasta
ST_reference_FIRSTPASS.fasta
ST_reference_FIRSTPASS_resected.fasta
pUC18_L09136.fasta
pUC18_L09136_resected.fasta
phix_AF176034.fasta
phix_AF176034_resected.fasta
pipeline.txt
ref_decoder.csv
resect_genome.py
```