# Preparation of FASTA Reference Sequences

## Data Source (Bacterial Genomes)
Bacterial genome files are generated via Unicycler hybrid assemblies of Illumina (MiSeq) and ONT MinION data. Two oddities must be noted:

1. The _E. coli_ reference includes the ΦX174 control track. The full genome of ΦX174 starts at line 3 (after the single contig).
```
tail -n +3 assembly.fasta > NEB_phix.fasta
head -2 assembly.fasta > EC_reference.fasta
```
Please check on above (and check for identity using diff)...

2. The _S. enterica_ serovar Typhimurium reference contains a plasmid and is still fragmented (owing to fragmentation during ONT prep?)

## Data Source (Short Genomes)
Viral (ΦX174) and spike-in (pUC18) genomes are obtained from NCBI with accession numbers noted in the name of the file.

## Workflow for Splicing Genomes for Mapping Steps in DeepSeq Pipeline
First activate a conda environment with python 3.x:
```
prompt$ conda activate
```

### Resect Small Genomes

1. pUC18 (spike-in):
```
(base) prompt$ ./resect_genome.py pUC18_L09136.fasta
Genome resected. It now begins at position 1344
```

2. ΦX174 (viral):
```
(base) prompt$ ./resect_genome.py phix_AF176034.fasta
Genome resected. It now begins at position 2694
```

### Create Resected Bacterial References
1. Do this by skipping the ΦX174 genome and replacing with resected version
```
(base) prompt$ head -2 EC_reference.fasta > EC_reference_resected.fasta
(base) prompt$ cat phix_AF176034_resected.fasta >> EC_reference_resected.fasta
```

2. For the ST genome:
```
(base) prompt$ cat ST_reference_FIRSTPASS.fasta > ST_reference_FIRSTPASS_resected.fasta
(base) prompt$ cat phix_AF176034_resected.fasta >> ST_reference_FIRSTPASS_resected.fasta
```

3. Don't forget also to concatenate the standard ΦX174 to the ST reference:
```
(base) prompt$ cat phix_AF176034.fasta >> ST_reference_FIRSTPASS.fasta
```
