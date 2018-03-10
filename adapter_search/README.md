# Command line for adapter like sequence search:

Search for adapter like sequences (writes to adapterlike.txt):
```
./adapterlike_detect.py AF176034.fasta
```

Summarise copy numbers to file:
```
tail -n+2 adapterlike.txt | cut -f1 | sort | uniq -c > adapter_counts.txt
```
