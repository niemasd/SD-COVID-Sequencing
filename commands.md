# Map Reads and Sort
```bash
minimap2 -t 16 -a -x map-ont NC045512.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools sort --threads 16 -o SORTED.BAM
```

# Generate Pile-Up from Sorted BAM
```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference NC045512.fas SORTED.BAM > PILEUP.TXT
```
