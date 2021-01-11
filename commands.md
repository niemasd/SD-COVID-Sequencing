# Map Reads and Sort
## Individual Command
```bash
minimap2 -t 64 -a -x map-ont ../ref/NC045512.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools sort --threads 64 -o SORTED.BAM
```

## Batch Command
```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( minimap2 -t 64 -a -x map-ont ../ref/NC045512.fas.mmi $s*.fastq.gz | samtools sort --threads 64 -o $s.sorted.bam ) ; } 2> $s.log.1.map.log ; done
```

# Generate Pile-Up from Sorted BAM
## Individual Command
```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference NC045512.fas SORTED.BAM | pigz -9 -p 64 > PILEUP.TXT.GZ
```

## Batch Command
```bash
{ time parallel --jobs 64 samtools mpileup -A -aa -d 0 -Q 0 --reference ../ref/NC045512.fas {}.sorted.bam "|" pigz -9 -p 64 ">" {}.sorted.pileup.txt.gz ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq) ; } 2> pileup.time.log
```
