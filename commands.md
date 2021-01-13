# General Tips/Resources
* This [cookbook](https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb) might be useful
* The [iVar Manual](https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md) might also be useful

# Step 0: Preparation
## Index the Reference Genome
```bash
minimap2 -t THREADS -d REFERENCE_GENOME.FAS.MMI REFERENCE_GENOME.FAS
```

This analysis is using [NC045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta) as the reference genome, so specifically:

```bash
minimap2 -t THREADS -d NC045512.fas.mmi NC045512.fas
```

## Generate Primer Coordinate BED
Example from the [cookboox](https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb) was the following:

```bash
bwa mem -k 5 -T 16 db/PRV.fa db/zika_primers.fa | samtools view -b -F 4 > db/zika_primers.bam
bedtools bamtobed -i db/zika_primers.bam > db/zika_primers.bed
```

I wonder if Minimap2 would be fine instead of BWA-MEM? So I propose the following:

```bash
minimap2 -t THREADS -a -x map-ont REFERENCE_GENOME.FAS.MMI PRIMERS.FAS | samtools view -b -F 4 > PRIMERS.BAM
bedtools bamtobed -i PRIMERS.BAM > PRIMERS.BED
```

# Step 1: Map Reads and Sort
* **Input:** FASTQ (or FASTQ.GZ) file(s) (`X.fastq` or `X.fastq.gz`)
* **Output:** Sorted Untrimmed BAM (`X.sorted.bam`)

## Individual Command
```bash
minimap2 -t THREADS -a -x map-ont ../ref/NC_045512.2.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools sort --threads THREADS -o SORTED.BAM
```

## Batch Command (64 threads)
```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( minimap2 -t 64 -a -x map-ont ../ref/NC_045512.2.fas.mmi $s*.fastq.gz | samtools sort --threads 64 -o $s.sorted.bam ) ; } 2> $s.log.1.map.log ; done
```

# Step 2: Trim Sorted BAM (resulting in unsorted trimmed BAM)
* **Input:** Sorted Untrimmed BAM (`X.sorted.bam`)
* **Output:** Unsorted Trimmed BAM (`X.trimmed.bam`)

## Individual Command
```bash
ivar trim -e -i SORTED.BAM -b PRIMERS.bed -p TRIMMED_PREFIX
```
* `TRIMMED_PREFIX` is the output file minus the `.bam` extension (e.g. `-p sample_name.trimmed` would result in `sample_name.trimmed.bam`)
* `-e` is "Include reads with no primers"
    * In the `Snakefile`, `trim_reads_illumina` has it, but `trim_reads_ont` doesn't have it
        * `trim_reads_illumina` has a comment `# Add -e if nextera used`
    * I don't see why not to include reads with no primers, so I'll keep `-e` unless someone tells me otherwise
* `-q` is the minimum quality score
    * In the `Snakefile`, `trim_reads_illumina` uses the default (whichi s 20), but `trim_reads_ont` uses it's 5
* The `Snakefile` and cookbook both say to sort the `BAM` after, but if it was sorted before, shouldn't it be sorted after trimming...?
    * Maybe the length trimmed can vary, so if read *x* started earlier than read *y* but more was cut off the beginning, it should actually come after read *y*?
    * If so, maybe I'll output to a temporary file (e.g. `-p sample_name.trimmed.bam`) and then have the final output be a different file (e.g. `sample_name.trimmed.sorted.bam`)
        * If I do this, I should probably just delete the unsorted trimmed BAM to save space (the sorted trimmed BAM has all the same info + it's sorted)
* Also, why does `ivar trim` need a sorted BAM?
    * Sorting doesn't take much extra time (I just pipe `minimap2`'s output to `samtools sort`), so not a big deal, but unclear why that pipe is a step

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" ivar trim -e -i {}.sorted.bam -b ../primers/swift/sarscov2_v2_primers.bed -p {}.trimmed ")" ";" "}" ">" {}.log.2.trim.log "2>&1" ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

# Step 3: Sort Trimmed BAM
* **Input:** Unsorted Trimmed BAM (`X.trimmed.bam`)
* **Output:** Sorted Trimmed BAM (`X.trimmed.sorted.bam`)

## Individual Command
```bash
samtools sort --threads THREADS -o TRIMMED_SORTED.BAM TRIMMED_PREFIX.BAM && rm TRIMMED_PREFIX.BAM
```
* I'm deleting the unsorted trimmed BAM to save space (the sorted trimmed BAM has all the same info + it's sorted)

## Batch Command (64 threads)
```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( samtools sort --threads 64 -o $s.trimmed.sorted.bam $s.trimmed.bam && rm $s.trimmed.bam ) ; } 2> $s.log.3.sorttrimmed.log ; done
```

# Step 4: Generate Pile-Up from Trimmed Sorted BAM
* **Input:** Sorted Trimmed BAM (`X.trimmed.sorted.bam`)
* **Output:** Pile-up (`X.pileup.txt` or `X.pileup.txt.gz`)

## Individual Command
```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference REFERENCE.FAS TRIMMED_SORTED.BAM > PILEUP.TXT
```
* I think everything that depends on the pile-up will stream it from standard input
    * Thus, to save space, gzip compression makes sense (e.g. `pigz -9 -p THREADS`)
* Because there's such huge variance in how long the pile-up generation takes across files, it makes sense to use a 64 core machine and do all the pile-ups in parallel, and once most of them are done, compress some of the finished ones using multithreaded `pigz` on the idling cores

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" samtools mpileup -A -aa -d 0 -Q 0 --reference ../ref/NC_045512.2.fas {}.trimmed.sorted.bam ")" ";" "}" ">" {}.trimmed.sorted.pileup.txt "2>" {}.log.4.pileup.log ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```
* There was quite a bit of variance between the runtimes, so rather than have the finished cores idle, I ran `pigz -9 -p THREADS` to compress all of the pile-up files that had already completed running as the last ones were finishing up

# Step 5: Call Variants from Pile-Up
* **Input:** Pile-up (`X.pileup.txt` or `X.pileup.txt.gz`)
* **Output:** Variants (`X.variants.tsv`)

## Individual Command
```bash
cat PILEUP.TXT | ivar variants -r REFERENCE.FAS -g REFERENCE.GFF -p VARIANTS.TSV -m 10
```
* If pile-up files are gzipped, use `zcat PILEUP.TXT.GZ` instead of `cat PILEUP.TXT`
* `-m` is the minimum read depth to call variants
    * The default is 0, but the only place in the `Snakefile` where `ivar variants` is called uses 10
* Potentially gzip the output TSV files to save space?
    * They're tiny (less than 1 MB each), so not worth it

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" zcat {}.trimmed.sorted.pileup.txt.gz "|" ivar variants -r ../ref/NC_045512.2.fas -g ../ref/NC_045512.2.gff3 -p {}.trimmed.sorted.pileup.variants.tsv -m 10 ")" ";" "}" "2>" {}.log.5.variants.log ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

# Step 6: Call Consensus Sequence from Pile-Up
* **Input:** Pile-up (`X.pileup.txt` or `X.pileup.txt.gz`)
* **Ouptut:** Consensus Sequence (`X.consensus.fas` or `X.consensus.fas.gz`)

## Individual Command
```bash
cat PILEUP.TXT | ivar consensus -p CONSENSUS.FAS -m 10 -n N -t 0.5
```
* `-m` is the minimum depth to call consensus
    * Default is 10, so not sure why they explicitly specify `-m 10` in the `Snakefile`
* `-n` is the ambiguous character (`N` or `-`) to print in sites that have lower than `-m` coverage
    * Default is `N`, so not sure why they explicitly specify `-n N` in the `Snakefile`
* `-t` is the minimum frequency threshold to call the consensus
    * In the `Snakefile`, `call_consensus_illumina` uses 0.5, but `call_consensus_ont` uses the default (which is 0)
* `-q` is the minimum quality score threshold to count a given base
    * In the `Snakefile`, `call_consensus_illumina` uses the default (which is 20), but `call_consensus_ont` uses 5
* This outputs a log file to standard output
    * When I run it in batch, I should redirect standard error to standard output (`2>&1`) to get `time` output and the `ivar consensus` log into the same log file

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" zcat {}.trimmed.sorted.pileup.txt.gz "|" ivar consensus -p {}.trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ")" ";" "}" ">" {}.log.6.consensus.log "2>&1" ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

# Step 7: Call Depth (supplemental summary stats)
* **Input:** Sorted Trimmed BAM (`X.trimmed.sorted.bam`)
* **Output:** Depth (`X.trimmed.sorted.depth.txt`)

## Individual Command
```bash
samtools depth -d 0 -Q 0 -q 0 -aa TRIMMED_SORTED.BAM > DEPTH.TXT
```
* Potentially gzip the output TXT files to save space?
    * They're tiny (less than 1 MB each), so not worth it

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" samtools depth -d 0 -Q 0 -q 0 -aa {}.trimmed.sorted.bam ")" ";" "}" ">" {}.trimmed.sorted.depth.txt "2>" {}.log.7.depth.log ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

# Step 8: Qualimap Report (supplemental summary stats)
* **Input:** Sorted BAM (`X.sorted.bam`) and Sorted Trimmed BAM (`X.trimmed.sorted.bam`)
* **Output:** Qualimap Reports (`X.sorted.stats.tar.gz` and `X.trimmed.sorted.stats.tar.gz`)

## Individual Command
```bash
qualimap bamqc -bam MAPPING.BAM -nt THREADS --java-mem-size=RAM -outdir STATS_DIR && tar c STATS_DIR | pigz -9 -p THREADS > STATS_DIR.tar.gz && rm -rf STATS_DIR
```

## Batch Command (64 threads)
```
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do for f in sorted trimmed.sorted ; do { time ( qualimap bamqc -bam $s.$f.bam -nt 64 --java-mem-size=4G -outdir $s.$f.stats && tar c $s.$f.stats | pigz -9 -p 64 > $s.$f.stats.tar.gz && rm -rf $s.$f.stats ) ; } > $s.log.8.qualimap.$f.log 2>&1 ; done ; done
```

## Consolidating Stats
I also wrote a [script](https://github.com/niemasd/tools/blob/master/qualimap_targz_to_TSV.py) to convert all the `X.trimmed.sorted.stats.tar.gz` files into a single TSV with all of the summary statistics:

```bash
qualimap_targz_to_TSV.py *.trimmed.sorted.stats.tar.gz > YYYY-MM-DD.trimmed.sorted.stats.tsv
```
