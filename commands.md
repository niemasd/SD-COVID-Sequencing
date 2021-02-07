# General Tips/Resources
* This [cookbook](https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb) might be useful
* The [iVar Manual](https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md) might also be useful
* See [`pipeline.sh`](pipeline.sh) for an implemented Bash script pipeline
    * To run it on all samples in the current folder: `parallel --jobs THREADS ./pipeline.sh {} ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)`

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
minimap2 -t THREADS -a -x sr REFERENCE_GENOME.FAS.MMI PRIMERS.FAS | samtools view -b -F 4 > PRIMERS.BAM
bedtools bamtobed -i PRIMERS.BAM > PRIMERS.BED
```
* I used to call `minimap2` with `-x map-ont`, but Karthik said I should use `-x sr`

# Step 1: Map Reads and Sort
* **Input:** FASTQ (or FASTQ.GZ) file(s) (`X.fastq` or `X.fastq.gz`)
* **Output:** Sorted Untrimmed BAM (`X.sorted.bam`)

## Individual Command
```bash
minimap2 -t THREADS -a -x sr ../ref/NC_045512.2.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools sort --threads THREADS -o SORTED.BAM
```
* * I used to call `minimap2` with `-x map-ont`, but Karthik said I should use `-x sr`

## Batch Command (64 threads per Minimap2, 1 at a time)
```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( minimap2 -t 64 -a -x sr ../ref/NC_045512.2.fas.mmi $s*.fastq.gz | samtools sort --threads 64 -o $s.sorted.bam ) ; } 2> $s.log.1.map.log ; done
```

## Batch Command (1 thread per Minimap2, 64 at a time)
```bash
parallel --jobs 64 "{" time "(" minimap2 -t 1 -a -x sr ../ref/NC_045512.2.fas.mmi {}*.fastq.gz "|" samtools sort --threads 1 -o {}.sorted.bam  ")" ";" "}" "2>" {}.log.1.map.log ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

## Speed Up by Truncating at Number of Mapped Reads
To speed up the entire workflow, we could theoretically just stop mapping reads once we hit a desired number of mapped reads (e.g. 250,000). If the output of Minimap2 is piped to `samtools view -F 4` to filter only for mapped reads, and then the SAM output of that is piped to `head -NUM_READS` to truncate at a certain number of reads, `head` should close the pipe (thus killing Minimap2 early), and downstream analyses will have much smaller number of reads to deal with.

For a single file, I could do the following:

```bash
minimap2 -t THREADS -a -x sr ../ref/NC_045512.2.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools view -h -F 4 | head -NUM_READS_PLUS_3 | samtools sort --threads THREADS -o SORTED.BAM
```
* You need to add 3 to the desired number of reads because the `-h` flag of `samtools view` outputs the header (which is 3 lines long)

For a batch command, I could do the following (I've hardcoded 64 threads and 250000 mapped reads):

```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( minimap2 -t 64 -a -x sr ../ref/NC_045512.2.fas.mmi $s*.fastq.gz | samtools view -h -F 4 | head -250003 | samtools sort --threads 64 -o $s.sorted.bam ) ; } 2> $s.log.1.map.log ; done
```

# Step 2: Trim Sorted BAM (resulting in unsorted trimmed BAM)
* **Input:** Sorted Untrimmed BAM (`X.sorted.bam`)
* **Output:** Unsorted Trimmed BAM (`X.trimmed.bam`)

## Individual Command
```bash
ivar trim -x 5 -e -i SORTED.BAM -b PRIMERS.bed -p TRIMMED_PREFIX
```
* `-x 5` is a new flag Karthik told me to add (requires iVar 1.3.1 or higher); something about extra bases added to the beginning/end of primers?
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
parallel --jobs 64 "{" time "(" ivar trim -x 5 -e -i {}.sorted.bam -b ../primers/swift/sarscov2_v2_primers.bed -p {}.trimmed ")" ";" "}" ">" {}.log.2.trim.log "2>&1" ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

## Using Swift's [`primerclip`](https://github.com/swiftbiosciences/primerclip)
At the time of writing this comment (2021-02-04), iVar Trim has a bug that is being fixed. As a quick-and-dirty alternative, Kristian Andersen recommended using Swift's [`primerclip`](https://github.com/swiftbiosciences/primerclip) tool for sequences generated using Swift primers. The usage from the README is as follows:

```bash
primerclip masterfile.txt alignmentfile.sam outputfilename.sam
```

Note that the SAMs need to be **name-sorted**, NOT coordinate-sorted. Also, it doesn't seem to natively support BAM files, but I should be able to use a FIFO to do so using `<()` for input and `>()` for output:

```bash
primerclip PRIMERS.txt <(samtools sort -n SORTED.bam | samtools view -h) >(samtools view -S -b > TRIMMED.bam)
```

The above SHOULD have worked, but `primerclip` is a bit janky and complains that there's no header... Might give up and wait for iVar Trim to get fixed.

# Step 3: Sort Trimmed BAM
* **Input:** Unsorted Trimmed BAM (`X.trimmed.bam`)
* **Output:** Sorted Trimmed BAM (`X.trimmed.sorted.bam`)

## Individual Command
```bash
samtools sort --threads THREADS -o TRIMMED_SORTED.BAM TRIMMED_PREFIX.BAM && rm TRIMMED_PREFIX.BAM
```
* I'm deleting the unsorted trimmed BAM to save space (the sorted trimmed BAM has all the same info + it's sorted)

## Batch Command (64 threads per samtools, 1 at a time)
```bash
for s in $(ls *.trimmed.bam | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq); do { time ( samtools sort --threads 64 -o $s.trimmed.sorted.bam $s.trimmed.bam && rm $s.trimmed.bam ) ; } 2> $s.log.3.sorttrimmed.log ; done
```

## Batch Command (1 thread per samtools, 64 at a time)
```bash
parallel --jobs 64 "{" time "(" samtools sort --threads 1 -o {}.trimmed.sorted.bam {}.trimmed.bam "&&" rm {}.trimmed.bam  ")" ";" "}" "2>" {}.log.3.sorttrimmed.log ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
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
    * However, given that I'll now be zipping everything at the end anyways, doesn't make sense to waste compute time compressing here as well

## Batch Command (64 threads)
```bash
parallel --jobs 64 "{" time "(" samtools mpileup -A -aa -d 0 -Q 0 --reference ../ref/NC_045512.2.fas {}.trimmed.sorted.bam ")" ";" "}" ">" {}.trimmed.sorted.pileup.txt "2>" {}.log.4.pileup.log ::: $(ls *.trimmed.sorted.bam | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
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
parallel --jobs 64 "{" time "(" cat {}.trimmed.sorted.pileup.txt.gz "|" ivar variants -r ../ref/NC_045512.2.fas -g ../ref/NC_045512.2.gff3 -p {}.trimmed.sorted.pileup.variants.tsv -m 10 ")" ";" "}" "2>" {}.log.5.variants.log ::: $(ls *.trimmed.sorted.pileup.txt.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
```

# Step 6: Call Consensus Sequence from Pile-Up
* **Input:** Pile-up (`X.pileup.txt` or `X.pileup.txt.gz`)
* **Ouptut:** Consensus Sequence (`X.consensus.fas` or `X.consensus.fas.gz`)

## Individual Command
```bash
cat PILEUP.TXT | ivar consensus -p CONSENSUS.FAS -m 10 -n N -t 0.5
```
* If pile-up files are gzipped, use `zcat PILEUP.TXT.GZ` instead of `cat PILEUP.TXT`
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
parallel --jobs 64 "{" time "(" cat {}.trimmed.sorted.pileup.txt.gz "|" ivar consensus -p {}.trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ")" ";" "}" ">" {}.log.6.consensus.log "2>&1" ::: $(ls *.trimmed.sorted.pileup.txt.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
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
parallel --jobs 64 "{" time "(" samtools depth -d 0 -Q 0 -q 0 -aa {}.trimmed.sorted.bam ")" ";" "}" ">" {}.trimmed.sorted.depth.txt "2>" {}.log.7.depth.log ::: $(ls *.trimmed.sorted.bam | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq)
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
parallel --jobs 64 "{" time "(" qualimap bamqc -bam {1}.{2}.bam -nt 1 --java-mem-size=4G -outdir {1}.{2}.stats "&&" tar c {1}.{2}.stats "|" pigz -9 -p 1 ">" {1}.{2}.stats.tar.gz "&&" rm -rf {1}.{2}.stats ")" ";" "}" ">" {1}.log.8.qualimap.{2}.log "2>&1" ::: $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq) ::: sorted trimmed.sorted
```

# When All Samples Finish
There are some things I do when all samples finish running, to archive as well as to compute summary stats.

## Archiving Each Sample's Files
It seems as though Google Shared Drives have a [400,000 file limit](https://support.google.com/a/answer/7338880?hl=en), and in general, there are tons of tiny files (which result in slow file transfer), so it makes sense to archive the files from each sample into a single file (e.g. zip). Now that I'm zipping things overall, I no longer need to gzip the pile-up files (because they'll get compressed when storing anyways).

### Individual Command
```bash
zip -9 SAMPLE.zip SAMPLE*
```

### Batch Command
```bash
for s in $(ls *.fastq.gz | sed 's/_R[12]_/./g' | cut -d'.' -f1 | sort | uniq) ; do zip -9 $s.zip $s* ; done
```

## Consolidating Qualimap Stats
I also wrote a [script](https://github.com/niemasd/tools/blob/master/qualimap_targz_to_TSV.py) to convert all the `X.trimmed.sorted.stats.tar.gz` files into a single TSV with all of the summary statistics.

### Command
```bash
qualimap_targz_to_TSV.py *.stats.tar.gz > output/qualimap.tsv
```

## Consolidating Consensus Sequences, Variant Calls, and Depths
People will likely mainly just be interested in consensus sequences, variant calls, and depths, so I zip them separately as well.

### Command
```bash
zip -9 consensus.zip *.consensus.fa && mv consensus.zip output/
zip -9 variants.zip *.variants.tsv && mv variants.zip output/
zip -9 depth.zip *.depth.txt && mv depth.zip output/
```

## Mapping Depth Distributions
It is useful to visualize the distributions of per-site mapping depth across the samples. I wrote a [script](https://github.com/niemasd/tools/blob/master/samtools_depth_violinplot.py) to generate a violin plot across all samples as well as a [script](https://github.com/niemasd/tools/blob/master/samtools_depth_lineplot.py) to generate line plots across all samples.

### Command
```bash
samtools_depth_violinplot.py *.depth.txt && mv depth_violin.pdf output/
samtools_depth_lineplot.py *.depth.txt && mv depth_lineplot.pdf output/
```
