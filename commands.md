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

## Batch Command
```bash
TODO
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

## Batch Command
```bash
TODO
```

# Step 7: Call Depth (supplemental summary stats)
* **Input:** Sorted Trimmed BAM (`X.trimmed.sorted.bam`)
* **Output:** Depth (`X.trimmed.sorted.depth.txt`)

## Individual Command
```bash
samtools depth -d 0 -Q 0 -q 0 -aa TRIMMED_SORTED.BAM > DEPTH.TXT
```
* Potentially gzip the output TXT files to save space?

## Batch Command
```bash
TODO
```

# Original Snakefile
```python
import os
from datetime import datetime
import time
import re
from shutil import copyfile

import pandas as pd
import re

configfile: "config.json"

illumina_reference = config["illumina"]["reference"]
ref_fasta = config["ref_fasta"]
gff = config["gff3_file"]
ont_reference = config["ont"]["reference"]
bed_file = config["bed_file"]
out_dir = config["out_dir"]
samples_path = config["samples_path"]

df = pd.read_csv(samples_path, sep="\t")
df_grp = df.groupby(["sample", "sequencing_tech"]).apply(lambda x: x.sum()) #.sum() fails for one row
lib_delim = config["library_delimiter"]

df_grp["sample_library"] = df_grp["sample_library"].apply(lambda x: x.split(lib_delim)[0] +"_" + "_".join(re.findall("L[0-9]{3}", x)))

# Illumina functions
def get_bams(wildcards):
    sample_name = df_grp[(df_grp.index.get_level_values("sequencing_tech") == "illumina") & (df_grp["sample_library"] == wildcards.sample)].index.get_level_values("sample")[0]
    libraries = df[(df["sequencing_tech"] == "illumina") & (df["sample"] == sample_name)]["sample_library"]
    return libraries.apply(lambda x: os.path.join(out_dir, "aligned_bams", "illumina", x +".sorted.bam")).tolist()

def get_forward(wildcards):
    sample_name = df_grp[(df_grp.index.get_level_values("sequencing_tech") == "illumina") & (df_grp["sample_library"] == wildcards.sample)].index.get_level_values("sample")[0]
    forward_reads = df[(df["sequencing_tech"] == "illumina") & (df["sample"] == sample_name)]["forward"]
    return forward_reads

def get_reverse(wildcards):
    sample_name = df_grp[(df_grp.index.get_level_values("sequencing_tech") == "illumina") & (df_grp["sample_library"] == wildcards.sample)].index.get_level_values("sample")[0]
    reverse_reads = df[(df["sequencing_tech"] == "illumina") & (df["sample"] == sample_name)]["reverse"]
    return reverse_reads

# Nanopore functions
def get_bams_ont(wildcards):
    sample_name = df_grp[(df_grp.index.get_level_values("sequencing_tech") == "ont") & (df_grp["sample_library"] == wildcards.sample)].index.get_level_values("sample")[0]
    libraries = df[(df["sequencing_tech"] == "ont") & (df["sample"] == sample_name)]["sample_library"]
    return libraries.apply(lambda x: os.path.join(out_dir, "aligned_bams", "ont", x +".sorted.bam")).tolist()

def get_forward_ont(wildcards):
    sample_name = df_grp[(df_grp.index.get_level_values("sequencing_tech") == "ont") & (df_grp["sample_library"] == wildcards.sample)].index.get_level_values("sample")[0]
    forward_reads = df[(df["sequencing_tech"] == "ont") & (df["sample"] == sample_name)]["forward"]
    return forward_reads

current_date_str = datetime.now().strftime("%Y-%m-%d")

rule all:
    input:
        "{out_dir}/msa/{current_date}_msa.fa".format(current_date = current_date_str, out_dir = out_dir),
        expand("{out_dir}/variants/{seq_tech}/{sample}.tsv", out_dir = out_dir, sample = df_grp["sample_library"], seq_tech = df_grp.index.get_level_values(1).unique()),
        expand("{out_dir}/depth/{seq_tech}/{sample}.depth", out_dir = out_dir, sample = df_grp["sample_library"], seq_tech = df_grp.index.get_level_values(1).unique()),
        expand("{out_dir}/barcode_counts/{seq_tech}/{sample}.tsv", out_dir = out_dir, sample = df_grp["sample_library"], seq_tech = df_grp.index.get_level_values(1).unique()),
        expand("{out_dir}/consensus_sequences/{seq_tech}/{sample}.fa", out_dir = out_dir, sample = df_grp["sample_library"], seq_tech = df_grp.index.get_level_values(1).unique()),
        "{out_dir}/trimmed_bams/illumina/reports/coverage_report.tsv".format(out_dir = out_dir),
        "{out_dir}/merged_aligned_bams/illumina/reports/mapped_unmapped_report.tsv".format(out_dir = out_dir),
        "{out_dir}/barcode_counts/illumina/contamination_report.html".format(out_dir = out_dir)


rule align_consensus_genomes:
    input:
        "{out_dir}/trimmed_bams/illumina/reports/coverage_report.tsv"
    output:
        temp("{out_dir}/msa/{current_date}.fa"),
        "{out_dir}/msa/{current_date}_msa.fa"
    threads: 4
    shell:
        """
        cat $(awk '$2 >= 70 && $1!="SAMPLE"' {input} | cut -f 1 | cut -f 1 -d . | sed 's|^|{out_dir}/consensus_sequences/illumina/|g' | sed 's/$/.fa/g') > {output[0]}
        mafft --auto --thread {threads} {output[0]} > {output[1]}
        """

rule call_depth_illumina:
    input:
        "{out_dir}/trimmed_bams/illumina/{sample}.trimmed.sorted.bam"
    output:
        "{out_dir}/depth/illumina/{sample}.depth"
    shell:
        """
        samtools depth -d 0 -Q 0 -q 0 -aa {input} > {output}
        """

rule call_variants_illumina:
    input:
        "{out_dir}/trimmed_bams/illumina/{sample}.trimmed.sorted.bam"
    output:
        "{out_dir}/variants/illumina/{sample}.tsv"
    params:
        ref= "{ref}".format(ref = ref_fasta),
        gff= "{gff}".format(gff = gff)
    shell:
        """
        samtools mpileup -A -aa -d 0 -Q 0 --reference {params.ref} {input} | ivar variants -r {params.ref}  -g {params.gff} -p {output} -m 10
        """

rule call_consensus_illumina:
    input:
        "{out_dir}/trimmed_bams/illumina/{sample}.trimmed.sorted.bam"
    output:
        "{out_dir}/consensus_sequences/illumina/{sample}.fa"
    log: "{out_dir}/logs/consensus_sequences/{sample}.log"
    shell:
        """
        samtools mpileup -aa -A -Q 0 -d 0 {input} | ivar consensus -p {output} -m 10 -n N -t 0.5 > {log}
        """

rule compute_coverage_stats:
    input:
        expand("{out_dir}/trimmed_bams/illumina/{sample}.stats", out_dir = out_dir, sample = df_grp["sample_library"])
    output:
        report="{out_dir}/trimmed_bams/illumina/reports/coverage_report.tsv",
        fig="{out_dir}/trimmed_bams/illumina/reports/coverage_report.png"
    shell:
        """
        echo -e "SAMPLE\tCOVERAGE\tAVG_DEPTH\tMIN\tMAX\tZERO_DEPTH" > {output.report}
        cat {input} | sort -n -k 2 >> {output.report}
        gnuplot -e "filename='{output.report}';ofilename='{output.fig}'" ~/code/hCoV19/plot_coverage.gnu
        """

rule trim_reads_illumina:
    input:
        "{out_dir}/merged_aligned_bams/illumina/{sample}.sorted.bam"
    output:
        bam="{out_dir}/trimmed_bams/illumina/{sample}.trimmed.sorted.bam",
        cov_stats="{out_dir}/trimmed_bams/illumina/{sample}.stats",
        bam_stats="{out_dir}/trimmed_bams/illumina/reports/{sample}_bam.stats",
        tmpbam=temp("{out_dir}/trimmed_bams/illumina/{sample}.trimmed.bam")
    params:
        bed="{bed}".format(bed = bed_file),
	plots="{out_dir}/trimmed_bams/illumina/plots/{sample}"
    log: "{out_dir}/logs/trimmed/{sample}.log"
    shell:
        """
        ivar trim -e -i {input} -b {params.bed} -p {output.tmpbam} > {log}  # Add -e if nextera used
        samtools sort -T {wildcards.sample}.trim -o {output.bam} {output.tmpbam}
        samtools index {output.bam}
	/home/gk/code/hCoV19/compute_coverage.sh {output.bam} > {output.cov_stats}
	samtools stats {output.bam} > {output.bam_stats}
        """

rule mapped_unmapped_report:
    input:
        expand("{out_dir}/merged_aligned_bams/illumina/{sample}.stats", out_dir = out_dir, sample = df_grp["sample_library"])
    output:
        "{out_dir}/merged_aligned_bams/illumina/reports/mapped_unmapped_report.tsv"
    shell:
        """
        echo -e "SAMPLE\tmapped\tunmapped" > {output}
        ls {input} | xargs -n 1 bash -c 'echo -e $0"\t"$(grep "mapped" $0 | head -n 1 | cut -f 3)"\t"$(grep "unmapped" $0 | head -n 1 | cut -f 3)' >> {output}
        """

rule demultiplex_barcodes:
    input:
        forward_barcode="{out_dir}/merged_fastq/illumina/{sample}_barcode_R1.fastq.gz",
        reverse_barcode="{out_dir}/merged_fastq/illumina/{sample}_barcode_R2.fastq.gz"
    output:
        "{out_dir}/barcode_counts/illumina/{sample}.tsv"
    params:
        forward_barcode_ref="/home/gk/code/hCoV19/db/barcodes/forward_barcodes_combined.fa",
        reverse_barcode_ref="/home/gk/code/hCoV19/db/barcodes/reverse_barcodes_combined.fa",
        barcode_dir="{out_dir}/demultiplexed_barcodes/illumina/{sample}"
    log: "{out_dir}/logs/demultiplex_barcodes/{sample}.log"
    shell:
        """
        mkdir -p {params.barcode_dir}
        cutadapt -O 15 -e 0.1 -g file:{params.forward_barcode_ref} -G file:{params.reverse_barcode_ref} -o {params.barcode_dir}/{{name1}}-{{name2}}.1.fastq.gz -p {params.barcode_dir}/{{name1}}-{{name2}}.2.fastq.gz {input.forward_barcode} {input.reverse_barcode} > {log}
        echo -e "forward_barcode\treverse_barcode\tpaired_read_count" > {output}
        find {params.barcode_dir} -name "*.fastq.gz" | sort | xargs -n 2 bash -c 'n1=$(basename $0 | sed "s/.1.fastq.gz//g" | sed "s/-/\\t/g");echo -e "${{n1}}\t"$(($(zcat $0 | wc -l)/4))' | sort -k 1 >> {output}
        """

rule extract_barcode_reads:
    input:
        forward="{out_dir}/merged_fastq/illumina/{sample}_R1.fastq.gz",
        reverse="{out_dir}/merged_fastq/illumina/{sample}_R2.fastq.gz"
    output:
        forward_barcode="{out_dir}/merged_fastq/illumina/{sample}_barcode_R1.fastq.gz",
        reverse_barcode="{out_dir}/merged_fastq/illumina/{sample}_barcode_R2.fastq.gz",
        tmpbam=temp("{out_dir}/demultiplexed_barcodes/illumina/{sample}_barcode.bam"),
        tmpbam_coord=temp("{out_dir}/demultiplexed_barcodes/illumina/{sample}_barcode.coord.bam"),
        tmpbam_name=temp("{out_dir}/demultiplexed_barcodes/illumina/{sample}_barcode.name.bam"),
        depth="{out_dir}/merged_fastq/illumina/depth/{sample}_barcode.depth"
    params:
        insert_reference="/home/gk/code/hCoV19/db/zeamays"
    shell:
        """
        bwa mem {params.insert_reference} {input.forward} {input.reverse} | samtools view -F 4 -b > {output.tmpbam}
        samtools sort -o {output.tmpbam_coord} {output.tmpbam}
        samtools index {output.tmpbam_coord}
        samtools depth -aa -d 0 {output.tmpbam_coord} > {output.depth}
        samtools sort -n -o {output.tmpbam_name} {output.tmpbam}
        samtools fastq -1 {output.forward_barcode} -2 {output.reverse_barcode} -0 /dev/null -s /dev/null {output.tmpbam_name}
        """

rule merge_multiple_libraries_illumina:
    input:
        bams = get_bams,
        forward = get_forward,
        reverse = get_reverse
    output:
        bam="{out_dir}/merged_aligned_bams/illumina/{sample}.sorted.bam",
        stats=temp("{out_dir}/merged_aligned_bams/illumina/{sample}.stats"),
        tmpbam=temp("{out_dir}/merged_aligned_bams/illumina/{sample}.with_host.sorted.bam"),
        tmp=temp("{out_dir}/merged_aligned_bams/illumina/{sample}.bam"),
        forward_merged_fastq=temp("{out_dir}/merged_fastq/illumina/{sample}_R1.fastq.gz"),
        reverse_merged_fastq=temp("{out_dir}/merged_fastq/illumina/{sample}_R2.fastq.gz")
    shell:
        """
        samtools merge {output.tmp} {input.bams}
        samtools sort -T {wildcards.sample}.merge -o {output.tmpbam} {output.tmp}
	samtools stats {output.tmpbam} > {output.stats}
        samtools view -b -F 4 {output.tmpbam} > {output.bam}
        cat {input.forward} > {output.forward_merged_fastq}
        cat {input.reverse} > {output.reverse_merged_fastq}
        """

rule align_reads_illumina:
    input:
        lambda wildcards: df[(df["sequencing_tech"] == "illumina") & (df["sample_library"]==wildcards.sample)][["forward", "reverse"]].values[0].tolist()
    output:
        temp("{out_dir}/aligned_bams/illumina/{sample}.sorted.bam")
    params:
        ref= "{ref}".format(ref = illumina_reference),
        tmp="{out_dir}/aligned_bams/illumina/{sample}.sorted.tmp.bam"
    threads: 1
    shell:
        """
        bwa mem -t {threads} {params.ref} {input[0]} {input[1]} | samtools view -Sb | samtools sort -T {wildcards.sample}.align -o {params.tmp}
        samtools addreplacerg -r "ID:{wildcards.sample}" -o {output} {params.tmp}
        rm {params.tmp}
        """

rule call_consensus_ont:
    input:
        "{out_dir}/trimmed_bams/ont/{sample}.trimmed.sorted.bam"
    output:
        "{out_dir}/consensus_sequences/ont/{sample}.fa"
    shell:
        """
        samtools mpileup -A -Q 0 -d 0 {input} | ivar consensus -p {output} -m 10 -n N -q 5
        """

rule trim_reads_ont:
    input:
        "{out_dir}/merged_aligned_bams/ont/{sample}.sorted.bam"
    output:
        "{out_dir}/trimmed_bams/ont/{sample}.trimmed.sorted.bam"
    params:
        bed="{bed}".format(bed = bed_file),
        tmp="{out_dir}/trimmed_bams/{sample}.trimmed.bam"
    shell:
        """
        ivar trim -i {input} -b {params.bed} -p {params.tmp} -q 5
        samtools sort -T {wildcards.sample}.trim -o {output} {params.tmp}
        rm {params.tmp}
        """

rule merge_multiple_libraries_ont:
    input:
        bams = get_bams_ont,
        fastq = get_forward_ont
    output:
        bam="{out_dir}/merged_aligned_bams/ont/{sample}.sorted.bam",
        fastq="{out_dir}/merged_fastq/ont/{sample}.fastq.gz"
    params:
        tmp="{out_dir}/merged_aligned_bams/ont/{sample}.bam"
    shell:
        """
        samtools merge {params.tmp} {input.bams}
        samtools sort -T {wildcards.sample}.merge -o {output.bam} {params.tmp}
        rm {params.tmp}
        cp {input.fastq} {output.fastq}
        """

rule adapter_trimming_ont:
    input:
        lambda wildcards: df[(df["sequencing_tech"] == "ont") & (df["sample_library"]==wildcards.sample)][["forward"]].values[0]
    output:
        "{out_dir}/adapter_trimmed/ont/{sample}.fastq.gz"
    shell:
        """
        porechop -i {input} -o {output}
        """

rule align_reads_ont:
    input:
        temp("{out_dir}/adapter_trimmed/ont/{sample}.fastq.gz")
    output:
        temp("{out_dir}/aligned_bams/ont/{sample}.sorted.bam")
    params:
        ref= "{ref}".format(ref = ont_reference),
        tmp="{out_dir}/aligned_bams/ont/{sample}.sorted.tmp.bam"
    threads: 1
    shell:
        """
        minimap2 -t {threads} -a -x map-ont {params.ref} {input} | samtools view -bS | samtools sort -T {wildcards.sample}.align -o {params.tmp}
        samtools addreplacerg -r "ID:{wildcards.sample}" -o {output} {params.tmp}
        rm {params.tmp}
        """

# create tmp file containing filenames for samples 
# use process substitution
rule analyse_contamination:
    input: 
        expand("{out_dir}/barcode_counts/illumina/{sample}.tsv", out_dir = out_dir, sample = df_grp["sample_library"])
    output: 
        "{out_dir}/barcode_counts/illumina/contamination_report.html"
    script: 
        "scripts/analyse_contamination.py"
```
