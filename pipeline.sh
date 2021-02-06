#!/usr/bin/env bash

# default constants
THREADS=1
REF_FAS='../ref/NC_045512.2.fas'
REF_MMI='../ref/NC_045512.2.fas.mmi'
REF_GFF='../ref/NC_045512.2.gff3'
PRIMER_BED='../primers/swift/sarscov2_v2_primers.bed'

# check usage
if [[ "$#" -ne 1 && "$#" -ne 2 ]] ; then
    echo "USAGE: $0 sample_prefix [num_threads]" ; exit 1
elif [[ "$#" -eq 2 ]] ; then
    THREADS=$2
fi

# Step 1: Map Reads + Sort
{ time ( minimap2 -t $THREADS -a -x sr "$REF_MMI" $1*.fastq.gz | samtools sort --threads $THREADS -o $1.sorted.bam ) ; } 2> $1.log.1.map.log && \

# Step 2: Trim Sorted BAM
{ time ( ivar trim -x 5 -e -i $1.sorted.bam -b "$PRIMER_BED" -p $1.trimmed ) ; } > $1.log.2.trim.log 2>&1 && \

# Step 3: Sort Trimmed BAM
{ time ( samtools sort --threads $THREADS -o $1.trimmed.sorted.bam $1.trimmed.bam && rm $1.trimmed.bam ) ; } 2> $1.log.3.sorttrimmed.log && \

# Step 4: Generate Pile-Up
{ time ( samtools mpileup -A -aa -d 0 -Q 0 --reference "$REF_FAS" $1.trimmed.sorted.bam ) ; } > $1.trimmed.sorted.pileup.txt 2> $1.log.4.pileup.log && \

# Step 5: Call Variants
{ time ( cat $1.trimmed.sorted.pileup.txt | ivar variants -r "$REF_FAS" -g "$REF_GFF" -p $1.trimmed.sorted.pileup.variants.tsv -m 10 ) ; } 2> $1.log.5.variants.log && \

# Step 6: Call Consensus
{ time ( cat $1.trimmed.sorted.pileup.txt | ivar consensus -p $1.trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ) ; } > $1.log.6.consensus.log 2>&1 && \

# Step 7: Call Depth
{ time ( samtools depth -d 0 -Q 0 -q 0 -aa $1.trimmed.sorted.bam ) ; } > $1.trimmed.sorted.depth.txt 2> $1.log.7.depth.log && \

# Step 8: Qualimap
for x in sorted trimmed.sorted ; do
    { time ( qualimap bamqc -bam $1.$x.bam -nt $THREADS --java-mem-size=4G -outdir $1.$x.stats && tar c $1.$x.stats | pigz -9 -p $THREADS > $1.$x.stats.tar.gz && rm -rf $1.$x.stats ) ; } > $1.log.8.qualimap.$x.log 2>&1
done && \

# Step 9: Zip
zip -9 $1.zip $1*
