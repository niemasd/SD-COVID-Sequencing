# Phylogenetic Inference
Given a bunch of assembled SARS-CoV-2 genomes, this document will provide step-by-step commands to align the genomes and infer a phylogeny.

# Step 1: Multiple Sequence Alignment (MSA)
* **Input:** FASTA file containing unaligned genomes (`X.fas`), which should include the [outgroup sequence (RmYN02)](reference_genome/RmYN02.fas)
* **Output:** FASTA file containing the MSA (`X.aln`)

## Command
```bash
ViralMSA.py -s UNALIGNED.FAS -r SARS-CoV-2 -e EMAIL_ADDRESS -o MSA_OUTPUT_DIR -t THREADS
```

## Example
```bash
ViralMSA.py -s consensus.fas -r SARS-CoV-2 -e niemamoshiri@gmail.com -o viralmsa_out -t 32
```

# Step 2: Remove Ends of the MSA
* **Input:** FASTA file containing the MSA (`X.aln`)
* **Output:** FASTA file containing the MSA with the ends removed (`X.trimmed.aln`)

The reason for this is that the beginning and end of the SARS-CoV-2 assemblies may be error-prone, so we trim the ends out to only include high-confidence variants. According to [NextStrain](https://github.com/nextstrain/ncov/blob/b61864fba9c5cfd5b5b9a52518f9096a9e631a6e/defaults/parameters.yaml#L75), it's reasonable to cut the first 100 and last 50 bases with respect to the reference genome. I wrote a script to help with this: [`trim_msa.py`](https://github.com/niemasd/tools/blob/master/trim_msa.py)

## Command
```bash
trim_msa.py -i ALIGNED.ALN -s TRIM_FROM_START -e TRIM_FROM_END -o ALIGNED.TRIMMED.ALN
```

## Example
```bash
trim_msa.py -i viralmsa_out/consensus.fas.aln -s 100 -e 50 -o consensus.trimmed.aln
```

# Step 3: Infer Phylogeny
* **Input:** FASTA file containing the MSA with ends removed (`X.trimmed.aln`)
* **Output:** Newick file containing the unrooted phylogeny (`X.unrooted.nwk`)

## Command
```bash
iqtree2 -T THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s ALIGNED.TRIMMED.ALN
```
* Optionally add `--fast` to use the less accurate but much faster tree search algorithm

## Example
```bash
iqtree2 -T 32 -m GTR+F+G4 --polytomy -blmin 1e-9 -s consensus.trimmed.aln
```

# Rooting
This gets us an unrooted tree, but we actually need to root the tree. I'll try to get this written up ASAP, but basically, we should be able to do outgroup rooting using [RmYN02](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211627/) as the outgroup.

1. Add the [outgroup RmYN02 sequence](reference_genome/RmYN02.fas) to the collection of sequences before Multiple Sequence Alignment
2. Do the whole phylogenetic inference workflow (Steps 1-3)
3. Use [FastRoot](https://github.com/uym2/MinVar-Rooting) to root the tree on the outgroup

# Step 4: Pangolin Lineage Assignment
* **Input:** FASTA file containing unaligned genomes (`X.fas`)
* **Output:** Assigned Pangolin lineages (`X.pangolin.csv`)

I believe we would want to run this on the original unaligned sequences, rather than on the MSA or trimmed MSA. Thus, this can happen in parallel to Steps 1-3.

## How to Install
Here's how I installed Pangolin in my environment (rather than using their conda environment):

```bash
git clone https://github.com/cov-lineages/pangolin.git && cd pangolin && python setup.py install && cd .. && rm -rf pangolin
git clone https://github.com/cov-lineages/pangoLEARN.git && cd pangoLEARN && python setup.py install && cd .. && rm -rf pangoLEARN
git clone https://github.com/cov-ert/datafunk.git && cd datafunk && python setup.py install && cd .. && rm -rf datafunk
pip install snakemake
```

## Command
```bash
pangolin --update
pangolin -t THREADS UNALIGNED.FAS
```

## Example
```bash
pangolin --update
pangolin -t 32 consensus.fas
```
