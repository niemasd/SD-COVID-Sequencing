# Phylogenetic Inference
Given a bunch of assembled SARS-CoV-2 genomes, this document will provide step-by-step commands to align the genomes and infer a phylogeny.

# Step 1: Multiple Sequence Alignment (MSA)
* **Input:** FASTA file containing unaligned genomes (`X.fas`)
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
trim_msa.py -i viralmsa_out/consensus.fas.aln -o consensus.trimmed.aln -s 100 -e 50
```
