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
