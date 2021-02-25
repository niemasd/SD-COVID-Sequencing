# Phylogenetic Inference
Given a bunch of assembled SARS-CoV-2 genomes, this document will provide step-by-step commands to align the genomes, infer an unrooted ML phylogeny, root the phylogeny, and assign the genomes to pangolin lineages.

# Step 1: Multiple Sequence Alignment (MSA)
* **Input:** FASTA file containing unaligned genomes (`X.fas`), which should include the [outgroup sequence (RmYN02)](reference_genome/RmYN02.fas)
* **Output:** FASTA file containing the MSA (`X.aln`)
* **Tool(s):** [ViralMSA](https://github.com/niemasd/ViralMSA) and [Minimap2](https://github.com/lh3/minimap2)

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
* **Tool(s):** [`trim_msa.py`](https://github.com/niemasd/tools/blob/master/trim_msa.py)

The reason for this is that the beginning and end of the SARS-CoV-2 assemblies may be error-prone, so we trim the ends out to only include high-confidence variants. According to [NextStrain](https://github.com/nextstrain/ncov/blob/b61864fba9c5cfd5b5b9a52518f9096a9e631a6e/defaults/parameters.yaml#L75), it's reasonable to cut the first 100 and last 50 bases with respect to the reference genome. I wrote a script to help with this: [`trim_msa.py`](https://github.com/niemasd/tools/blob/master/trim_msa.py)

## Command
```bash
trim_msa.py -i ALIGNED.ALN -s TRIM_FROM_START -e TRIM_FROM_END -o ALIGNED.TRIMMED.ALN
```

## Example
```bash
trim_msa.py -i viralmsa_out/consensus.fas.aln -s 100 -e 50 -o consensus.trimmed.aln
```

# Step 3.1: Infer Unrooted Phylogeny
* **Input:** FASTA file containing the MSA with ends removed (`X.trimmed.aln`)
* **Output:** Newick file containing the unrooted phylogeny (`X.unrooted.nwk`)
* **Tool(s):** [IQ-TREE 2](http://www.iqtree.org/)

## Command
```bash
iqtree2 -T THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s ALIGNED.TRIMMED.ALN
```
* Optionally add `--fast` to use the less accurate but much faster tree search algorithm

## Example
```bash
iqtree2 -T 32 -m GTR+F+G4 --polytomy -blmin 1e-9 -s consensus.trimmed.aln
```

# Step 3.2: Root the Tree
* **Input:** Newick file containing the unrooted phylogeny (`X.unrooted.nwk`)
* **Output:** Newick file containing the rooted phylogeny (`X.rooted.nwk`)
* **Tool(s):** [FastRoot](https://github.com/uym2/MinVar-Rooting)

## Command
```bash
FastRoot.py -i UNROOTED_TREE.NWK -m OG -g OUTGROUP
```

## Example
```bash
FastRoot.py -i consensus.trimmed.aln.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"
```

# Step 3.3: Rename the Leaves
The leaves will have ugly names that include a bunch of info beyond the sample name, but we may want to strip out the extra info to just have labels be sample names (Yoshiki from Rob's lab asked for this). I wrote the following script (`rename_samples.py`) to do so:

```python
#!/usr/bin/env python3
from sys import argv
if len(argv) != 3:
    print("USAGE: %s <input_tree> <output_tree>" % argv[0]); exit(1)
import re
pattern = re.compile('_(L00\d|S\d{1,4})[_\.].*')
from treeswift import read_tree_newick
tree = read_tree_newick(argv[1])
for node in tree.traverse_leaves():
    match = pattern.search(node.label)
    if match:
        node.label = node.label[:match.start(0)]
f = open(argv[2],'w'); f.write(tree.newick()); exit()
```

## Command
```bash
./rename_samples ROOTED_TREE.NWK RENAMED_ROOTED_TREE.NWK
```

# Step 4: Pangolin Lineage Assignment
* **Input:** FASTA file containing unaligned genomes (`X.fas`)
* **Output:** Assigned Pangolin lineages (`X.pangolin.csv`)
* **Tool(s):** [pangolin](https://github.com/cov-lineages/pangolin)

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
