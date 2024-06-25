# Bayesian Bootstrap workflow


## Introduction
Felsenstein's bootstrap is the most commonly used method to measure branch support in phylogenetics. Current sequencing technologies can result in massive sampling of taxa (e.g. SARS-CoV-2). In this case, the sequences are very close, the trees are short, and the branches correspond to a small number of mutations (possibly 0). Nevertheless, these trees contain a strong signal, with unresolved parts but a low rate of false branches. With such data, Felsenstein's bootstrap is not satisfactory. Due to the frequentist nature of bootstrap sampling, the expected support of a branch corresponding to a single mutation is ~63%, even though it is highly likely to be correct.

Here we propose a Bayesian version of the phylogenetic bootstrap in which sites are assigned uninformative prior probabilities. The branch support can then be interpreted as a posterior probability. We do not view the alignment as a small subsample of a large sample of sites, but rather as containing all available information (e.g., as with complete viral genomes, which are becoming routine). We give formulas for expected supports under the assumption of perfect phylogeny, in both the frequentist and Bayesian frameworks, where a branch corresponding to a single mutation now has an expected support of ~90%. Simulation results show that these theoretical results are robust for realistic data. Results with viral and non-viral low-homoplasy datasets show that Bayesian bootstrap supports are easier to interpret, with high supports for branches very likely to be correct. As homoplasy increases, the two supports become closer and strongly correlated.

This repository contains BBOOT, the workflow that computes baysien supports from a multiple sequence alignment (MSA). To do so, it first generates weight vectors (following a Dirichlet distribution), then infers reference and bootstrap trees using PhyML, and finally computes Bayesian supports. In addition, it produces a file containing several metrics and outputs:

1. MSA alphabet (nucleotides or proteins);
2. MSA length;
3. Minimal number of mutations in the MSA (see text);
4. Reference tree size (sum of branch lengths, in number of mutations);
5. Level of homoplasy.

Finally, it produces a figure showing the cumulative distribution of branch lengths.


## Prerequisites

Java and Apptainer/Singularity must be installed.

## Configuration file

If you want to execute BBOOT locally on your linux / macos machine, it should work as is.
However, if you want to execute BBOOT on a specific environment, you may modify the configuration file, especially the "profiles" section. Several profiles are already defined:

1. standard: Standard configuration of processes (containers, etc.)
2. local: To run the workflow locally (no slurm / HPC)
3. slurm: To run the workflow on a SLURM HPC (you need to modify the "queue" and "clusterOptions" values
4. singularity: To run the workflow using singularity images.
5. arm64: If you are on an arm64 architecture, you may use this option

## Executing the workflow

Execute the workflow with the following command:

```
nextflow run main.nf --msa <MSA>
      --results <OUT DIR: results> 
      --nboot <# BOOT REPLICATES: 200> 
      --collapse <COLLAPSE THRESHOLD: 0.1>
      -profile singularity,local
```

## Example dataset

You can try the workflow with the example dataset:

```
nextflow run main.nf --msa sample_data/sample.fasta.gz
      --results results
      --nboot 1000
      --collapse 0.1
      -profile standard,singularity,local
```
