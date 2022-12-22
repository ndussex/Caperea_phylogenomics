# Caperea phylogenomics

This repository contains the data related to the publication:

Ludovic Dutoit,† Kieren J. Mitchell†, Nicolas Dussex, †, Catherine M. Kemper, Petter Larsson, Love Dalén, Nicolas J. Rawlence and , Felix G. Marx. *Convergent evolution of skim feeding in baleen whales.*


It is meant to ensure reproducibility and trace analyses in the paper. Please get in touch with nicolas.dussex@gmail.com  and/or ludovic.dutoit@otago.ac.nz for any questions. 

## Genome consensus building

The [bcftools_consensus_genome.bsh](bcftools_consensus_genome.bsh) script uses a bcf file and an assembly to build a new consensus genome

## SNP calling

This is recorded [SNPCalling_clean.md](SNPCalling_clean.md) that includes masking and splitting all the way to one file per window with all species.

## Tree building

Saved as one file in [build_trees_clean.md](build_trees_clean.md) that includes the astral analysis. Astral visualisations are not included, but scripts are availabe upon request.

## PSMC

The [bcftools_consensus_genome.bsh](bcftools_consensus_genome.bsh) is used to split a bam file into diplois.fastq file which is then used as input in the [Pygmy_whale_PSMC_bootstrap.bsh](Pygmy_whale_PSMC_bootstrap.bsh) to run the PSMC.
