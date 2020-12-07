# ONTseq pipeline

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for processing of Oxford
Nanopore Technology (ONT) sequencing data (base-calling and de-multiplexing) followed by some basic
quality control assessment.

While the pipeline should be universal for all types of ONT sequencing data, it has so far only been
tested for data generated on the MinION flowcell *FLO-MIN106* with the ligation sequencing kit
*SQK-LSK109*. However, you can specify the correct profile for base-calling with *guppy* in the
config file.

## Prerequisites

  1. Conda with Python3

  2. Python packages

  3. ONT guppy


