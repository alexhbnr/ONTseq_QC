# ONTseq pipeline

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for processing Oxford
Nanopore Technology (ONT) sequencing data (base-calling and de-multiplexing) followed by some basic
quality control assessment.

While the pipeline should be universal for all types of ONT sequencing data, it has so far only been
tested for data generated on the MinION flowcell *FLO-MIN106* with the ligation sequencing kit
*SQK-LSK109*. However, you can specify the correct profile for base-calling with *guppy* in the
Snakemake file.

The workflow assumes that it is started and run on a Unix-based system (Linux, Mac OS X).

## Prerequisites

  1. Conda with Python3

  In order to be able to run the workflow, a working installation of *conda* and *Python3* is
  required.  This is easiest way of obtaining it is to download the **Miniconda3 Python3
  installer** from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
  for your operating system and to follow the instructions on the website.

  2. Python packages

  Next to *conda* and *Python3*, there are two more Python modules that have to be installed prior
  to usage: [snakemake](https://snakemake.readthedocs.io/) and
  [pyfastx](https://pyfastx.readthedocs.io/en/latest/).

  You can install these either using *conda*

  `conda install snakemake pyfastx`

  or via the Python module manager *pip*

  `pip install snakemake pyfastx`

  3. ONT guppy

  Finally, the workflow uses the Oxford Nanopore Technology (ONT) software *guppy*, which is
  proprietary and only available via the [ONT community
  website](https://community.nanoporetech.com/downloads). Please download the binaries suitable for
  your operating system and unpack them.  

## Workflow configuration

  This *snakemake* workflow uses a configuration file in JSON format that we will provide as input
  to *snakemake* using the `--configfile` parameter. A template with the required fields is
  available under `ONTseq_QC-config.template`.

  The configuration file expects the following values:

  - `projname`: name of the project that will be used to create a sub-folder in temporary directory
  - `datadir`: directory that is used to search for Fast5 files generate by the ONT sequencer
  - `tmpdir`: directory in which the temporary output of the analysis will be stored
  - `projdir`: directory in which the final output of the analysis will be stored
  - `expectedbarcodes`: list of the barcodes, one per line, that were used during library
     preparation and are therefore expected to be observed, e.g. *barcode01*
  - `basecalltype`: base calling mode to be used in *guppy*: hac or fast
  - `guppydir`: directory to which the files downloaded from the ONT community website were
    installed to; should contain both the *bin* as well as the *data* sub-folders

  Make a local copy of the template configuration file and fill in the information that fit your
  data.

## Running the workflow

  After setting up the required programs and generating the configuration file, we can start our
  workflow. For the execution on a local machine, we can simply run:

  ```
  snakemake -s ONTseq_QC.Snakefile \
            --configfile <config file> \
            --use-conda \
            --cores <number of cores>
  ```

  Please change the path of the config file to the config file that you just created and adjust the
  number of cores to the desired number available on your system.

  The option `--use-conda` will force *snakemake* to create a *conda* environment for this pipeline
  with the name `ONTseq_QC` and install the additional software prerequisites on the fly. These
  additional tools are:

  - [NanoStat](https://github.com/wdecoster/nanostat)
  - [FastQC](https://github.com/s-andrews/FastQC)
  - [MultiQC](https://multiqc.info/)
  - [pycoQC](https://github.com/a-slide/pycoQC)

  If you want to avoid having to create a new temporary *conda* environment over and over again, you
  can additional specify the path in which this conda environment is created using `--conda-prefix`.

  In case a computing cluster with a scheduling software such as SLURM is available, we can enable
  the submission to the cluster by using:

  ```
  snakemake -s ONTseq_QC.Snakefile \
            --configfile <config file> \
            --cluster-config ONTseq_QC-SLURM.json \
            --cluster "sbatch --mem {cluster.mem} -p {cluster.partition} -o {cluster.out} -e
            {cluster.err} -c {threads}" \
            --use-conda \
            --cores <number of cores>
  ```

  The actual configuration of scheduler might likely differ due to your local set-up and can be
  adjusted by altering the JSON configuration file `ONTseq_QC-SLURM.json`.

## Output

  During its runtime, the workflow will write all temporary output to the temporary folder inside
  the sub-folder of the project name, both specified in the configuration file. The final output
  will be written to the project folder specified in the same file.

  Here is the overview of the output that is generated in the project folder:

  1. `fastq`: contains the base-called, de-multiplexed FastQ files, one file per expected barcode
  2. `logs/nreads_per_barcode.txt`: a table with the number of demultiplexed reads for every barcode
     observed
  3. `logs/pycoqc`: the HTML and JSON report generated by PycoQC across the whole sequencing run
  4. `logs/nanostat`: the simple text summary of the sequencing results per barcode produced by
     NanoStat
  5. `logs/fastqc`: the output of FastQC for each expected barcode
  6. `logs/multiqc_report.html`: the summary HTML report of the FastQC results by MultiQC across all
     barcodes

  The temporary output folder will not be deleted automatically in order to provide the possibility
  for inspecting the output of the individual programs. If this temporary output is not required any
  longer, it can be just simply deleted.
