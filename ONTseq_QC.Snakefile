################################################################################
# Snakemake pipeline for base-calling, de-multiplexing and quality control of
# raw MinION Nanopore data
#
# Steps:
#   1. Basecalling with Guppy either using the fast or high-accuracy mode
#   2. De-multiplexing with Guppy
#   3. Create a single FastQ file for expected barcodes (specified via CSV file)
#   4. Summary of the quality of the sequencing run (summary using MultiQC)
#     a. FastQC
#     b. pycoQC
#     c. NanoStat
#   5. Taxonomic profile
#
# Alex Huebner, 10/11/2020
################################################################################

from glob import glob
import os

import pyfastx

workdir: config['tmpdir']

#### Set-up workspace ##########################################################
DATADIR = config['datadir']
PROJNAME = config['projname']
PROJDIR = config['projdir']
BCTYPE = config['basecalltype']
GUPPYDIR = config['guppydir']
################################################################################

#### Fast5 #####################################################################
FAST5FILES = glob(f"{DATADIR}/**/*.fast5", recursive=True)
FAST5DIRS = list(set([os.path.dirname(fn) for fn in FAST5FILES]))
BARCODES = [line.rstrip() for line in open(config['expectedbarcodes'], "rt")]
################################################################################

os.makedirs(PROJDIR, exist_ok=True)
os.makedirs("cluster_logs", exist_ok=True)

rule all:
    input: 
        expand("{projdir}/fastq/{barcode}.fastq.gz", projdir=[PROJDIR], barcode=BARCODES),
        f"{PROJDIR}/logs/pycoqc/pycoqc_summary.html",
        f"{PROJDIR}/logs/multiqc_report.html",
        f"{PROJDIR}/logs/nreads_per_barcode.txt"


rule guppy_basecalling:
    output:
        touch("{projname}/guppy_basecalling_{bctype}/{i}_basecalling.done")
    message: "Perform basecalling with guppy on folder {wildcards.i} of project {wildcards.projname} with basecaller type {wildcards.bctype}"
    params:
        fast5dir = lambda wildcards: FAST5DIRS[int(wildcards.i)],
        fastqdir = "{projname}/guppy_basecalling_{bctype}/{i}",
        config = lambda wildcards: f"{GUPPYDIR}/data/dna_r9.4.1_450bps_{wildcards.bctype}.cfg" 
    threads: 16
    shell:
        """
        PATH={GUPPYDIR}/bin:$PATH \
        guppy_basecaller -i {params.fast5dir} \
                         -r \
                         -s {params.fastqdir} \
                         -c {params.config} \
                         --compress_fastq \
                         --num_callers 4
        """

rule guppy_demultiplexing:
    input:
        [f"{PROJNAME}/guppy_basecalling_{BCTYPE}/{i:04d}_basecalling.done" for i in range(len(FAST5DIRS))]
    output:
        touch("{projname}/guppy_demultiplexer.done")
    message: "Demultiplex output with Guppy"
    params:
        dir = f"{PROJNAME}/guppy_basecalling_{BCTYPE}",
        outdir = f"{PROJNAME}/guppy_demultiplexing"
    threads: 16
    shell:
        """
        PATH={GUPPYDIR}/bin:$PATH \
        guppy_barcoder -i {params.dir} \
                       -r \
                       -s {params.outdir} \
                       --compress_fastq \
                       --trim_barcodes \
                       -t {threads}
        """

rule noreads_per_barcode:
    input:
        f"{PROJNAME}/guppy_demultiplexer.done"
    output:
        "{projdir}/logs/nreads_per_barcode.txt"
    message: "Determine the number of reads per barcode"
    params:
        dir = f"{PROJNAME}/guppy_demultiplexing"
    run:
        with open(output[0], "wt") as outfile:
            print("barcode\tnReads", file = outfile)
            for dn in glob(f'{params.dir}/barcode*') + [f'{params.dir}/unclassified']:
                i = 0
                for fn in glob(f'{dn}/*.fastq.gz'):
                    for _, _, _ in pyfastx.Fastq(fn, build_index=False):
                        i += 1
                print(f"{os.path.basename(dn)}\t{i}", file = outfile)

rule concat_expected_barcodes:
    input:
        f"{PROJNAME}/guppy_demultiplexer.done"
    output:
        "{projdir}/fastq/{barcode}.fastq.gz"
    message: "Concatenate FastQ file of barcode {wildcards.barcode}"
    params:
        dir = f"{PROJNAME}/guppy_demultiplexing"
    shell:
        """
        cat {params.dir}/{wildcards.barcode}/*.fastq.gz > {output}
        """

rule nanostat:
    input:
        "{projdir}/fastq/{barcode}.fastq.gz"
    output:
        "{projdir}/logs/nanostat/{barcode}.nanostat.txt"
    message: "Infer quality statistics for barcode {wildcards.barcode}"
    conda: f"{workflow.basedir}/ONTseq_QC.yaml"
    params:
        outdir = "{projdir}/logs/nanostat",
        outname = "{barcode}.nanostat.txt"
    threads: 4
    shell:
        """
        NanoStat --fastq {input} \
                 --outdir {params.outdir} \
                 --outname {params.outname} \
                 --threads {threads}
        """

rule fastqc:
    input:
        "{projdir}/fastq/{barcode}.fastq.gz"
    output:
        "{projdir}/logs/fastqc/{barcode}/{barcode}_fastqc.html"
    message: "Run FastQC on barcode {wildcards.barcode}"
    conda: f"{workflow.basedir}/ONTseq_QC.yaml"
    params:
        outdir = "{projdir}/logs/fastqc/{barcode}"
    threads: 4
    shell:
        """
        fastqc -o {params.outdir} -t {threads} {input}
        """

rule concat_summaryreports:
    input:
        [f"{PROJNAME}/guppy_basecalling_{BCTYPE}/{i:04d}_basecalling.done" for i in range(len(FAST5DIRS))]
    output:
        f"{PROJNAME}/logs/pycoqc/sequencing_summary.txt"
    message: "Concatenate the sequencing summary files of Guppy basecalling for pycoQC analysis"
    run:
        for i, dummyfn in enumerate(input):
            summaryfn = dummyfn.replace("_basecalling.done", "/sequencing_summary.txt")
            if i == 0:
                with open(output[0], "wb") as outfile:
                    for line in open(summaryfn, "rb"):
                        outfile.write(line)
            else:
                with open(output[0], "ab") as outfile:
                    for j, line in enumerate(open(summaryfn, "rb")):
                        if j > 0:
                            outfile.write(line)

rule pycoqc:
    input:
        f"{PROJNAME}/logs/pycoqc/sequencing_summary.txt"
    output:
        html = f"{PROJDIR}/logs/pycoqc/pycoqc_summary.html",
        json = f"{PROJDIR}/logs/pycoqc/pycoqc_summary.json"
    message: "Run PycoQC to analyse the quality of the run"
    conda: f"{workflow.basedir}/ONTseq_QC.yaml"
    params:
        barcode = f"{PROJNAME}/guppy_demultiplexing/barcoding_summary.txt"
    shell:
        """
        pycoQC -f {input} \
               -b {params.barcode} \
               -o {output.html} \
               -j {output.json} \
               --min_pass_len 250
        """

rule multiqc:
    input:
        expand("{projdir}/logs/fastqc/{barcode}/{barcode}_fastqc.html", projdir=[PROJDIR], barcode=BARCODES)
    output:
        "{projdir}/logs/multiqc_report.html"
    message: "Summarise the FastQC reports using MultiQC"
    conda: f"{workflow.basedir}/ONTseq_QC.yaml"
    params:
        dir = "{projdir}/logs/fastqc"
    shell:
        """
        cd $(dirname {params.dir})
        multiqc --force {params.dir}
        """

