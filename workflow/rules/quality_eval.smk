#### Auxilliary functions ######################################################

def expected_nread_files(wildcards):
    barcodes = [line.strip() for line in open(checkpoints.summarise_expected_barcodes.get(**wildcards).output[0], "rt")]
    return [f"{config['tmpdir']}/nreads/{barcode}.n" for barcode in sorted(barcodes)]


def expected_nanostat_files(wildcards):
    barcodes = [line.strip() for line in open(checkpoints.summarise_expected_barcodes.get(**wildcards).output[0], "rt")]
    return [f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/logs/nanostat/{barcode}.nanostat.txt" for barcode in barcodes]


def expected_fastqc_files(wildcards):
    barcodes = [line.strip() for line in open(checkpoints.summarise_expected_barcodes.get(**wildcards).output[0], "rt")]
    return [f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/logs/fastqc/{barcode}/{barcode}_fastqc.html" for barcode in barcodes]

################################################################################

rule quality_evalution:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/logs/number_reads.tsv",
        "{projdir}/guppy_basecalling_{accuracymode}/logs/pycoqc/pycoqc_summary.html",
        "{projdir}/guppy_basecalling_{accuracymode}/logs/multiqc_report.html",
        expected_nanostat_files
    output:
        touch("{projdir}/guppy_basecalling_{accuracymode}/quality_evaluation.done")


#### Number of reads per barcode ###############################################

rule count_reads:
    input:
        f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/demultiplexing.done"
    output:
        "{tmpdir}/nreads/{barcode}.n"
    message: "Count the number of reads: {wildcards.barcode}"
    conda: "../envs/ENVS_bioawk.yaml"
    params:
        fastq = lambda wildcards: glob(f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/fastq/{wildcards.barcode}/*.fastq.gz")[0]
    shell:
        """
        echo -e "{wildcards.barcode}\t$(bioawk -c fastx 'END{{print NR}}' {params.fastq})" > {output}
        """

rule summarise_nreads:
    input:
        expected_nread_files
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/logs/number_reads.tsv"
    message: "Concatenate the number of reads into a single file"
    shell:
        """
        echo -e "barcode\tnReads\n" > {output}
        cat {input} >> {output}
        """

################################################################################

#### Nanostat ##################################################################

rule nanostat:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/{barcode}/{barcode}.fastq.gz"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/logs/nanostat/{barcode}.nanostat.txt"
    message: "Infer quality statistics for barcode {wildcards.barcode}"
    conda: "../envs/ENVS_nanostat.yaml"
    params:
        outdir = "{projdir}/guppy_basecalling_{accuracymode}/logs/nanostat",
        outname = "{barcode}.nanostat.txt"
    threads: 4
    shell:
        """
        NanoStat --fastq {input} \
                 --outdir {params.outdir} \
                 --name {params.outname} \
                 --threads {threads}
        """

################################################################################

#### FastQC ####################################################################

rule fastqc:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/{barcode}/{barcode}.fastq.gz"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/logs/fastqc/{barcode}/{barcode}_fastqc.html"
    message: "Run FastQC on barcode {wildcards.barcode}"
    conda: "../envs/ENVS_fastqc.yaml"
    params:
        outdir = "{projdir}/guppy_basecalling_{accuracymode}/logs/fastqc/{barcode}"
    threads: 4
    shell:
        """
        fastqc -o {params.outdir} -t {threads} {input}
        """

################################################################################

#### PycoQC ####################################################################

rule pycoqc:
    input:
        summary = "{projdir}/guppy_basecalling_{accuracymode}/guppy/sequencing_summary.txt",
        barcode = "{projdir}/guppy_basecalling_{accuracymode}/fastq/barcoding_summary.txt"
    output:
        html = "{projdir}/guppy_basecalling_{accuracymode}/logs/pycoqc/pycoqc_summary.html",
        json = "{projdir}/guppy_basecalling_{accuracymode}/logs/pycoqc/pycoqc_summary.json"
    message: "Run PycoQC to analyse the quality of the run"
    conda: "../envs/ENVS_pycoqc.yaml"
    params:
    shell:
        """
        pycoQC -f {input.summary} \
               -b {input.barcode} \
               -o {output.html} \
               -j {output.json} \
               --min_pass_len 200
        """

################################################################################

#### MultiQC ###################################################################

rule multiqc:
    input:
        expected_fastqc_files
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/logs/multiqc_report.html"
    message: "Summarise the FastQC reports using MultiQC"
    conda: "../envs/ENVS_multiqc.yaml"
    params:
        dir = "{projdir}/guppy_basecalling_{accuracymode}/logs"
    shell:
        """
        cd $(dirname {params.dir})
        multiqc --force {params.dir} -o {params.dir} --exclude pycoqc
        """

################################################################################
