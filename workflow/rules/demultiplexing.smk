def return_renamed_fastqs(wildcards):
    barcodes = [line.strip() for line in open(checkpoints.summarise_expected_barcodes.get(**wildcards).output[0], "rt")]
    return [f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/fastq/{barcode}/{barcode}.fastq.gz" for barcode in barcodes]

localrules: concat_fastq

rule demultiplexing:
    input:
        return_renamed_fastqs
    output:
        touch("{projdir}/guppy_basecalling_{accuracymode}/demultiplexing.done")

rule guppy_demultiplexing:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/guppy/sequencing_summary.txt"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/barcoding_summary.txt"
    message: "Demultiplex output with Guppy"
    params:
        guppydir = config['guppydir'],
        dir = "{projdir}/guppy_basecalling_{accuracymode}/guppy/pass",
        outdir = "{projdir}/guppy_basecalling_{accuracymode}/fastq",
        barcodekit = config['barcodekit']
    threads: 16
    shell:
        """
        PATH={params.guppydir}/bin:$PATH \
        guppy_barcoder -i {params.dir} \
                       -r \
                       -s {params.outdir} \
                       --compress_fastq \
                       --enable_trim_barcodes \
                       --trim_adapters \
                       -t {threads} \
                       --barcode_kits {params.barcodekit}
        """

checkpoint summarise_expected_barcodes:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/barcoding_summary.txt"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/expected_barcodes.txt"
    message: "Write expected barcodes to file"
    run:
        with open(output[0], "wt") as outfile:
            barcodes = pd.read_csv(input[0], sep="\t")['barcode_arrangement'] \
                .unique().tolist()
            for barcode in sorted(barcodes):
                outfile.write(barcode + "\n")


rule concat_fastq:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/expected_barcodes.txt"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/{barcode}/{barcode}.fastq.gz"
    message: "Rename the FastQ to unified name: {wildcards.barcode}"
    params:
        fastq = lambda wildcards: " ".join(glob(f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/fastq/{wildcards.barcode}/*.fastq.gz")),
        n_fastq = lambda wildcards: len(glob(f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/fastq/{wildcards.barcode}/*.fastq.gz")) 
    shell:
        """
        if [[ {params.n_fastq} -eq 1 ]]; then
            mv {params.fastq} {output}
        else
            cat {params.fastq} > {output} && rm {params.fastq}
        fi
        """
