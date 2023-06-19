#### Auxilliary functions ######################################################

def expected_porechop_files(wildcards):
    barcodes = [line.strip() for line in open(checkpoints.summarise_expected_barcodes.get(**wildcards).output[0], "rt")]
    return [f"{config['projdir']}/guppy_basecalling_{config['accuracymode']}/logs/porechop/{barcode}.porechop.log" for barcode in barcodes]

################################################################################

#### Evaluate chimeric reads with porechop #####################################

rule porechop:
    input:
        "{projdir}/guppy_basecalling_{accuracymode}/demultiplexing/{barcode}/{barcode}.fastq.gz"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/fastq/{barcode}.fastq.gz"
    message: "Identify remaining barcodes and adapter sequences after guppy processing: {wildcards.barcode}"
    conda: "../envs/ENVS_porechop.yaml"
    log: "{projdir}/guppy_basecalling_{accuracymode}/logs/porechop/{barcode}.porechop.log"
    threads: 8
    shell:
        """
        porechop_abi --ab_initio -i {input} -o {output} -v 2 -t {threads} > {log}
        """

################################################################################
