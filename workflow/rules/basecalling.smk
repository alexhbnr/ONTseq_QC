rule guppy_basecalling:
    input:
        f"{config['tmpdir']}/prepare_data.done"
    output:
        "{projdir}/guppy_basecalling_{accuracymode}/guppy/sequencing_summary.txt"
    message: "Perform basecalling with guppy"
    resources:
        mem = 24,
        cores = 24,
        cuda = 1
    params:
        guppydir = config['guppydir'],
        rawseqdir = f"{config['tmpdir']}/seqfiles_dirs",
        fastqdir = "{projdir}/guppy_basecalling_{accuracymode}/guppy",
        config = f"{config['guppydir']}/data/dna_r{config['flowcell_generation']}_{config['motor_enzyme']}_{config['accuracymode']}.cfg",
        gpu_mode = "-x cuda:all" if config['guppymode'] == "gpu" else ""
    threads: 24
    shell:
        """
        PATH={params.guppydir}/bin:$PATH \
        guppy_basecaller -i {params.rawseqdir} \
                         -r {params.gpu_mode} \
                         -s {params.fastqdir} \
                         -c {params.config} \
                         --compress_fastq \
                         --num_callers 4 \
                         --cpu_threads_per_caller {threads}
        """
