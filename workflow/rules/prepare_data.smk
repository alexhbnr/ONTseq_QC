localrules: link_seqfile

if config['fileformat'] == "fast5" or config['fileformat'] == "pod5":

    seqfilenames = glob(f"{config['datadir']}/**/*.{config['fileformat']}", recursive=True)
    seqids = {os.path.basename(fn).replace(f".{config['fileformat']}", ""): fn for fn in seqfilenames}

    rule prepare_data:
        input:
            expand("{tmpdir}/seqfiles_dirs/{seqid}.{ff}",
                   tmpdir=[config['tmpdir']],
                   seqid=seqids,
                   ff=[config['fileformat']])
        output:
            touch(f"{config['tmpdir']}/prepare_data.done")

    rule link_seqfile:
        output:
            "{tmpdir}/seqfiles_dirs/{seqid}.{ff}"
        message: "Link sequencing file {wildcards.seqid} into temporary folder"
        params:
            fn = lambda wildcards: seqids[wildcards.seqid]
        shell:
            """
            ln -s {params.fn} {output}
            """

else:
    print(f"Unknown file format {config['fileformat']}. Exiting",
          file=sys.stderr)
    sys.exit(1)
