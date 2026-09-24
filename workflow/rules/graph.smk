rule vg_stage_reference:
    input:
        validated=PANEL_VALIDATION_OK,
        reference=REFERENCE,
    output:
        reference=VG_REFERENCE,
        fai=VG_REFERENCE_FAI,
    params:
        outdir=lambda wildcards, output: str(Path(output.reference).parent),
        logdir=f"{LOG_DIR}/graph",
    log:
        f"{LOG_DIR}/graph/stage_reference.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        ln -sfn "$(realpath {input.reference:q})" {output.reference:q}
        samtools faidx {output.reference:q} > {log:q} 2>&1
        """


rule vg_stage_panel:
    input:
        validated=PANEL_VALIDATION_OK,
        vcf=PANEL_VCF,
    output:
        vcf=VG_PANEL_VCF,
        tbi=VG_PANEL_TBI,
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        logdir=f"{LOG_DIR}/graph",
    threads:
        cfg_threads("vcf", 2)
    log:
        f"{LOG_DIR}/graph/stage_panel.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        bcftools view --threads {threads} -Oz -o {output.vcf:q} {input.vcf:q} 2> {log:q}
        bcftools index --threads {threads} --tbi -o {output.tbi:q} {output.vcf:q} 2>> {log:q}
        """


rule vg_construct_graph:
    input:
        reference=VG_REFERENCE,
        fai=VG_REFERENCE_FAI,
        vcf=VG_PANEL_VCF,
        tbi=VG_PANEL_TBI,
    output:
        VG_GRAPH,
    params:
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/graph",
    threads:
        cfg_threads("graph_construct", 8)
    resources:
        mem_mb=cfg_mem("graph_construct", 32000),
    log:
        f"{LOG_DIR}/graph/construct.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg construct -a -t {threads} -r {input.reference:q} -v {input.vcf:q} \
            > {output:q} 2> {log:q}
        """


rule vg_xg_index:
    input:
        VG_GRAPH,
    output:
        VG_XG,
    params:
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/graph",
    threads:
        cfg_threads("graph_index", 8)
    resources:
        mem_mb=cfg_mem("graph_index", 32000),
    log:
        f"{LOG_DIR}/graph/xg.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg index -L -t {threads} -x {output:q} {input:q} > {log:q} 2>&1
        """


rule vg_gbwt_index:
    input:
        xg=VG_XG,
        vcf=VG_PANEL_VCF,
        tbi=VG_PANEL_TBI,
    output:
        VG_GBWT,
    params:
        tmp=f"{WORK_DIR}/graph/tmp/gbwt",
        logdir=f"{LOG_DIR}/graph",
    threads:
        cfg_threads("graph_index", 8)
    resources:
        mem_mb=cfg_mem("graph_index", 32000),
    log:
        f"{LOG_DIR}/graph/gbwt.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.tmp:q} {params.logdir:q}
        vg gbwt -x {input.xg:q} -v {input.vcf:q} -o {output:q} \
            --num-jobs {threads} -d {params.tmp:q} > {log:q} 2>&1
        """


rule vg_gbz_index:
    input:
        xg=VG_XG,
        gbwt=VG_GBWT,
    output:
        VG_GBZ,
    params:
        logdir=f"{LOG_DIR}/graph",
    threads:
        cfg_threads("graph_index", 8)
    resources:
        mem_mb=cfg_mem("graph_index", 32000),
    log:
        f"{LOG_DIR}/graph/gbz.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.logdir:q}
        vg gbwt -x {input.xg:q} -g {output:q} --gbz-format {input.gbwt:q} \
            > {log:q} 2>&1
        """
