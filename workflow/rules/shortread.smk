rule vg_shortread_indexes:
    input:
        validated=SHORT_VALIDATION_OK,
        gbz=VG_GBZ,
    output:
        dist=SHORT_DIST,
        minimizer=SHORT_MIN,
        zipcodes=SHORT_ZIPCODES,
    params:
        outdir=lambda wildcards, output: str(Path(output.dist).parent),
        prefix=lambda wildcards, output: str(Path(output.dist).with_suffix("")),
        logdir=f"{LOG_DIR}/shortread",
    threads:
        cfg_threads("shortread_index", 8)
    resources:
        mem_mb=cfg_mem("shortread_index", 32000),
    log:
        f"{LOG_DIR}/shortread/index.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg index -j {output.dist:q} {input.gbz:q} > {log:q} 2>&1
        vg minimizer -t {threads} -d {output.dist:q} \
            -o {output.minimizer:q} -z {output.zipcodes:q} \
            {input.gbz:q} >> {log:q} 2>&1
        """


rule vg_shortread_map:
    input:
        validated=SHORT_VALIDATION_OK,
        gbz=VG_GBZ,
        dist=SHORT_DIST,
        minimizer=SHORT_MIN,
        zipcodes=SHORT_ZIPCODES,
        read1=short_read1,
        read2=short_read2,
    output:
        temp(SHORT_GAM) if not config.get("shortread", {}).get("keep_intermediates", False) else SHORT_GAM,
    params:
        preset=str(config.get("shortread", {}).get("preset", "default")),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/shortread/map",
    threads:
        cfg_threads("shortread_giraffe", 16)
    resources:
        mem_mb=cfg_mem("shortread_giraffe", 32000),
    log:
        f"{LOG_DIR}/shortread/map/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg giraffe -t {threads} -b {params.preset:q} \
            -Z {input.gbz:q} -d {input.dist:q} -m {input.minimizer:q} \
            -z {input.zipcodes:q} -f {input.read1:q} -f {input.read2:q} \
            > {output:q} 2> {log:q}
        """


rule vg_shortread_pack:
    input:
        xg=VG_XG,
        gam=SHORT_GAM,
    output:
        temp(SHORT_PACK) if not config.get("shortread", {}).get("keep_intermediates", False) else SHORT_PACK,
    params:
        min_mapq=int(config.get("shortread", {}).get("min_mapq", 5)),
        expected_coverage=int(config.get("shortread", {}).get("expected_coverage", 128)),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/shortread/pack",
    threads:
        cfg_threads("pack", 8)
    resources:
        mem_mb=cfg_mem("pack", 16000),
    log:
        f"{LOG_DIR}/shortread/pack/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg pack -t {threads} -x {input.xg:q} -g {input.gam:q} \
            -Q {params.min_mapq} -e {params.expected_coverage} \
            -o {output:q} > {log:q} 2>&1
        """


rule vg_shortread_genotype:
    input:
        xg=VG_XG,
        pack=SHORT_PACK,
        panel=VG_PANEL_VCF,
        panel_tbi=VG_PANEL_TBI,
    output:
        SHORT_RAW_VCF,
    params:
        ploidy=int(config.get("vg", {}).get("ploidy", 2)),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/shortread/genotype",
    threads:
        cfg_threads("genotype", 8)
    resources:
        mem_mb=cfg_mem("genotype", 16000),
    log:
        f"{LOG_DIR}/shortread/genotype/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg call {input.xg:q} -k {input.pack:q} -v {input.panel:q} \
            -s {wildcards.sample:q} -d {params.ploidy} -t {threads} \
            > {output:q} 2> {log:q}
        """


rule vg_shortread_compress_vcf:
    input:
        SHORT_RAW_VCF,
    output:
        vcf=SHORT_SAMPLE_VCF,
        tbi=SHORT_SAMPLE_TBI,
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        logdir=f"{LOG_DIR}/shortread/vcf",
    threads:
        cfg_threads("vcf", 2)
    resources:
        mem_mb=cfg_mem("vcf", 2048),
    log:
        f"{LOG_DIR}/shortread/vcf/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        bcftools view --threads {threads} -Oz -o {output.vcf:q} {input:q} 2> {log:q}
        bcftools index --threads {threads} --tbi -o {output.tbi:q} {output.vcf:q} 2>> {log:q}
        """


rule vg_shortread_merge_cohort:
    input:
        vcfs=expand(SHORT_SAMPLE_VCF, sample=SHORT_SAMPLES),
        tbis=expand(SHORT_SAMPLE_TBI, sample=SHORT_SAMPLES),
    output:
        vcf=SHORT_COHORT_VCF,
        tbi=SHORT_COHORT_TBI,
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        logdir=f"{LOG_DIR}/shortread/cohort",
    threads:
        cfg_threads("vcf", 2)
    resources:
        mem_mb=cfg_mem("vcf", 2048),
    log:
        f"{LOG_DIR}/shortread/cohort/merge.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        bcftools merge --force-single -m none --threads {threads} \
            -Oz -o {output.vcf:q} {input.vcfs:q} 2> {log:q}
        bcftools index --threads {threads} --tbi -o {output.tbi:q} {output.vcf:q} 2>> {log:q}
        """
