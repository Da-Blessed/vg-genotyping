rule vg_longread_indexes:
    input:
        validated=LONG_VALIDATION_OK,
        gbz=VG_GBZ,
    output:
        dist=LONG_DIST,
        minimizer=LONG_MIN,
        zipcodes=LONG_ZIPCODES,
    params:
        prefix=lambda wildcards, output: str(Path(output.dist).with_suffix("")),
        tmp=f"{WORK_DIR}/longread/tmp/autoindex",
        target_mem=str(config.get("longread", {}).get("autoindex_target_mem", "64G")),
        outdir=lambda wildcards, output: str(Path(output.dist).parent),
        logdir=f"{LOG_DIR}/longread",
    threads:
        cfg_threads("longread_autoindex", 16)
    resources:
        mem_mb=cfg_mem("longread_autoindex", 64000),
    log:
        f"{LOG_DIR}/longread/index.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.tmp:q} {params.logdir:q}
        vg autoindex --workflow lr-giraffe \
            --prefix {params.prefix:q} --gbz {input.gbz:q} \
            --threads {threads} --tmp-dir {params.tmp:q} \
            --target-mem {params.target_mem:q} > {log:q} 2>&1
        """


rule vg_longread_map:
    input:
        validated=LONG_VALIDATION_OK,
        gbz=VG_GBZ,
        dist=LONG_DIST,
        minimizer=LONG_MIN,
        zipcodes=LONG_ZIPCODES,
        read=long_read,
    output:
        temp(LONG_GAM) if not config.get("longread", {}).get("keep_intermediates", False) else LONG_GAM,
    params:
        platform=long_platform,
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/longread/map",
    threads:
        cfg_threads("longread_giraffe", 16)
    resources:
        mem_mb=cfg_mem("longread_giraffe", 64000),
    log:
        f"{LOG_DIR}/longread/map/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg giraffe -t {threads} -b {params.platform:q} \
            -Z {input.gbz:q} -d {input.dist:q} -m {input.minimizer:q} \
            -z {input.zipcodes:q} -f {input.read:q} \
            > {output:q} 2> {log:q}
        """


rule vg_longread_pack:
    input:
        xg=VG_XG,
        gam=LONG_GAM,
    output:
        temp(LONG_PACK) if not config.get("longread", {}).get("keep_intermediates", False) else LONG_PACK,
    params:
        min_mapq=int(config.get("longread", {}).get("min_mapq", 5)),
        expected_coverage=int(config.get("longread", {}).get("expected_coverage", 64)),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/longread/pack",
    threads:
        cfg_threads("pack", 8)
    resources:
        mem_mb=cfg_mem("pack", 16000),
    log:
        f"{LOG_DIR}/longread/pack/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg pack -t {threads} -x {input.xg:q} -g {input.gam:q} \
            -Q {params.min_mapq} -e {params.expected_coverage} \
            -o {output:q} > {log:q} 2>&1
        """


rule vg_longread_genotype:
    input:
        xg=VG_XG,
        pack=LONG_PACK,
        panel=VG_PANEL_VCF,
        panel_tbi=VG_PANEL_TBI,
    output:
        LONG_RAW_VCF,
    params:
        ploidy=int(config.get("vg", {}).get("ploidy", 2)),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/longread/genotype",
    threads:
        cfg_threads("genotype", 8)
    resources:
        mem_mb=cfg_mem("genotype", 16000),
    log:
        f"{LOG_DIR}/longread/genotype/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        vg call {input.xg:q} -k {input.pack:q} -v {input.panel:q} \
            -s {wildcards.sample:q} -d {params.ploidy} -t {threads} \
            > {output:q} 2> {log:q}
        """


rule vg_longread_compress_vcf:
    input:
        LONG_RAW_VCF,
    output:
        vcf=LONG_SAMPLE_VCF,
        tbi=LONG_SAMPLE_TBI,
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        logdir=f"{LOG_DIR}/longread/vcf",
    threads:
        cfg_threads("vcf", 2)
    resources:
        mem_mb=cfg_mem("vcf", 2048),
    log:
        f"{LOG_DIR}/longread/vcf/{{sample}}.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        bcftools view --threads {threads} -Oz -o {output.vcf:q} {input:q} 2> {log:q}
        bcftools index --threads {threads} --tbi -o {output.tbi:q} {output.vcf:q} 2>> {log:q}
        """


rule vg_longread_merge_cohort:
    input:
        vcfs=expand(LONG_SAMPLE_VCF, sample=LONG_SAMPLES),
        tbis=expand(LONG_SAMPLE_TBI, sample=LONG_SAMPLES),
    output:
        vcf=LONG_COHORT_VCF,
        tbi=LONG_COHORT_TBI,
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        logdir=f"{LOG_DIR}/longread/cohort",
    threads:
        cfg_threads("vcf", 2)
    resources:
        mem_mb=cfg_mem("vcf", 2048),
    log:
        f"{LOG_DIR}/longread/cohort/merge.log",
    conda:
        "../../envs/vg.yaml"
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        bcftools merge --force-single -m none --threads {threads} \
            -Oz -o {output.vcf:q} {input.vcfs:q} 2> {log:q}
        bcftools index --threads {threads} --tbi -o {output.tbi:q} {output.vcf:q} 2>> {log:q}
        """
