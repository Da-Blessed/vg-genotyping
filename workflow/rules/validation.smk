rule validate_panel:
    input:
        reference=REFERENCE,
        vcf=PANEL_VCF,
    output:
        touch(PANEL_VALIDATION_OK),
    params:
        script=str(PROJECT_ROOT / "scripts" / "validate_inputs.py"),
        strict="--strict" if config.get("validation", {}).get("strict", True) else "",
        allow_missing=(
            "--allow-missing-genotypes"
            if config.get("validation", {}).get("allow_missing_genotypes", False)
            else ""
        ),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/validation",
    log:
        f"{LOG_DIR}/validation/panel.log",
    conda:
        "../../envs/common.yaml"
    resources:
        mem_mb=cfg_mem("validation", 1024),
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        python {params.script:q} \
            --reference {input.reference:q} \
            --vcf {input.vcf:q} \
            --panel-only {params.strict} {params.allow_missing} > {log:q} 2>&1
        """


rule validate_shortread_inputs:
    input:
        reference=REFERENCE,
        vcf=PANEL_VCF,
        panel=PANEL_VALIDATION_OK,
        samples=SHORT_SAMPLES_TSV,
        reads=lambda wildcards: [
            path
            for sample in SHORT_SAMPLES
            for path in (
                SHORT_RECORDS[sample]["read1"],
                SHORT_RECORDS[sample]["read2"],
            )
        ],
    output:
        touch(SHORT_VALIDATION_OK),
    params:
        script=str(PROJECT_ROOT / "scripts" / "validate_inputs.py"),
        strict="--strict" if config.get("validation", {}).get("strict", True) else "",
        allow_missing=(
            "--allow-missing-genotypes"
            if config.get("validation", {}).get("allow_missing_genotypes", False)
            else ""
        ),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/validation",
    log:
        f"{LOG_DIR}/validation/shortread.log",
    conda:
        "../../envs/common.yaml"
    resources:
        mem_mb=cfg_mem("validation", 1024),
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        python {params.script:q} \
            --reference {input.reference:q} \
            --vcf {input.vcf:q} \
            --samples {input.samples:q} \
            --read-type short {params.strict} {params.allow_missing} > {log:q} 2>&1
        """


rule validate_longread_inputs:
    input:
        reference=REFERENCE,
        vcf=PANEL_VCF,
        panel=PANEL_VALIDATION_OK,
        samples=LONG_SAMPLES_TSV,
        reads=lambda wildcards: [LONG_RECORDS[sample]["read"] for sample in LONG_SAMPLES],
    output:
        touch(LONG_VALIDATION_OK),
    params:
        script=str(PROJECT_ROOT / "scripts" / "validate_inputs.py"),
        strict="--strict" if config.get("validation", {}).get("strict", True) else "",
        allow_missing=(
            "--allow-missing-genotypes"
            if config.get("validation", {}).get("allow_missing_genotypes", False)
            else ""
        ),
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        logdir=f"{LOG_DIR}/validation",
    log:
        f"{LOG_DIR}/validation/longread.log",
    conda:
        "../../envs/common.yaml"
    resources:
        mem_mb=cfg_mem("validation", 1024),
    shell:
        r"""
        mkdir -p {params.outdir:q} {params.logdir:q}
        python {params.script:q} \
            --reference {input.reference:q} \
            --vcf {input.vcf:q} \
            --samples {input.samples:q} \
            --read-type long {params.strict} {params.allow_missing} > {log:q} 2>&1
        """
