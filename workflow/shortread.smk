configfile: "config/config.yaml"

include: "rules/common.smk"


rule all:
    input:
        SHORT_COHORT_VCF,
        SHORT_COHORT_TBI,


include: "rules/validation.smk"
include: "rules/graph.smk"
include: "rules/shortread.smk"
