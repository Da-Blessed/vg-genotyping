configfile: "config/config.yaml"

include: "rules/common.smk"


rule all:
    input:
        LONG_COHORT_VCF,
        LONG_COHORT_TBI,


include: "rules/validation.smk"
include: "rules/graph.smk"
include: "rules/longread.smk"
