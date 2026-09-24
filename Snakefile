configfile: "config/config.yaml"

include: "workflow/rules/common.smk"


rule all:
    input:
        SHORT_COHORT_VCF,
        SHORT_COHORT_TBI,
        LONG_COHORT_VCF,
        LONG_COHORT_TBI,


include: "workflow/rules/validation.smk"
include: "workflow/rules/graph.smk"
include: "workflow/rules/shortread.smk"
include: "workflow/rules/longread.smk"
