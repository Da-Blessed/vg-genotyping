import csv
import re
from pathlib import Path


PROJECT_ROOT = Path(workflow.basedir).resolve()
if not (PROJECT_ROOT / "config").is_dir():
    PROJECT_ROOT = PROJECT_ROOT.parent


def required_config(key):
    value = config.get(key)
    if value is None or str(value).strip() == "":
        raise WorkflowError(f"Missing required configuration key: {key}")
    return str(value)


REFERENCE = required_config("reference")
PANEL_VCF = required_config("panel_vcf")
SHORT_SAMPLES_TSV = required_config("shortread_samples")
LONG_SAMPLES_TSV = required_config("longread_samples")
OUTPUT_DIR = str(config.get("output_dir", "results")).rstrip("/")
WORK_DIR = str(config.get("work_dir", "work")).rstrip("/")
LOG_DIR = str(config.get("log_dir", "logs")).rstrip("/")

SAFE_SAMPLE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
LONGREAD_PLATFORMS = {"hifi", "r10"}


def load_sample_sheet(path, required_columns):
    records = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", *required_columns}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            columns = ", ".join(["sample", *required_columns])
            raise WorkflowError(f"{path} must contain tab-separated columns: {columns}")
        for line_number, row in enumerate(reader, start=2):
            record = {column: (row.get(column) or "").strip() for column in required}
            sample = record["sample"]
            if not SAFE_SAMPLE.fullmatch(sample):
                raise WorkflowError(f"Unsafe sample name at {path}:{line_number}: {sample!r}")
            if any(not record[column] for column in required):
                raise WorkflowError(f"Empty value in {path}, line {line_number}")
            if sample in records:
                raise WorkflowError(f"Duplicate sample '{sample}' in {path}")
            records[sample] = record
    if not records:
        raise WorkflowError(f"No samples found in {path}")
    return records


SHORT_RECORDS = load_sample_sheet(SHORT_SAMPLES_TSV, ["read1", "read2"])
LONG_RECORDS = load_sample_sheet(LONG_SAMPLES_TSV, ["read", "platform"])
for sample, record in LONG_RECORDS.items():
    record["platform"] = record["platform"].lower()
    if record["platform"] not in LONGREAD_PLATFORMS:
        raise WorkflowError(
            f"Unsupported long-read platform for {sample}: {record['platform']}; "
            "use hifi or r10"
        )

SHORT_SAMPLES = list(SHORT_RECORDS)
LONG_SAMPLES = list(LONG_RECORDS)


def short_read1(wildcards):
    return SHORT_RECORDS[wildcards.sample]["read1"]


def short_read2(wildcards):
    return SHORT_RECORDS[wildcards.sample]["read2"]


def long_read(wildcards):
    return LONG_RECORDS[wildcards.sample]["read"]


def long_platform(wildcards):
    return LONG_RECORDS[wildcards.sample]["platform"]


def cfg_threads(name, default):
    return max(1, int(config.get("threads", {}).get(name, default)))


def cfg_mem(name, default):
    return max(1, int(config.get("resources", {}).get(name, default)))


PANEL_VALIDATION_OK = f"{WORK_DIR}/validation/panel.ok"
SHORT_VALIDATION_OK = f"{WORK_DIR}/validation/shortread.ok"
LONG_VALIDATION_OK = f"{WORK_DIR}/validation/longread.ok"

VG_REFERENCE = f"{WORK_DIR}/graph/reference/reference.fa"
VG_REFERENCE_FAI = f"{WORK_DIR}/graph/reference/reference.fa.fai"
VG_PANEL_VCF = f"{WORK_DIR}/graph/panel/panel.vcf.gz"
VG_PANEL_TBI = f"{WORK_DIR}/graph/panel/panel.vcf.gz.tbi"
VG_GRAPH = f"{WORK_DIR}/graph/index/panel.vg"
VG_XG = f"{WORK_DIR}/graph/index/panel.xg"
VG_GBWT = f"{WORK_DIR}/graph/index/panel.gbwt"
VG_GBZ = f"{WORK_DIR}/graph/index/panel.gbz"

SHORT_DIST = f"{WORK_DIR}/shortread/index/panel.dist"
SHORT_MIN = f"{WORK_DIR}/shortread/index/panel.withzip.min"
SHORT_ZIPCODES = f"{WORK_DIR}/shortread/index/panel.zipcodes"
SHORT_GAM = f"{WORK_DIR}/shortread/mapping/{{sample}}.gam"
SHORT_PACK = f"{WORK_DIR}/shortread/pack/{{sample}}.pack"
SHORT_RAW_VCF = f"{WORK_DIR}/shortread/genotypes/{{sample}}.vcf"
SHORT_SAMPLE_VCF = f"{OUTPUT_DIR}/shortread/vcf/{{sample}}.vcf.gz"
SHORT_SAMPLE_TBI = f"{OUTPUT_DIR}/shortread/vcf/{{sample}}.vcf.gz.tbi"
SHORT_COHORT_VCF = f"{OUTPUT_DIR}/shortread/cohort/vg.shortread.cohort.vcf.gz"
SHORT_COHORT_TBI = f"{OUTPUT_DIR}/shortread/cohort/vg.shortread.cohort.vcf.gz.tbi"

LONG_INDEX_PREFIX = f"{WORK_DIR}/longread/index/panel"
LONG_DIST = f"{LONG_INDEX_PREFIX}.dist"
LONG_MIN = f"{LONG_INDEX_PREFIX}.longread.withzip.min"
LONG_ZIPCODES = f"{LONG_INDEX_PREFIX}.longread.zipcodes"
LONG_GAM = f"{WORK_DIR}/longread/mapping/{{sample}}.gam"
LONG_PACK = f"{WORK_DIR}/longread/pack/{{sample}}.pack"
LONG_RAW_VCF = f"{WORK_DIR}/longread/genotypes/{{sample}}.vcf"
LONG_SAMPLE_VCF = f"{OUTPUT_DIR}/longread/vcf/{{sample}}.vcf.gz"
LONG_SAMPLE_TBI = f"{OUTPUT_DIR}/longread/vcf/{{sample}}.vcf.gz.tbi"
LONG_COHORT_VCF = f"{OUTPUT_DIR}/longread/cohort/vg.longread.cohort.vcf.gz"
LONG_COHORT_TBI = f"{OUTPUT_DIR}/longread/cohort/vg.longread.cohort.vcf.gz.tbi"
