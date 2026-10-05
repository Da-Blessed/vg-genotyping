#!/usr/bin/env python3
"""Validate panel and short- or long-read sample inputs for vg genotyping."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
SAMPLE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
LONGREAD_PLATFORMS = {"hifi", "r10"}


class ValidationError(RuntimeError):
    pass


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if is_gzip(path):
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open() as handle:
            yield handle


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise ValidationError(f"{label} does not exist or is not a file: {path}")
    if path.stat().st_size == 0:
        raise ValidationError(f"{label} is empty: {path}")


def read_fasta_contigs(path: Path) -> dict[str, int]:
    require_file(path, "Reference FASTA")
    if is_gzip(path):
        raise ValidationError("The reference FASTA must be uncompressed")
    contigs: dict[str, int] = {}
    current: str | None = None
    length = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    contigs[current] = length
                current = line[1:].strip().split()[0]
                if not current:
                    raise ValidationError(f"Empty FASTA name at line {line_number}")
                if current in contigs:
                    raise ValidationError(f"Duplicate FASTA sequence name: {current}")
                length = 0
            else:
                if current is None:
                    raise ValidationError("FASTA sequence appears before the first header")
                length += len(line.strip())
    if current is not None:
        contigs[current] = length
    if not contigs:
        raise ValidationError("Reference FASTA contains no sequences")
    return contigs


def validate_read(path: Path, label: str) -> None:
    require_file(path, label)
    with open_text(path) as handle:
        first = handle.readline()
        if first.startswith("@"):
            sequence = handle.readline().rstrip("\n\r")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n\r")
            if not sequence or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValidationError(f"Malformed first FASTQ record in {path}")
        elif first.startswith(">"):
            if not handle.readline().strip():
                raise ValidationError(f"Malformed first FASTA record in {path}")
        else:
            raise ValidationError(f"Reads must be FASTA or FASTQ: {path}")


def validate_samples(path: Path, read_type: str) -> int:
    require_file(path, "Sample sheet")
    data_columns = ["read1", "read2"] if read_type == "short" else ["read", "platform"]
    required = {"sample", *data_columns}
    seen: set[str] = set()
    count = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValidationError(
                f"Sample sheet columns must be: {', '.join(['sample', *data_columns])}"
            )
        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            if not SAMPLE_NAME.fullmatch(sample):
                raise ValidationError(f"Invalid sample name at {path}:{line_number}: {sample!r}")
            if sample in seen:
                raise ValidationError(f"Duplicate sample name: {sample}")
            seen.add(sample)
            if read_type == "short":
                for column in ("read1", "read2"):
                    value = (row.get(column) or "").strip()
                    if not value:
                        raise ValidationError(f"Empty {column} at {path}:{line_number}")
                    validate_read(Path(value), f"{sample} {column}")
            else:
                read = (row.get("read") or "").strip()
                platform = (row.get("platform") or "").strip().lower()
                if platform not in LONGREAD_PLATFORMS:
                    raise ValidationError(
                        f"Long-read platform must be hifi or r10 at {path}:{line_number}"
                    )
                validate_read(Path(read), f"{sample} read")
            count += 1
    if count == 0:
        raise ValidationError("Sample sheet contains no samples")
    return count


def validate_gt(
    gt: str,
    sample: str,
    locus: str,
    strict: bool,
    allow_missing: bool,
) -> None:
    if not strict:
        return
    alleles = re.split(r"[|/]", gt)
    has_missing = any(allele == "." for allele in alleles)
    if not alleles or (has_missing and not allow_missing):
        raise ValidationError(f"Missing panel genotype for {sample} at {locus}")
    # A completely missing genotype carries no phase information. Partially
    # called and fully called diploid genotypes must remain phased.
    if (
        len(alleles) > 1
        and "|" not in gt
        and not all(allele == "." for allele in alleles)
    ):
        raise ValidationError(f"Unphased panel genotype for {sample} at {locus}: {gt}")
    if not all(
        allele.isdigit() or (allow_missing and allele == ".") for allele in alleles
    ):
        raise ValidationError(f"Invalid panel genotype for {sample} at {locus}: {gt}")


def validate_vcf(
    path: Path,
    contigs: dict[str, int],
    strict: bool,
    allow_missing: bool,
) -> tuple[int, int]:
    require_file(path, "Panel VCF")
    if is_gzip(path):
        raise ValidationError("Provide the panel as an uncompressed VCF")
    samples: list[str] = []
    variants = 0
    header_seen = False
    contig_order = {name: index for index, name in enumerate(contigs)}
    last_contig_index = -1
    last_contig: str | None = None
    last_position = 0
    previous_end: dict[str, int] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n\r").split("\t")
                if len(columns) < 10:
                    raise ValidationError("Panel VCF must contain at least one sample")
                samples = columns[9:]
                header_seen = True
                continue
            if line.startswith("#"):
                continue
            if not header_seen:
                raise ValidationError("VCF records appear before the #CHROM header")
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9 + len(samples):
                raise ValidationError(f"Malformed VCF record at line {line_number}")
            chrom, pos_text, _identifier, ref, alt = fields[:5]
            locus = f"{chrom}:{pos_text}"
            if chrom not in contigs:
                raise ValidationError(f"VCF contig is absent from reference: {chrom}")
            try:
                pos = int(pos_text)
            except ValueError as error:
                raise ValidationError(f"Invalid VCF position at {locus}") from error
            if pos < 1 or pos + len(ref) - 1 > contigs[chrom]:
                raise ValidationError(f"VCF record lies outside the reference at {locus}")
            alleles = [ref, *alt.split(",")]
            if any(not allele or allele in {".", "*"} or not DNA.fullmatch(allele) for allele in alleles):
                raise ValidationError(f"Explicit sequence-resolved alleles are required at {locus}")
            current_contig_index = contig_order[chrom]
            if current_contig_index < last_contig_index:
                raise ValidationError(f"VCF contigs are not in reference order at {locus}")
            if chrom == last_contig and pos < last_position:
                raise ValidationError(f"VCF positions are not sorted at {locus}")
            last_contig_index = current_contig_index
            last_contig = chrom
            last_position = pos
            end = pos + len(ref) - 1
            if pos <= previous_end.get(chrom, 0):
                raise ValidationError(f"Overlapping VCF records at {locus}")
            previous_end[chrom] = end
            formats = fields[8].split(":")
            if "GT" not in formats:
                raise ValidationError(f"GT is absent from FORMAT at {locus}")
            gt_index = formats.index("GT")
            for sample, sample_field in zip(samples, fields[9:]):
                values = sample_field.split(":")
                if gt_index >= len(values):
                    raise ValidationError(f"GT is absent for {sample} at {locus}")
                validate_gt(values[gt_index], sample, locus, strict, allow_missing)
            variants += 1
    if not header_seen or variants == 0:
        raise ValidationError("Panel VCF must contain a #CHROM header and variants")
    return len(samples), variants


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--samples", type=Path)
    parser.add_argument("--read-type", choices=("short", "long"))
    parser.add_argument("--panel-only", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--allow-missing-genotypes", action="store_true")
    args = parser.parse_args()
    try:
        contigs = read_fasta_contigs(args.reference)
        panel_samples, variants = validate_vcf(
            args.vcf,
            contigs,
            args.strict,
            args.allow_missing_genotypes,
        )
        if args.panel_only:
            target_samples = None
        else:
            if args.samples is None or args.read_type is None:
                parser.error("--samples and --read-type are required unless --panel-only is used")
            target_samples = validate_samples(args.samples, args.read_type)
    except (OSError, UnicodeError, gzip.BadGzipFile, ValidationError) as error:
        parser.error(str(error))
    summary = (
        f"Validated {len(contigs)} contigs, {variants} variants, "
        f"and {panel_samples} panel samples"
    )
    if target_samples is not None:
        summary += f", plus {target_samples} {args.read_type}-read target samples"
    print(summary + ".")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
