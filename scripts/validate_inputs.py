#!/usr/bin/env python

import argparse
import gzip
import re


CONTIG_PATTERN = re.compile(r"^##contig=<(?P<fields>.*)>$")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Validate that an alignment and VCF use contigs compatible with "
            "the indexed reference genome."
        )
    )
    parser.add_argument(
        "--alignment-header", required=True, help="SAM header produced by samtools view -H."
    )
    parser.add_argument("--variants", required=True, help="Input VCF or VCF.GZ file.")
    parser.add_argument("--reference-index", required=True, help="Target reference FASTA index.")
    args = parser.parse_args()

    reference_contigs = read_fai(args.reference_index)
    alignment_contigs = read_alignment_header(args.alignment_header)
    validate_alignment_contigs(alignment_contigs, reference_contigs)

    declared_contigs, record_contigs = read_vcf_contigs(args.variants)
    validate_vcf_contigs(
        declared_contigs, record_contigs, alignment_contigs, reference_contigs
    )


def read_fai(path):
    contigs = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid FASTA index line {line_number}: {line.rstrip()}")
            add_contig(contigs, fields[0], int(fields[1]), "reference FASTA index")
    if not contigs:
        raise ValueError("The reference FASTA index contains no contigs.")
    return contigs


def read_alignment_header(path):
    contigs = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("@SQ\t"):
                continue
            values = {}
            for field in line.rstrip("\n").split("\t")[1:]:
                key, separator, value = field.partition(":")
                if separator:
                    values[key] = value
            if "SN" not in values or "LN" not in values:
                raise ValueError(f"Malformed alignment @SQ line: {line.rstrip()}")
            add_contig(contigs, values["SN"], int(values["LN"]), "alignment header")
    if not contigs:
        raise ValueError("The alignment header contains no @SQ sequence records.")
    return contigs


def read_vcf_contigs(path):
    declared = {}
    records = set()
    open_vcf = gzip.open if str(path).lower().endswith(".gz") else open
    with open_vcf(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            match = CONTIG_PATTERN.match(line.rstrip("\n"))
            if match:
                attributes = parse_contig_attributes(match.group("fields"))
                if "ID" in attributes:
                    length = int(attributes["length"]) if "length" in attributes else None
                    add_contig(declared, attributes["ID"], length, "VCF header")
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t", 2)
            if len(fields) < 2:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")
            records.add(fields[0])
    if not records:
        raise ValueError("The VCF contains no variant records.")
    return declared, records


def parse_contig_attributes(text):
    attributes = {}
    for field in text.split(","):
        key, separator, value = field.partition("=")
        if separator:
            attributes[key] = value.strip('"')
    return attributes


def add_contig(contigs, name, length, source):
    if name in contigs and contigs[name] != length:
        raise ValueError(f"Conflicting lengths for contig '{name}' in the {source}.")
    contigs[name] = length


def validate_alignment_contigs(alignment_contigs, reference_contigs):
    missing = sorted(set(alignment_contigs) - set(reference_contigs))
    if missing:
        raise ValueError(
            "Alignment contigs absent from the reference FASTA: " + ", ".join(missing)
        )
    mismatched = sorted(
        name
        for name, length in alignment_contigs.items()
        if reference_contigs[name] != length
    )
    if mismatched:
        details = ", ".join(
            f"{name} (alignment={alignment_contigs[name]}, reference={reference_contigs[name]})"
            for name in mismatched
        )
        raise ValueError("Alignment/reference contig length mismatch: " + details)


def validate_vcf_contigs(
    declared_contigs, record_contigs, alignment_contigs, reference_contigs
):
    missing_reference = sorted(record_contigs - set(reference_contigs))
    if missing_reference:
        raise ValueError(
            "VCF record contigs absent from the reference FASTA: "
            + ", ".join(missing_reference)
        )
    missing_alignment = sorted(record_contigs - set(alignment_contigs))
    if missing_alignment:
        raise ValueError(
            "VCF record contigs absent from the alignment header: "
            + ", ".join(missing_alignment)
        )

    mismatched = []
    for name, length in declared_contigs.items():
        if length is None:
            continue
        if name in reference_contigs and reference_contigs[name] != length:
            mismatched.append(
                f"{name} (VCF={length}, reference={reference_contigs[name]})"
            )
        if name in alignment_contigs and alignment_contigs[name] != length:
            mismatched.append(
                f"{name} (VCF={length}, alignment={alignment_contigs[name]})"
            )
    if mismatched:
        raise ValueError("VCF contig length mismatch: " + ", ".join(sorted(mismatched)))


if __name__ == "__main__":
    main()
