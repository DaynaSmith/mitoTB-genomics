#!/usr/bin/env python3
"""
Build the unified PE metadata table for Analysis 2 (L vs non-L DE) and
Analysis 3 (TB case/control DE).

Joins:
  - metadata.csv (clinical/demographic data, FILE column transformed to
    match your sample ID scheme)
  - haplogroup_qc_summary.tsv (haplogroup, haplogroup_quality — already
    covers all 520 samples including PE)
  - master_excluded_samples.txt (samples dropped from every analysis)

Usage:
    python 29_build_pe_metadata.py \
        --metadata metadata.csv \
        --haplogroup haplogroup_qc_summary.tsv \
        --excluded master_excluded_samples.txt \
        --bam-manifest bam_manifest.txt \
        --out pe_metadata_unified.tsv
"""

import argparse
import csv


def transform_file_id(file_str):
    """Convert metadata.csv's FILE column format to the project's sample_id scheme."""
    if not file_str:
        return None
    parts = file_str.split("_")
    if len(parts) != 5:
        return None
    date, run, flowcell_full, lane, tail = parts
    if not flowcell_full.startswith("FC") or len(flowcell_full) < 7:
        return None
    flowcell_short = flowcell_full[3:-4]
    if "-" not in tail:
        return None
    idx = tail.split("-")[-1]
    return f"{date}_{run}.{flowcell_short}.{lane}.{idx}"


def is_missing(val):
    """Treat None, empty string, and common missing-value markers as missing —
    not just Python falsiness, since the source CSV may encode blanks as the
    literal string 'NA' (non-empty, so 'not val' alone would miss it)."""
    if val is None:
        return True
    return str(val).strip().upper() in ("", "NA", "N/A", "NAN", "NULL")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--haplogroup", required=True)
    p.add_argument("--excluded", required=True)
    p.add_argument("--bam-manifest", required=True)
    p.add_argument("--out", default="pe_metadata_unified.tsv")
    args = p.parse_args()

    # ---- Load master exclusion list ----
    excluded = set()
    with open(args.excluded) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            excluded.add(line.split("\t")[0])

    # ---- Load PE sample list from bam_manifest ----
    pe_samples = set()
    with open(args.bam_manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == "PE":
                pe_samples.add(parts[0])

    # ---- Load haplogroup data ----
    haplogroup = {}
    with open(args.haplogroup) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            haplogroup[row["sample"]] = (row["haplogroup"], row["quality"])

    # ---- Load and transform clinical metadata ----
    clinical = {}
    with open(args.metadata) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = transform_file_id(row.get("FILE"))
            if sid is None:
                continue
            clinical[sid] = row

    # ---- Build unified table ----
    rows = []
    no_clinical_data = []
    no_haplogroup = []

    for sample_id in sorted(pe_samples):
        if sample_id in excluded:
            continue  # master-excluded, skip entirely

        clin = clinical.get(sample_id)
        hap = haplogroup.get(sample_id)

        if hap is None:
            no_haplogroup.append(sample_id)
            haplogroup_val, quality_val = "NA", "NA"
        else:
            haplogroup_val, quality_val = hap

        if clin is None or is_missing(clin.get("STUDY_GROUP")):
            no_clinical_data.append(sample_id)
            condition, sex, age = "NA", "NA", "NA"
        else:
            condition = clin.get("STUDY_GROUP") or "NA"
            sex = clin.get("GENDER") or "NA"
            age = clin.get("AGE") or "NA"

        rows.append({
            "sample_id": sample_id,
            "cohort": "PE",
            "condition": condition,
            "sex": sex,
            "age": age,
            "haplogroup": haplogroup_val,
            "haplogroup_quality": quality_val,
        })

    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "sample_id", "cohort", "condition", "sex", "age",
            "haplogroup", "haplogroup_quality"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_usable = sum(1 for r in rows if r["condition"] != "NA")
    print(f"Wrote {len(rows)} PE samples to {args.out}")
    print(f"  {n_usable} have complete condition/sex/age data — usable for Analysis 2/3")
    print(f"  {len(no_clinical_data)} have NO clinical metadata match "
          f"(sequenced but not linked to phenotype records)")
    print(f"  {len(no_haplogroup)} have no haplogroup call (unexpected — investigate if >0)")

    if no_clinical_data:
        with open("pe_samples_missing_clinical_data.txt", "w") as f:
            for s in no_clinical_data:
                f.write(s + "\n")
        print(f"\n  Full list written to pe_samples_missing_clinical_data.txt")


if __name__ == "__main__":
    main()
