#!/usr/bin/env python3
"""
Build the unified SE metadata table for Analysis 2 (L vs non-L DE) and
Analysis 3 (TB case/control DE).

Joins:
  - se_clinical_metadata.csv (exported from "Records for analysis" sheet of
    ALL_COMBINED_RNASEQ_SAMPLES.xlsx). Sample Key for RNASEQ matches your
    sample IDs after stripping the trailing _X suffix
    (e.g. RSEQ006_E -> RSEQ006).
  - haplogroup_qc_summary.tsv (haplogroup, haplogroup_quality)
  - master_excluded_samples.txt (samples dropped from every analysis)

Case/Control classification (per project decision):
  Case:    Disease Category in {"TB only", "TB-only"}
  Control: Disease Category == "Healthy Control"
  Excluded from case/control comparison (but NOT from Analysis 2/haplogroup
  work): "DM only", "TB-DM", anything else

Usage:
    python 34_build_se_metadata.py \
        --clinical se_clinical_metadata.csv \
        --haplogroup haplogroup_qc_summary.tsv \
        --excluded master_excluded_samples.txt \
        --bam-manifest bam_manifest.txt \
        --out se_metadata_unified.tsv
"""

import argparse
import csv


def is_missing(val):
    if val is None:
        return True
    return str(val).strip().upper() in ("", "NA", "N/A", "NAN", "NULL")


def classify_condition(disease_category):
    if is_missing(disease_category):
        return "NA"
    dc = disease_category.strip()
    if dc in ("TB only", "TB-only"):
        return "Case"
    if dc == "Healthy Control":
        return "Control"
    return "NA"   # DM only, TB-DM, or anything else -> excluded from case/control comparison


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--clinical", required=True)
    p.add_argument("--haplogroup", required=True)
    p.add_argument("--excluded", required=True)
    p.add_argument("--bam-manifest", required=True)
    p.add_argument("--out", default="se_metadata_unified.tsv")
    args = p.parse_args()

    # ---- Load master exclusion list ----
    excluded = set()
    with open(args.excluded) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            excluded.add(line.split("\t")[0])

    # ---- Load SE sample list from bam_manifest ----
    se_samples = set()
    with open(args.bam_manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == "SE":
                se_samples.add(parts[0])

    # ---- Load haplogroup data ----
    haplogroup = {}
    with open(args.haplogroup) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            haplogroup[row["sample"]] = (row["haplogroup"], row["quality"])

    # ---- Load and index clinical metadata by base ID (strip trailing _X) ----
    # encoding="utf-8-sig" strips a BOM if Excel's CSV export added one — a
    # very common cause of every-row-fails-to-match when the header text
    # looks correct but silently isn't.
    clinical = {}
    with open(args.clinical, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [fn.strip() for fn in reader.fieldnames]  # strip stray whitespace too
        for row in reader:
            key = row.get("Sample Key for RNASEQ")
            if is_missing(key):
                continue
            clinical[key.strip()] = row

    # ---- Build unified table ----
    rows = []
    no_clinical_data = []
    no_haplogroup = []
    excluded_disease_category = []

    for sample_id in sorted(se_samples):
        if sample_id in excluded:
            continue

        base_id = sample_id.split("_")[0]
        clin = clinical.get(base_id)
        hap = haplogroup.get(sample_id)

        if hap is None:
            no_haplogroup.append(sample_id)
            haplogroup_val, quality_val = "NA", "NA"
        else:
            haplogroup_val, quality_val = hap

        if clin is None:
            no_clinical_data.append(sample_id)
            condition, sex, age, disease_category = "NA", "NA", "NA", "NA"
        else:
            disease_category = clin.get("Disease Category") or "NA"
            condition = classify_condition(disease_category)
            if condition == "NA" and not is_missing(disease_category):
                excluded_disease_category.append((sample_id, disease_category))
            sex = clin.get("Gender") or "NA"
            age = clin.get("Age at TB diagnosis") or "NA"

        rows.append({
            "sample_id": sample_id,
            "cohort": "SE",
            "condition": condition,
            "sex": sex,
            "age": age,
            "haplogroup": haplogroup_val,
            "haplogroup_quality": quality_val,
            "disease_category_raw": disease_category,
        })

    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "sample_id", "cohort", "condition", "sex", "age",
            "haplogroup", "haplogroup_quality", "disease_category_raw"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_case = sum(1 for r in rows if r["condition"] == "Case")
    n_control = sum(1 for r in rows if r["condition"] == "Control")
    print(f"Wrote {len(rows)} SE samples to {args.out}")
    print(f"  Case: {n_case}, Control: {n_control} — usable for Analysis 3 (TB case/control)")
    print(f"  All {len(rows)} have haplogroup data — usable for Analysis 2 (L vs non-L)")

    if no_clinical_data:
        print(f"\n{len(no_clinical_data)} samples with NO clinical metadata match:")
        print("  " + ", ".join(no_clinical_data[:10]) + (" ..." if len(no_clinical_data) > 10 else ""))
    if no_haplogroup:
        print(f"\n{len(no_haplogroup)} samples with no haplogroup call (unexpected — investigate if >0)")
    if excluded_disease_category:
        print(f"\n{len(excluded_disease_category)} samples excluded from case/control due to Disease Category "
              f"(DM only / TB-DM / other):")
        cats = {}
        for _, dc in excluded_disease_category:
            cats[dc] = cats.get(dc, 0) + 1
        for dc, n in cats.items():
            print(f"    {dc}: {n}")


if __name__ == "__main__":
    main()
