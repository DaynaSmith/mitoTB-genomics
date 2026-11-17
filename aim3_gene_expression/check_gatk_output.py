#!/usr/bin/env python3
"""
Sweep gatk_out/ and summarize every sample's variant calling result.
Flags: missing output, empty VCFs, unusually low/high variant counts.

Usage:
    python 20_check_gatk_output.py --manifest bam_manifest.txt \
        --outdir gatk_out --out gatk_qc_summary.tsv
"""

import argparse
import gzip
import os


def count_variants(vcf_path):
    """Count non-header lines (variant records) in a (possibly gzipped) VCF."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    n = 0
    n_pass = 0
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n += 1
            fields = line.split("\t")
            if len(fields) > 6 and (fields[6] == "PASS" or fields[6] == "."):
                n_pass += 1
    return n, n_pass


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="bam_manifest.txt used for the run")
    p.add_argument("--outdir", default="gatk_out")
    p.add_argument("--out", default="gatk_qc_summary.tsv")
    args = p.parse_args()

    samples = []
    with open(args.manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            sample_id = parts[0]
            lib_type = parts[2] if len(parts) > 2 else "?"
            samples.append((sample_id, lib_type))

    rows = []
    for sample_id, lib_type in samples:
        vcf_path = os.path.join(args.outdir, sample_id, f"{sample_id}.filtered.vcf.gz")
        if not os.path.exists(vcf_path):
            rows.append({"sample": sample_id, "library_type": lib_type,
                         "status": "MISSING", "n_variants": 0, "n_pass": 0})
            continue

        try:
            n_total, n_pass = count_variants(vcf_path)
        except Exception as e:
            rows.append({"sample": sample_id, "library_type": lib_type,
                         "status": f"READ_ERROR: {e}", "n_variants": 0, "n_pass": 0})
            continue

        if n_total == 0:
            status = "EMPTY"
        else:
            status = "OK"

        rows.append({"sample": sample_id, "library_type": lib_type,
                     "status": status, "n_variants": n_total, "n_pass": n_pass})

    with open(args.out, "w") as f:
        f.write("sample\tlibrary_type\tstatus\tn_variants\tn_pass\n")
        for r in rows:
            f.write(f"{r['sample']}\t{r['library_type']}\t{r['status']}\t"
                     f"{r['n_variants']}\t{r['n_pass']}\n")

    n_total = len(rows)
    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_missing = sum(1 for r in rows if r["status"] == "MISSING")
    n_empty = sum(1 for r in rows if r["status"] == "EMPTY")
    n_error = sum(1 for r in rows if r["status"].startswith("READ_ERROR"))

    print(f"Wrote {n_total} samples to {args.out}")
    print(f"  OK:      {n_ok}")
    print(f"  MISSING: {n_missing}")
    print(f"  EMPTY:   {n_empty}  (VCF exists but has zero variant records — worth checking)")
    print(f"  ERROR:   {n_error}")

    if n_ok:
        pass_counts = [r["n_pass"] for r in rows if r["status"] == "OK"]
        print(f"\nPASS variant count across OK samples: "
              f"min={min(pass_counts)}, max={max(pass_counts)}, "
              f"mean={sum(pass_counts)/len(pass_counts):.1f}")

    if n_missing or n_empty or n_error:
        print("\nSamples needing attention:")
        for r in rows:
            if r["status"] != "OK":
                print(f"  {r['sample']} ({r['library_type']}): {r['status']}")


if __name__ == "__main__":
    main()
