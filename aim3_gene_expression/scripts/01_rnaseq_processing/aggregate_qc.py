#!/usr/bin/env python3
"""
Aggregate per-sample mtDNA coverage summaries (from 04_mtdna_coverage.pbs)
and STAR Log.final.out mapping stats into one master QC table.

Usage:
    python 05_aggregate_qc.py \
        --coverage-dir coverage_summaries \
        --pe-dir aligned_PE \
        --se-dir aligned_SE \
        --out qc_summary.tsv
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd


def parse_star_log(log_path):
    """Pull key stats out of a STAR Log.final.out file."""
    stats = {}
    with open(log_path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = line.split("|")
            key = key.strip()
            val = val.strip().replace("%", "")
            if key == "Number of input reads":
                stats["input_reads"] = int(val)
            elif key == "Uniquely mapped reads %":
                stats["uniquely_mapped_pct"] = float(val)
            elif key == "% of reads mapped to too many loci":
                stats["multimapped_pct"] = float(val)
            elif key == "% of reads unmapped: too short":
                stats["unmapped_too_short_pct"] = float(val)
    return stats


def collect_star_stats(aligned_dir, library_type):
    rows = []
    for log_path in glob.glob(os.path.join(aligned_dir, "*", "*_Log.final.out")):
        sample = os.path.basename(log_path).replace("_Log.final.out", "")
        stats = parse_star_log(log_path)
        stats["sample"] = sample
        stats["library_type"] = library_type
        rows.append(stats)
    return rows


def collect_coverage_stats(coverage_dir):
    rows = []
    for cov_path in glob.glob(os.path.join(coverage_dir, "*_chrM_coverage.txt")):
        with open(cov_path) as f:
            line = f.readline().strip()
        if not line:
            continue
        sample, lib, reads, breadth1, breadth5, breadth10, mean_depth = line.split("\t")
        rows.append({
            "sample": sample,
            "chrM_reads": int(reads),
            "chrM_breadth_1x": float(breadth1),
            "chrM_breadth_5x": float(breadth5),
            "chrM_breadth_10x": float(breadth10),
            "chrM_mean_depth": float(mean_depth),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coverage-dir", required=True)
    p.add_argument("--pe-dir", required=True)
    p.add_argument("--se-dir", required=True)
    p.add_argument("--out", default="qc_summary.tsv")
    p.add_argument("--breadth-threshold", type=float, default=0.90,
                    help="Minimum breadth at 5x to flag as PASS (default 0.90)")
    args = p.parse_args()

    star_rows = collect_star_stats(args.pe_dir, "PE") + collect_star_stats(args.se_dir, "SE")
    cov_rows = collect_coverage_stats(args.coverage_dir)

    if not star_rows:
        print("WARNING: no STAR Log.final.out files found — check --pe-dir/--se-dir paths", file=sys.stderr)
    if not cov_rows:
        print("WARNING: no coverage summary files found — check --coverage-dir path", file=sys.stderr)

    star_df = pd.DataFrame(star_rows)
    cov_df = pd.DataFrame(cov_rows)

    merged = star_df.merge(cov_df, on="sample", how="outer", indicator=True)

    missing_cov = merged[merged["_merge"] == "left_only"]["sample"].tolist()
    missing_star = merged[merged["_merge"] == "right_only"]["sample"].tolist()
    if missing_cov:
        print(f"NOTE: {len(missing_cov)} samples have STAR stats but no coverage summary yet: "
              f"{missing_cov[:5]}{'...' if len(missing_cov) > 5 else ''}", file=sys.stderr)
    if missing_star:
        print(f"NOTE: {len(missing_star)} samples have coverage summary but no STAR stats: "
              f"{missing_star[:5]}{'...' if len(missing_star) > 5 else ''}", file=sys.stderr)

    merged = merged.drop(columns=["_merge"])
    merged["qc_pass_mtdna"] = merged["chrM_breadth_5x"] >= args.breadth_threshold

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} samples to {args.out}")

    # ---- Summary printed to console ----
    print("\n=== Uniquely mapped % by library type ===")
    print(merged.groupby("library_type")["uniquely_mapped_pct"].describe())

    print("\n=== mtDNA breadth at 5x by library type ===")
    print(merged.groupby("library_type")["chrM_breadth_5x"].describe())

    n_fail = (~merged["qc_pass_mtdna"]).sum()
    print(f"\n{n_fail} / {len(merged)} samples fall below "
          f"{args.breadth_threshold:.0%} breadth-at-5x threshold on chrM.")
    if n_fail > 0:
        print("Samples failing mtDNA coverage QC:")
        print(merged.loc[~merged["qc_pass_mtdna"], ["sample", "library_type", "chrM_breadth_5x"]]
              .sort_values("chrM_breadth_5x").to_string(index=False))


if __name__ == "__main__":
    main()
