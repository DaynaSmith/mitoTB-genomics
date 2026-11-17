#!/usr/bin/env python3
"""
Build ok_samples.txt from gatk_qc_summary.tsv — one sample_id per line,
status=OK only. This feeds both the GVCF rerun array and the joint
genotyping sample map.

Usage:
    python 22_build_ok_samples_list.py --qc gatk_qc_summary.tsv --out ok_samples.txt
"""

import argparse
import csv

p = argparse.ArgumentParser()
p.add_argument("--qc", default="gatk_qc_summary.tsv")
p.add_argument("--out", default="ok_samples.txt")
args = p.parse_args()

samples = []
with open(args.qc) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["status"] == "OK":
            samples.append(row["sample"])

with open(args.out, "w") as out:
    for s in samples:
        out.write(s + "\n")

print(f"Wrote {len(samples)} OK samples to {args.out}")
