#!/usr/bin/env python3
"""
Split mito_nuclear_genes.bed into N roughly-equal chunks (by gene count) for
parallel GenotypeGVCFs, reusing the existing GenomicsDB workspace.

Usage:
    python 24b_split_bed_for_genotyping.py --bed mito_nuclear_genes.bed \
        --n-chunks 20 --outdir bed_chunks
"""

import argparse
import os

p = argparse.ArgumentParser()
p.add_argument("--bed", default="mito_nuclear_genes.bed")
p.add_argument("--n-chunks", type=int, default=20)
p.add_argument("--outdir", default="bed_chunks")
args = p.parse_args()

with open(args.bed) as f:
    lines = [l for l in f if l.strip()]

os.makedirs(args.outdir, exist_ok=True)
n = args.n_chunks
chunk_size = (len(lines) + n - 1) // n

for i in range(n):
    chunk_lines = lines[i * chunk_size:(i + 1) * chunk_size]
    if not chunk_lines:
        continue
    with open(os.path.join(args.outdir, f"chunk_{i+1:02d}.bed"), "w") as out:
        out.writelines(chunk_lines)

n_written = min(n, (len(lines) + chunk_size - 1) // chunk_size)
print(f"Split {len(lines)} gene intervals into {n_written} chunks in {args.outdir}/")
