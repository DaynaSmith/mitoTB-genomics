#!/usr/bin/env python3
"""
Aggregate nuclear (genome-wide) gene expression from the ORIGINAL STAR runs
(aligned_PE/aligned_SE), which already have --quantMode GeneCounts output --
no realignment needed, just gathering what already exists.

Needed for Analysis 1's nuclear co-expression + GO enrichment step (Phase F).

Uses column 4 (reverse-stranded), matching what was confirmed correct for
the mtDNA personalized remap earlier in this project.

Usage:
    python 37_aggregate_nuclear_counts.py \
        --manifest bam_manifest.txt \
        --pe-dir aligned_PE --se-dir aligned_SE \
        --out nuclear_expression_matrix.tsv
"""

import argparse
import os

STAR_SUMMARY_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--pe-dir", default="aligned_PE")
    p.add_argument("--se-dir", default="aligned_SE")
    p.add_argument("--out", default="nuclear_expression_matrix.tsv")
    p.add_argument("--strand-column", type=int, default=3)
    args = p.parse_args()

    samples = []
    with open(args.manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                samples.append((parts[0], parts[2]))  # sample_id, library_type

    all_genes = None
    matrix = {}
    missing = []

    for sample_id, lib_type in samples:
        base_dir = args.pe_dir if lib_type == "PE" else args.se_dir
        tab_path = os.path.join(base_dir, sample_id, f"{sample_id}_ReadsPerGene.out.tab")
        if not os.path.exists(tab_path):
            missing.append(sample_id)
            continue

        with open(tab_path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                gene_id = fields[0]
                if gene_id in STAR_SUMMARY_ROWS:
                    continue
                count = int(fields[args.strand_column])
                matrix.setdefault(gene_id, {})[sample_id] = count

        if all_genes is None:
            with open(tab_path) as f:
                all_genes = [line.split("\t")[0] for line in f
                             if line.split("\t")[0] not in STAR_SUMMARY_ROWS]

    sample_ids = [s for s, _ in samples if s not in missing]

    with open(args.out, "w") as out:
        out.write("gene_id\t" + "\t".join(sample_ids) + "\n")
        for gene_id in all_genes:
            row = [gene_id] + [str(matrix.get(gene_id, {}).get(s, "NA")) for s in sample_ids]
            out.write("\t".join(row) + "\n")

    print(f"Wrote {len(all_genes)} genes x {len(sample_ids)} samples to {args.out}")
    if missing:
        print(f"\n{len(missing)} samples missing ReadsPerGene.out.tab:")
        print("  " + ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))


if __name__ == "__main__":
    main()
