#!/usr/bin/env python3
"""
Build the final mtDNA gene expression matrix from personalized remap output.

For each sample: uses the MAIN pass (mismatch<=8) counts for protein-coding
and rRNA genes, and the tRNA pass (mismatch=0) counts specifically for the
22 tRNA genes — matching the paper's two-pass approach.

Uses column 4 (reverse-stranded) of STAR's ReadsPerGene.out.tab by default —
confirmed as the correct strand assignment from manual inspection (see
conversation notes). Override with --strand-column if your library differs.

Usage:
    python aggregate_mtdna_counts.py \
        --manifest bam_manifest.txt \
        --gtf chrM_genes.gtf \
        --remap-dir personalized_remap_out \
        --out mtdna_expression_matrix.tsv
"""

import argparse
import csv
import os
import re

STAR_SUMMARY_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}


def parse_gtf_genes(gtf_path):
    """Return {gene_id: (gene_name, gene_type)} from gene-level GTF lines."""
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    gene_type_re = re.compile(r'gene_type "([^"]+)"')
    genes = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gid_m = gene_id_re.search(fields[8])
            gname_m = gene_name_re.search(fields[8])
            gtype_m = gene_type_re.search(fields[8])
            if gid_m:
                gid = gid_m.group(1)
                gname = gname_m.group(1) if gname_m else gid
                gtype = gtype_m.group(1) if gtype_m else "unknown"
                genes[gid] = (gname, gtype)
    return genes


def read_star_counts(tab_path, strand_col):
    """Return {gene_id: count} from a STAR ReadsPerGene.out.tab file."""
    counts = {}
    with open(tab_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[0]
            if gene_id in STAR_SUMMARY_ROWS:
                continue
            counts[gene_id] = int(fields[strand_col])
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--remap-dir", default="personalized_remap_out")
    p.add_argument("--out", default="mtdna_expression_matrix.tsv")
    p.add_argument("--strand-column", type=int, default=3,
                    help="0-indexed column in ReadsPerGene.out.tab (default 3 = "
                         "4th column = reverse-stranded, confirmed correct for this data)")
    args = p.parse_args()

    genes = parse_gtf_genes(args.gtf)
    print(f"Loaded {len(genes)} genes from {args.gtf}")
    trna_genes = {gid for gid, (name, gtype) in genes.items() if gtype == "Mt_tRNA"}
    print(f"  of which {len(trna_genes)} are Mt_tRNA (will use tRNA-pass counts for these)")

    samples = []
    with open(args.manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            samples.append(parts[0])

    matrix = {}  # gene_id -> {sample_id: count}
    missing_main = []
    missing_trna = []

    for sample_id in samples:
        main_path = os.path.join(args.remap_dir, sample_id, "main_pass",
                                  f"{sample_id}_main_ReadsPerGene.out.tab")
        trna_path = os.path.join(args.remap_dir, sample_id, "trna_pass",
                                  f"{sample_id}_trna_ReadsPerGene.out.tab")

        main_counts = {}
        trna_counts = {}

        if os.path.exists(main_path):
            main_counts = read_star_counts(main_path, args.strand_column)
        else:
            missing_main.append(sample_id)

        if os.path.exists(trna_path):
            trna_counts = read_star_counts(trna_path, args.strand_column)
        else:
            missing_trna.append(sample_id)

        for gene_id in genes:
            source = trna_counts if gene_id in trna_genes else main_counts
            count = source.get(gene_id)  # None if missing
            matrix.setdefault(gene_id, {})[sample_id] = count

    # Write wide matrix: rows = genes, columns = samples
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        header = ["gene_id", "gene_name", "gene_type"] + samples
        writer.writerow(header)
        for gene_id in sorted(genes):
            gene_name, gene_type = genes[gene_id]
            row = [gene_id, gene_name, gene_type]
            row += [matrix[gene_id][s] if matrix[gene_id][s] is not None else "NA" for s in samples]
            writer.writerow(row)

    print(f"\nWrote matrix ({len(genes)} genes x {len(samples)} samples) to {args.out}")

    if missing_main:
        print(f"\n{len(missing_main)} samples missing MAIN pass output:")
        print("  " + ", ".join(missing_main[:10]) + (" ..." if len(missing_main) > 10 else ""))
    if missing_trna:
        print(f"\n{len(missing_trna)} samples missing tRNA pass output:")
        print("  " + ", ".join(missing_trna[:10]) + (" ..." if len(missing_trna) > 10 else ""))
    if not missing_main and not missing_trna:
        print("\nAll samples had both MAIN and tRNA pass output — no gaps.")


if __name__ == "__main__":
    main()
