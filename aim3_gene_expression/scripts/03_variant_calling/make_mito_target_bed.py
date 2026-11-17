#!/usr/bin/env python3
"""
Build a BED file of gene coordinates for the MitoCarta3.0 nuclear-encoded
mitochondrial gene set, using your existing GENCODE GTF for coordinates.

MitoCarta3.0 download: https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
(download the "Human MitoCarta3.0" Excel file, then export/save the gene
symbol column as a plain text file, one gene symbol per line — see notes
below for a quick way to do this without opening Excel on the HPC.)

Usage:
    python 17_make_mito_target_bed.py \
        --gene-list mitocarta3_human_symbols.txt \
        --gtf gencode.v50.annotation.gtf \
        --out mito_nuclear_genes.bed
"""

import argparse
import re


def load_gene_list(path):
    genes = set()
    with open(path) as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    return genes


def parse_gtf_genes(gtf_path, target_genes):
    """Extract chrom/start/end for gene-level features matching target_genes."""
    found = {}
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            m = gene_name_re.search(fields[8])
            if not m:
                continue
            gene_name = m.group(1)
            if gene_name in target_genes:
                chrom, start, end = fields[0], int(fields[3]) - 1, int(fields[4])  # BED is 0-based
                found[gene_name] = (chrom, start, end)
    return found


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gene-list", required=True, help="Plain text file, one gene symbol per line")
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", default="mito_nuclear_genes.bed")
    p.add_argument("--flank", type=int, default=0,
                    help="Optional bp of flanking sequence to add around each gene (default 0)")
    args = p.parse_args()

    target_genes = load_gene_list(args.gene_list)
    print(f"Loaded {len(target_genes)} target gene symbols")

    found = parse_gtf_genes(args.gtf, target_genes)
    print(f"Matched {len(found)} / {len(target_genes)} genes to GTF coordinates")

    missing = target_genes - set(found.keys())
    if missing:
        print(f"WARNING: {len(missing)} genes not found in GTF (possibly renamed/deprecated symbols):")
        print("  " + ", ".join(sorted(missing)[:20]) + (" ..." if len(missing) > 20 else ""))
        print("  Consider checking these against a symbol-alias table (e.g. HGNC) if the count is large.")

    with open(args.out, "w") as out:
        for gene_name, (chrom, start, end) in sorted(found.items()):
            start = max(0, start - args.flank)
            end = end + args.flank
            out.write(f"{chrom}\t{start}\t{end}\t{gene_name}\n")

    print(f"Wrote {len(found)} intervals to {args.out}")


if __name__ == "__main__":
    main()
