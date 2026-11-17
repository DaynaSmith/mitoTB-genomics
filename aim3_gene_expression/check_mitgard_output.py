#!/usr/bin/env python3
"""
Sweep mitgard_out/ and summarize every sample's mitogenome assembly result.
Flags: missing output, unusually short/fragmented assemblies, multi-sequence
files (should be exactly 1 sequence per sample).

Usage:
    python 15_check_mitgard_output.py --manifest fastq_manifest.txt \
        --outdir mitgard_out --out mitgard_qc_summary.tsv
"""

import argparse
import glob
import os

REF_LENGTH = 16569  # rCRS length


def parse_fasta_lengths(fasta_path):
    """Return list of (header, length) for each sequence in a FASTA."""
    seqs = []
    header = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, length))
                header = line[1:]
                length = 0
            else:
                length += len(line)
    if header is not None:
        seqs.append((header, length))
    return seqs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="fastq_manifest.txt used for the run")
    p.add_argument("--outdir", default="mitgard_out")
    p.add_argument("--out", default="mitgard_qc_summary.tsv")
    p.add_argument("--min-breadth-pct", type=float, default=90.0,
                    help="Flag samples whose assembled length is below this %% of rCRS length")
    args = p.parse_args()

    samples = []
    with open(args.manifest) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            sample_id = parts[0]
            lib_type = parts[3] if len(parts) > 3 else "?"
            samples.append((sample_id, lib_type))

    rows = []
    for sample_id, lib_type in samples:
        fasta_path = os.path.join(args.outdir, sample_id, f"{sample_id}_mitogenome.fa")
        if not os.path.exists(fasta_path):
            rows.append({
                "sample": sample_id, "library_type": lib_type,
                "status": "MISSING", "n_sequences": 0, "length": 0,
                "pct_of_rcrs": 0.0,
            })
            continue

        seqs = parse_fasta_lengths(fasta_path)
        if not seqs:
            rows.append({
                "sample": sample_id, "library_type": lib_type,
                "status": "EMPTY_FASTA", "n_sequences": 0, "length": 0,
                "pct_of_rcrs": 0.0,
            })
            continue

        total_length = sum(length for _, length in seqs)
        pct = 100.0 * total_length / REF_LENGTH
        n_seqs = len(seqs)

        if n_seqs > 1:
            status = "MULTI_SEQUENCE"
        elif pct < args.min_breadth_pct:
            status = "SHORT"
        else:
            status = "OK"

        rows.append({
            "sample": sample_id, "library_type": lib_type,
            "status": status, "n_sequences": n_seqs, "length": total_length,
            "pct_of_rcrs": round(pct, 2),
        })

    # Write full summary
    with open(args.out, "w") as f:
        f.write("sample\tlibrary_type\tstatus\tn_sequences\tlength\tpct_of_rcrs\n")
        for r in rows:
            f.write(f"{r['sample']}\t{r['library_type']}\t{r['status']}\t"
                     f"{r['n_sequences']}\t{r['length']}\t{r['pct_of_rcrs']}\n")

    # Console summary
    n_total = len(rows)
    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_missing = sum(1 for r in rows if r["status"] == "MISSING")
    n_short = sum(1 for r in rows if r["status"] == "SHORT")
    n_multi = sum(1 for r in rows if r["status"] == "MULTI_SEQUENCE")
    n_empty = sum(1 for r in rows if r["status"] == "EMPTY_FASTA")

    print(f"Wrote {n_total} samples to {args.out}")
    print(f"  OK:             {n_ok}")
    print(f"  MISSING:        {n_missing}  (job likely still running or failed)")
    print(f"  SHORT:          {n_short}  (assembled < {args.min_breadth_pct}% of rCRS length)")
    print(f"  MULTI_SEQUENCE: {n_multi}  (more than 1 contig — check assembly)")
    print(f"  EMPTY_FASTA:    {n_empty}")

    if n_missing or n_short or n_multi or n_empty:
        print("\nSamples needing attention:")
        for r in rows:
            if r["status"] != "OK":
                print(f"  {r['sample']} ({r['library_type']}): {r['status']} "
                      f"— {r['length']} bp ({r['pct_of_rcrs']}% of rCRS)")


if __name__ == "__main__":
    main()
