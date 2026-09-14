# 01_reference_panel: mtDNA Imputation Reference Panel Construction

**Objective:** To develop an mtDNA imputation reference panel using publicly available mitochondrial genome sequences.

The panel was constructed following the workflow outlined by McInerney et al. (2021), implemented here using the reference alignment and resources distributed with [MitoImpute](https://github.com/sjfandrews/MitoImpute) (Sherman et al.).

## Pipeline summary

```
sequence.fasta (68,670 raw sequences, GenBank/MitoMap, 18 May 2026)
        │
        ▼
[01_align.pbs]
  1a. Split into batches of 2,500 sequences
  1b. MAFFT v7.453 --add --keeplength --reorder --thread 16
      against hsapiensCRS7k.fasta (7,747-seq reference MSA, MitoImpute repo)
  1c. Merge aligned batches (reference seqs taken once; de-duplicated)
  1d. Remove alignment columns where rCRS has a gap (insertions relative to rCRS)
        │
        ▼  aligned/full_msa_rcrs_trimmed.fasta
[02_filter.pbs]
  Exclude sequences with ≥5 ambiguous bases OR ≥8 gap characters
  (rCRS always retained regardless of stats)
        │
        ▼  filtered/filtered_msa.fasta
[03_msa_to_vcf.pbs → msa_to_vcf.py]
  Memory-mapped column-wise MSA → VCF conversion
  rCRS defines REF; gaps/IUPAC codes → missing; invariant sites skipped
        │
        ▼  data/vcf/mtdna_raw.vcf
  bcftools: remove any residual invariant/fully-missing sites,
  bcftools norm -m -any (multiallelic → biallelic)
        │
        ▼
  bcftools --haploid2diploid → IMPUTE2 format (.hap.gz/.legend.gz/.samples)
  + zero-recombination map
        │
        ▼  impute2_panel/
HaploGrep v3 / PhyloTree v17.0 haplogroup assignment
        │
        ▼  haplogroups/
```

## Directory contents

```
01_reference_panel/
├── README.md
├── impute2_panel/                   # Final QC'd .hap.gz, .legend.gz + recombination map
├── haplogroups/
│   ├── haplogrep_output.txt         # HaploGrep v3 / PhyloTree v17.0 calls for panel samples
│   └── Reference_panel_haplogroup_composition_donut.png
└── scripts/
   ├── 01_align.pbs
   ├── 02_filter.pbs
   ├── 03_msa_to_vcf.pbs
   └── msa_to_vcf.py
```

## Methods

### Sequence identification and download

Complete human mtDNA sequences were identified in GenBank using the MitoMap search function with the following search terms:

```
(015400[SLEN]:016700[SLEN]) AND Homo[Organism] AND mitochondrion[FILT]
NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis
OR heidelbergensis OR consensus)
```

A total of 68,670 mtDNA sequences were identified and downloaded in FASTA format (`sequence.fasta`) on 18 May 2026.

### Alignment (`01_align.pbs`)

Raw sequences were split into batches of 2,500. Each batch was aligned against a 7,747-sequence reference MSA (`hsapiensCRS7k.fasta`, sourced from the [MitoImpute](https://github.com/sjfandrews/MitoImpute) repository; McInerney et al., 2021) using MAFFT v7.453 (Katoh & Frith, 2012), with the `--add --keeplength --reorder` options and 16 threads. The rCRS (NC_012920.1) was present in the reference MSA and used to anchor standard mitochondrial genome coordinates.

Aligned batches were merged: reference sequences were retained once (from the first batch output), and only newly aligned, non-duplicate sequences were appended from each subsequent batch (`full_msa.fasta`).

Alignment columns corresponding to insertions relative to the rCRS (i.e., columns where the rCRS carries a gap) were then removed, producing `full_msa_rcrs_trimmed.fasta`.

### Quality filtering (`02_filter.pbs`)

Sequences with **≥5 ambiguous (non-ACGT) bases or ≥8 gap characters** were excluded; the rCRS was retained regardless of its stats. This produced a filtered MSA of 58,808 sequences spanning 16,569 mitochondrial genome positions. The 8-gap threshold follows McInerney et al. (2021) and avoids excluding haplogroup B sequences, which carry a characteristic 9 bp deletion (Ramírez-Molina et al., 2026).

> **Note for methods write-up:** the script excludes sequences at *≥5* ambiguous bases (i.e. 5 or more), not *>5* ("more than five") — worth double-checking this wording against your thesis text so the threshold described matches the script exactly.

A per-sequence report (`filter_report.tsv`) and list of excluded accessions (`excluded_ids.txt`) are written alongside the filtered MSA.

### MSA → VCF conversion (`03_msa_to_vcf.pbs`, `msa_to_vcf.py`)

The filtered MSA was converted to VCF format using a custom, memory-mapped Python script (`msa_to_vcf.py`) rather than loading the full sequence matrix into memory, keeping RAM usage roughly constant (~60 KB per column) regardless of panel size. At each alignment column, the rCRS base defined the REF allele; gap characters and IUPAC ambiguous bases were coded as missing (`.`).

Where the rCRS itself carried an "N" (position 3,107), no VCF record was written for that site, but the rCRS coordinate counter was still incremented — preserving correct rCRS numbering at all downstream positions (Salas et al., 2012). Invariant sites (no alternate allele among samples) were also identified and skipped at this stage.

### Downstream QC and IMPUTE2 formatting

The remaining steps — removing any residual invariant/fully-missing sites and normalising multiallelic sites with `bcftools norm -m -any`, then converting to IMPUTE2 format (`.hap.gz`, `.legend.gz`, `.samples`) via `bcftools --haploid2diploid`, and generating a zero-recombination map — follow as previously described.

### Haplogroup assignment

Haplogroups for individual reference panel samples were obtained using HaploGrep v3 and PhyloTree v17.0. The resulting composition (see `haplogroups/Reference_panel_haplogroup_composition_donut.png`) is dominated by haplogroups H, L, and M, with broad representation across most major mtDNA haplogroups (A–Z), consistent with a globally diverse panel suited to imputation across multiple mitochondrial lineages.

## Software

- MAFFT v7.453 (Katoh & Frith, 2012)
- Biopython 1.81
- BCFtools v1.3 (Danecek et al., 2021)
- IMPUTE2 (Howie et al., 2009)
- HaploGrep v3
- PhyloTree v17.0

## References

- McInerney, T. et al. (2021). A globally diverse reference alignment and panel for imputation of mitochondrial DNA variants. *BMC Bioinformatics*.
- Katoh, K. & Frith, M.C. (2012). Adding unaligned sequences into an existing alignment using MAFFT and LAST. *Bioinformatics*.
- Salas, A. et al. (2012). A cautionary note on switching mitochondrial DNA reference sequences in forensic genetics. *Forensic Science International: Genetics*.
- Ramírez-Molina, D. et al. (2026). Assay for numerical and graphical detection of Native American mitochondrial haplogroup B with real-time PCR. *Molecular Biology Reports*.
- Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*.
- Howie, B.N. et al. (2009). A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. *PLoS Genetics*.
