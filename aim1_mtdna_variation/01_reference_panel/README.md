# 01_reference_panel: mtDNA Imputation Reference Panel Construction

**Objective:** To develop an mtDNA imputation reference panel using publicly available mitochondrial genome sequences.

## Methods

The mtDNA imputation reference panel was constructed using the workflow outlined by McInerney et al. (2021).

### Sequence identification and download

Complete human mtDNA sequences were identified in GenBank using the MitoMap search function with the following search terms:

```
(015400[SLEN]:016700[SLEN]) AND Homo[Organism] AND mitochondrion[FILT]
NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis
OR heidelbergensis OR consensus)
```

A total of 68,670 mtDNA sequences were identified and downloaded in FASTA format on 18 May 2026.

### Alignment

The downloaded sequences were aligned in batches of 2,500 to a reference multiple sequence alignment (MSA) comprising 7,747 complete mtDNA sequences (McInerney et al., 2021) using MAFFT v7.453 (Katoh & Frith, 2012). The `--add` and `--keeplength` options were used to add the sequences to the existing alignment without altering its gap structure. The rCRS was included in the reference MSA and used to maintain standard mitochondrial genome coordinates. Alignment columns corresponding to insertions relative to the rCRS were removed. The rCRS placeholder at position m.3107 was retained during alignment processing to preserve coordinate consistency across downstream positions (Salas et al., 2012).

### Filtering

Sequences containing more than five ambiguous bases or eight or more gap characters were excluded, resulting in a filtered MSA comprising 58,808 sequences spanning 16,569 mitochondrial genome positions. The eight-gap exclusion threshold was recommended by McInerney et al. (2021) to avoid excluding sequences belonging to haplogroup B, which is characterised by a nine base-pair deletion (Ramírez-Molina et al., 2026).

### MSA to VCF conversion

The filtered MSA was converted to VCF format using a custom Python script. At each alignment column, the rCRS base defined the reference (REF) allele, while gap characters and IUPAC ambiguous bases were coded as missing. Sites at which the rCRS carried an "N" (i.e., position 3,107) were excluded from the VCF; however, their position was retained in the coordinate system to preserve rCRS numbering at all downstream sites.

Invariant and fully missing sites were removed using BCFtools v1.3 (Danecek et al., 2021). Multi-allelic sites were converted to biallelic records using:

```
bcftools norm -m -any
```

### IMPUTE2 panel formatting

The reference panel VCF file was converted to IMPUTE2-compatible (Howie et al., 2009) formats using BCFtools, producing `.hap.gz`, `.legend.gz`, and `.samples` files. The `--haploid2diploid` flag was used to convert haploid mtDNA genotypes to phased homozygous diploid format to accommodate the input requirements for IMPUTE2.

A recombination map was created by assigning a recombination rate of zero to all variant sites, consistent with the non-recombining nature of the mitochondrial genome.

### Haplogroup assignment

Haplogroups for individual samples in the reference panel were obtained using HaploGrep v3 and PhyloTree v17.0.

## Pipeline summary

```
MitoMap/GenBank search (68,670 sequences, 18 May 2026)
        │
        ▼
MAFFT v7.453 alignment (--add --keeplength) to McInerney et al. (2021)
7,747-sequence reference MSA, batches of 2,500
        │
        ▼
Remove insertion columns relative to rCRS; retain m.3107 placeholder
        │
        ▼
Filter: >5 ambiguous bases OR ≥8 gap characters excluded
        → 58,808 sequences × 16,569 positions
        │
        ▼
Custom Python script: MSA → VCF (rCRS = REF; gaps/IUPAC codes = missing;
m.3107 excluded from VCF but retained in coordinate system)
        │
        ▼
BCFtools v1.3: remove invariant/fully missing sites;
bcftools norm -m -any (multiallelic → biallelic)
        │
        ▼
BCFtools --haploid2diploid → IMPUTE2 format
(.hap.gz, .legend.gz, .samples) + zero-recombination map
        │
        ▼
HaploGrep v3 / PhyloTree v17.0 haplogroup assignment
```

## Directory contents

```
01_reference_panel/
├── README.md
├── data/
│   ├── raw_sequences/          # 68,670 downloaded FASTA sequences (18 May 2026)
│   ├── reference_msa/          # McInerney et al. (2021) 7,747-sequence MSA + rCRS
│   └── filtered_msa/           # Filtered MSA: 58,808 sequences × 16,569 positions
├── scripts/
│   ├── 01_align_mafft.pbs            # MAFFT --add --keeplength, batches of 2,500
│   ├── 02_filter_sequences.pbs       # Ambiguous base / gap character filtering
│   └── 03_msa_to_vcf.py              # Custom MSA → VCF conversion script
├── vcf/                         # Reference panel VCF (post bcftools norm -m -any)
├── impute2_panel/               # .hap.gz, .legend.gz, .samples, recombination map
└── haplogroups/                 # HaploGrep3 / PhyloTree v17.0 haplogroup calls
```

## Inputs

- MitoMap/GenBank search results (FASTA sequences)
- McInerney et al. (2021) reference MSA (7,747 sequences) including rCRS

## Outputs

- Filtered MSA (58,808 sequences × 16,569 positions)
- Reference panel VCF (biallelic, invariant/missing sites removed)
- IMPUTE2-format reference panel (`.hap.gz`, `.legend.gz`, `.samples`) and recombination map
- Haplogroup assignments for reference panel samples

## Software

- MAFFT v7.453 (Katoh & Frith, 2012)
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
