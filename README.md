# *PlasPredict* Bacterial Taxonomic Partitioning Tool

A Flask web application for predicting the host range of plasmids from their DNA sequences. Given a plasmid FASTA file, the tool runs a bioinformatics pipeline to extract biological features and uses a machine learning model to classify the plasmid's likely host range.

The web application is hosted at: https://plaspart.serve.scilifelab.se/

---

## Overview

The model classifies plasmids into one of four host range categories:

| Category | Description |
|---|---|
| ≥ Phylum | Plasmid likely to proliferate across bacterial phyla |
| Order/Class | Plasmid likely restricted to one bacterial class or phylum |
| Genus/Family | Plasmid likely restricted to one bacterial family or order |
| ≤ Species | Plasmid unlikely to move beyond a single bacterial genus |


The prediction is based on a combination of sequence-derived features (k-mer frequencies, plasmid size), molecular typing (incompatibility groups, conjugation systems), antibiotic resistance gene content, and optional isolation source metadata.

---

## Dependencies

### Python packages

- Flask
- Biopython
- scikit-learn / joblib
- NumPy / Pandas

### External tools (must be available on `PATH`)

| Tool | Purpose |
|---|---|
| [Prodigal](https://github.com/hyacc/prodigal) | Protein-coding gene prediction |
| [MacSyFinder](https://macsyfinder.readthedocs.io) + CONJScan | Conjugation system detection |
| [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder) | Incompatibility (Inc) type typing |
| [RGI](https://github.com/arpcard/rgi) (CARD) | Antibiotic resistance gene identification |

---

## Feature extraction pipeline

For each submitted sequence, the pipeline:

1. **k-mer frequencies** — computes canonical 3-mer frequency distributions across the full plasmid sequence.
2. **Plasmid size** — records total sequence length.
3. **Conjugation systems** — calls Prodigal for gene prediction, then MacSyFinder (CONJScan/Plasmids model) to identify conjugation system type.
4. **Inc typing** — runs PlasmidFinder to detect incompatibility group markers.
5. **Resistance genes** — runs RGI against the CARD database to identify antibiotic resistance genes, drug classes, and resistance mechanisms. Overlapping hits are resolved by tier (Perfect > Strict > Loose) and bitscore.
6. **Isolation source** — optional one-hot encoded metadata feature.

All features are assembled into a single feature vector and passed to the pre-trained XGBoost model.

---

## File size limit

The maximum upload size is **50 MB**. Requests exceeding this limit receive a `413` error.

---

