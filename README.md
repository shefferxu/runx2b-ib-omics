# runx2b-ib-omics
Processed transcriptomic and proteomic analysis scripts for Runx2b-associated intermuscular bone development in Culter alburnus.
# Integrated transcriptomic and proteomic analysis scripts

This repository contains custom Python scripts used to organize processed transcriptomic/proteomic results and generate manuscript figures for the Runx2b-associated intermuscular bone study in Culter alburnus.

## Manuscript figure mapping

The final manuscript uses five main figures. Some internal script names reflect earlier working figure numbers.

| Final manuscript figure | Content | Main scripts |
|---|---|---|
| Fig. 2 | DEG/DEP overview and volcano plots | `scripts/Fig2_DEG_DEP_volcano_plots_make_overview.py` |
| Fig. 3 | mRNA-protein correlation, DEG-DEP overlap, and nine-quadrant analysis | `scripts/Fig3_mRNA_protein_correlation_and_overlap.py`; `scripts/Fig3_nine_quadrant_analysis.py` |
| Fig. 4 | GO/KEGG and biological theme enrichment | `scripts/Fig4_GO_KEGG_theme_enrichment.py`; `scripts/Fig4_theme_enrichment_reanalysis.py` |
| Fig. 5 | Candidate expression heatmaps and candidate network | `scripts/Fig5_candidate_expression_heatmaps.py`; `scripts/Fig5_candidate_mechanism_network.py` |

Fig. 1 contains experimental workflow, genotyping, and alizarin red staining images and was assembled from experimental image materials rather than generated solely by these scripts.

## Data

Processed figure source data and supplementary matrices are deposited in Figshare: [DOI to be added]. Raw RNA-seq reads and raw proteomics files should be deposited separately in NCBI SRA and ProteomeXchange/PRIDE, respectively.

## Environment

Python 3.11 or later is recommended. Install dependencies with:

```bash
pip install -r requirements.txt
```

## Notes

The scripts expect the same processed result directory layout used in the local analysis project. Large raw sequencing files and raw mass spectrometry files are not included in this code repository.
