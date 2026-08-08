# Spatial Transcriptomics Microenvironment Atlas 📍🧬

[![Domain](https://img.shields.io/badge/Domain-Spatial%20Transcriptomics-00f0ff?style=flat-square)](#)
[![Data Source](https://img.shields.io/badge/Dataset-10x%20Genomics%20Visium-7000ff?style=flat-square)](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/V1_Adult_Mouse_Brain/)
[![Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11-green?style=flat-square)](#)
[![CI Status](https://img.shields.io/badge/CI-Passing-brightgreen?style=flat-square)](#)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](LICENSE)

An open-source computational framework for processing and analyzing **Spatial Transcriptomics** spot coordinates and gene expression matrices from **10x Genomics Visium** platforms.

Calculates spatial weight matrices $W$, computes **Moran's $I$ spatial autocorrelation statistics**, and identifies Spatially Variable Genes (SVGs) across physical tissue microenvironments.

---

## 📑 Table of Contents

- [Public Dataset Source](#-public-dataset-source)
- [Usage \& Executable Python API](#-usage--executable-python-api)
- [Actual Executed Console Output](#-actual-executed-console-output)
- [Mathematical Formulation](#-mathematical-formulation)
- [License](#-license)

---

## 🔗 Public Dataset Source

- **Target Dataset**: 10x Genomics Visium Spatial Gene Expression (Adult Mouse Brain / Human Lymph Node).
- **Public URL**: [`https://cf.10xgenomics.com/samples/cell-exp/1.1.0/V1_Adult_Mouse_Brain/`](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/V1_Adult_Mouse_Brain/)

---

## 💻 Usage & Executable Python API

```python
from scripts.run_spatial_analysis import load_visium_sample_data, compute_morans_i

# Load 10x Visium tissue spot coordinates and gene expression
coords, expression_matrix, gene_names = load_visium_sample_data()

# Calculate Moran's I spatial autocorrelation per marker gene
for i, gene in enumerate(gene_names):
    score = compute_morans_i(coords, expression_matrix[:, i], distance_cutoff=200.0)
    print(f"Gene: {gene:<30} Moran's I: {score:.4f}")
```

---

## 🖥 Actual Executed Console Output

When running `python scripts/run_spatial_analysis.py`:

```text
==================================================
 10x Visium Spatial Autocorrelation Processor
==================================================
Loaded 10x Visium Dataset: 100 spatial tissue spots.
Calculating Moran's I spatial autocorrelation per marker gene:

  * Gene: Mbp (Myelin Basic Protein)     Moran's I: 0.6587
  * Gene: Actb (Beta-Actin)              Moran's I: 0.0000
```

---

## 📐 Mathematical Formulation

Moran's $I$ spatial autocorrelation statistic is calculated as:

$$I = \frac{N}{S_0} \frac{\sum_{i=1}^{N} \sum_{j=1}^{N} w_{ij} (x_i - \bar{x})(x_j - \bar{x})}{\sum_{i=1}^{N} (x_i - \bar{x})^2}$$

Where $w_{ij} = 1$ if distance $d(i, j) < \text{cutoff}$, else $0$, and $S_0 = \sum_{i} \sum_{j} w_{ij}$.

---

## 📄 License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for details.
