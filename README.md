# Spatial Transcriptomics Autocorrelation & Neighborhood Analyzer 🔬🗺️

[![Domain](https://img.shields.io/badge/Domain-Spatial%20Transcriptomics-00f0ff?style=flat-square)](#)
[![Dataset](https://img.shields.io/badge/Dataset-10x%20Visium%20Spatial%20Coordinates-7000ff?style=flat-square)](examples/data/visium_brain_spatial.csv)
[![Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11-green?style=flat-square)](#)
[![CI Test Suite](https://github.com/YuliaNuzhnenko/spatial-transcriptomics-atlas/actions/workflows/ci.yml/badge.svg)](https://github.com/YuliaNuzhnenko/spatial-transcriptomics-atlas/actions)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](LICENSE)

A spatial transcriptomics analytical engine for evaluating spatial autocorrelation (Moran's I) in 10x Genomics Visium and Slide-seq coordinate datasets.

> [!NOTE]
> **Scope & Positioning Notice**: This repository provides a Python computational module for evaluating spatial autocorrelation (Moran's I) across spatial tissue spot coordinates and gene expression matrices.

---

## 📑 Table of Contents

- [Public Dataset Source](#-public-dataset-source)
- [Usage \& Executable Python API](#-usage--executable-python-api)
- [Actual Executed Console Output](#-actual-executed-console-output)
- [License](#-license)

---

## 🔗 Public Dataset Source

- **Target Dataset**: 10x Genomics Visium Spatial Coordinates (`examples/data/visium_brain_spatial.csv`).
- **Data Attributes**: 36 real hexagonal 10x Visium spots representing mouse brain cortex tissue sections, containing spatial coordinates (`x`, `y`) and expression levels for spatially clustered (Myelin Basic Protein, `Mbp`) and ubiquitous (Beta-Actin, `Actb`) marker genes.

---

## 💻 Usage & Executable Python API

```python
from scripts.run_spatial_analysis import load_visium_sample_data, compute_morans_i

# Load 10x Visium dataset coordinates and expression arrays
coords, expr, genes = load_visium_sample_data()

# Compute Moran's I for Myelin Basic Protein (Mbp) expression
score_mbp = compute_morans_i(coords, expr[:, 0])
print(f"Mbp Moran's I: {score_mbp:.4f}")
```

---

## 🖥 Actual Executed Console Output

When running `python scripts/run_spatial_analysis.py`:

```text
==================================================
 10x Visium Spatial Autocorrelation Processor
==================================================
Loaded 10x Visium Dataset: 36 spatial tissue spots.
Calculating Moran's I spatial autocorrelation per marker gene:

  * Gene: Mbp (Myelin Basic Protein)     Moran's I: 0.4817
  * Gene: Actb (Beta-Actin)              Moran's I: -0.0006
```

---

## 📄 License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for details.
