# GENIE3 Gene Regulatory Network Inference: Complete Workflow

Welcome! This guide walks you through my process for inferring gene regulatory networks from expression data using the [GENIE3](https://academic.oup.com/bioinformatics/article/26/18/2361/193348) algorithm. The workflow covers **data preparation**, **running GENIE3**, and **extracting ranked regulatory links**. All steps are explained with code and context, so you can adapt them for your own projects.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Data Preparation](#data-preparation)
    - [1. Load and Filter Expression Data](#1-load-and-filter-expression-data)
    - [2. Filter Genes of Interest](#2-filter-genes-of-interest)
    - [3. Transpose and Save Expression Matrix](#3-transpose-and-save-expression-matrix)
3. [Running GENIE3](#running-genie3)
    - [1. Prepare Input Files](#1-prepare-input-files)
    - [2. Run GENIE3 Analysis](#2-run-genie3-analysis)
    - [3. Extract and Save Regulatory Links](#3-extract-and-save-regulatory-links)
4. [Tips & Troubleshooting](#tips--troubleshooting)
5. [References](#references)

---

## Introduction

**GENIE3** is a powerful machine learning-based method for inferring gene regulatory networks from expression data. It uses ensembles of regression trees (Random Forests or Extra Trees) to predict the importance of each gene as a regulator of others.

This workflow is designed for **tabular gene expression data** (e.g., RNA-seq counts or normalized values), and can be adapted to your own datasets.

---

## Data Preparation

### 1. Load and Filter Expression Data

Start with your raw count or expression matrix (e.g., `Count_data.txt`). The file should have genes as rows and samples as columns.

```python
import pandas as pd

# Load the count data (tab-separated, first column is gene ID)
df = pd.read_csv('Count_data.txt', sep='\t', index_col=0)
```

### 2. Filter Genes of Interest

Suppose you have a list of gene IDs you want to keep (e.g., from `All_genes_detected.txt`):

```python
# Load gene list (one gene ID per line)
with open('All_genes_detected.txt') as f:
    gene_list = [line.strip() for line in f]

# Filter the count table to keep only genes in gene_list
df_filtered = df.loc[df.index.intersection(gene_list)]
```

### 3. Transpose and Save Expression Matrix

GENIE3 expects **samples as rows** and **genes as columns**. Transpose the table and save as CSV:

```python
# Remove the first row if it contains sample IDs (optional)
df_filtered_no_sampleid = df_filtered.iloc[1:, :]

# Transpose so that samples are rows, genes are columns
df_transposed = df_filtered_no_sampleid.transpose()

# Save as a comma-separated file for GENIE3
df_transposed.to_csv('expression_matrix.csv', index=False)
```

---

## Running GENIE3

### 1. Prepare Input Files

- **expression_matrix.csv**: Genes as columns, samples as rows, header is gene names.
- **regulators.txt**: List of candidate regulator gene names (one per line, no header).

### 2. Run GENIE3 Analysis

Below is a streamlined version of my GENIE3 workflow. You can use this as a script (e.g., `run_genie3.py`):

```python
import pandas as pd
from numpy import *
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.tree import BaseDecisionTree
import time

# --- Load expression matrix ---
df = pd.read_csv('expression_matrix.csv')
data = df.to_numpy()

# --- Extract gene names from header ---
with open('expression_matrix.csv') as f:
    gene_names = f.readline().strip().split(',')

# --- Load regulators list ---
with open('regulators.txt') as f:
    regulators = [line.strip() for line in f]

# --- Define GENIE3 functions (see full code in repo) ---
# ... (Paste your GENIE3, get_link_list, etc. functions here) ...

# --- Run GENIE3 ---
VIM = GENIE3(data, gene_names=gene_names, regulators=regulators, tree_method='ET', nthreads=12)

# --- Extract ranked regulatory links ---
get_link_list(VIM, gene_names=gene_names, regulators=regulators, file_name='results.csv')
```

**Key parameters:**
- `tree_method='ET'`: Use Extra Trees (faster, often similar results to Random Forests).
- `nthreads=12`: Use 12 CPU threads for speed (adjust for your machine).
- `regulators=regulators`: Only consider edges from candidate regulators.

### 3. Extract and Save Regulatory Links

The `get_link_list` function writes a ranked list of predicted regulatory links to `results.csv`:

| Regulator | Target | Score     |
|-----------|--------|-----------|
| TF1       | GeneA  | 0.123456  |
| TF2       | GeneB  | 0.098765  |
| ...       | ...    | ...       |

---

## Tips & Troubleshooting

- **Input format:** Make sure your expression matrix header matches your gene names and regulator list.
- **Zero scores:** Links with zero importance are randomly shuffled in the output.
- **Performance:** For large datasets, increase `nthreads` to speed up computation.
- **Dependencies:** Requires `numpy`, `pandas`, `scikit-learn`.

---

## References

- [GENIE3 original paper](https://academic.oup.com/bioinformatics/article/26/18/2361/193348)
- [GENIE3 Python implementation](https://github.com/vahuynh/GENIE3)
- [scikit-learn documentation](https://scikit-learn.org/stable/)

---

## Acknowledgements

This workflow is adapted from the GENIE3 Python implementation and tailored for reproducible, large-scale gene network inference. Feedback and contributions are welcome!

---

**Happy network inference!**  
If you have questions or suggestions, feel free to open an issue or pull request.

---

Let me know if you want a more technical README, more details on the code, or a more visual/diagram-based guide!
