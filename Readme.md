
# PMD Electron Transfer Rate Calculator

## Overview

The **PMD Electron Transfer Rate Calculator** is a Python workflow for predicting biological electron transfer (ET) rates directly from protein structures.

The pipeline combines structural analysis of Protein Data Bank (PDB) files with the **Page–Moser–Dutton (PMD) empirical electron transfer model** to estimate electron transfer rate constants between redox-active amino acid residues.

Unlike centroid-based approaches, this implementation measures the **minimum edge-to-edge atomic distance** between donor and acceptor chromophores, providing a more realistic estimate of electron tunneling pathways.

The workflow automatically:

1. Reads protein structures
2. Identifies chromophoric amino acids
3. Computes all donor–acceptor edge distances
4. Builds pairwise distance matrices
5. Calculates Page–Moser–Dutton (PMD) electron transfer rates
6. Produces summary statistics
7. Identifies the fastest ET pathways
8. Generates heatmaps for visualization

---

# Features

- ✔ Automatic PDB parsing
- ✔ Identification of electron transfer chromophores

  - Tryptophan (TRP)
  - Tyrosine (TYR)
  - Histidine (HIS)
  - Methionine (MET)
  - Cysteine (CYS)
- ✔ Edge-to-edge atomic distance calculations
- ✔ Complete donor → acceptor pair analysis
- ✔ Page–Moser–Dutton electron transfer calculations
- ✔ CSV outputs
- ✔ Heatmap visualization
- ✔ Ranking of fastest electron transfer pathways
- ✔ Summary statistics

---

# Theory

Electron transfer in proteins occurs through quantum mechanical tunneling.

The Page–Moser–Dutton model estimates the logarithm of the ET rate constant using
## Equation

$$
\log_{10}(k_{ET}) = 13 - \beta(R-3.6) - \frac{3.1(\Delta G+\lambda)^2}{\lambda}
$$

where:
where

| Symbol       | Meaning                        |
| ------------ | ------------------------------ |
| $k_{ET}$   | Electron transfer rate (s⁻¹) |
| $R$        | Donor–acceptor distance (Å)  |
| $\beta$    | Distance decay constant        |
| $\Delta G$ | Free energy difference (eV)    |
| $\lambda$  | Reorganization energy (eV)     |

---

## Distance Decay Constant

The script uses

$$
\beta = 1.2 - 0.8\rho
$$

with

$$
\rho = 0.76
$$

therefore

$$
\beta = 0.592\ \AA^{-1}
$$

---

## Distance Measurement

For every donor–acceptor pair,

the script computes

$$
R =
\min
\left(
d_{ij}
\right)
$$

where

$$
d_{ij}
=
\|
\mathbf r_i-\mathbf r_j
\|
$$

is the Euclidean distance between every atom in the donor residue and every atom in the acceptor residue.

This represents the **minimum edge-to-edge distance**, rather than centroid-to-centroid separation.

---

## Electron Transfer Rate

After computing

$$
\log_{10}(k_{ET}),
$$

the script converts it into the physical rate

$$
\boxed{
k_{ET}=10^{\log_{10}(k_{ET})}
}
$$

with units

$$
s^{-1}
$$

---

# Chromophore Residues

The current implementation considers

| Residue | Role                    |
| ------- | ----------------------- |
| TRP     | Electron donor/acceptor |
| TYR     | Electron donor/acceptor |
| HIS     | Electron donor/acceptor |
| MET     | Electron donor/acceptor |
| CYS     | Electron donor/acceptor |

All atoms belonging to these residues are included when determining edge distances.

---

# ΔG Values

The script contains experimentally derived free energy differences for every donor–acceptor pair.

Example:

| Donor | Acceptor | ΔG (eV) |
| ----- | -------- | -------- |
| TRP   | TYR      | 0.3527   |
| TYR   | TRP      | −0.3527 |
| TRP   | MET      | 0.3459   |
| TYR   | HIS      | 0.0740   |

The complete lookup table is stored in

```python
DELTA_G_VALUES
```

---

# Reorganization Energies

Likewise, pair-specific reorganization energies are stored in

```python
LAMBDA_VALUES
```

Example

| Donor | Acceptor | λ (eV) |
| ----- | -------- | ------- |
| TRP   | TYR      | 2.290   |
| TRP   | MET      | 2.361   |
| TYR   | HIS      | 4.294   |

---

# Workflow

```text
PDB Structure
      │
      ▼
Read Protein
      │
      ▼
Find Chromophore Residues
      │
      ▼
Extract Atomic Coordinates
      │
      ▼
Compute Edge Distances
      │
      ▼
Construct Distance Matrix
      │
      ▼
Apply PMD Equation
      │
      ▼
Generate ET Rate Matrix
      │
      ▼
Summary Statistics
      │
      ▼
Top Electron Transfer Pairs
      │
      ▼
Heatmap
```

---

# Output Files

Running the script on

```text
protein.pdb
```

produces

## 1. Chromophore Atom List

```text
protein_chromophore_atom_list.csv
```

Contains

- residue
- chain
- residue number
- atom names
- number of atoms

---

## 2. Pairwise Distance Table

```text
protein_allatoms_edge_to_edge_pairwise_table.csv
```

Contains

- donor
- acceptor
- donor atom
- acceptor atom
- minimum distance

---

## 3. Distance Matrix

```text
protein_allatoms_edge_to_edge_r_matrix.csv
```

Square distance matrix

Rows

```text
Donors
```

Columns

```text
Acceptors
```

Values

```text
Minimum edge distance (Å)
```

---

## 4. PMD Rate Matrix

```text
protein_allatoms_edge_to_edge_r_matrix_pmd_rates.csv
```

Contains

- ET rate matrix
- summary statistics
- fastest pathways
- ET category counts

---

## 5. Heatmap

```text
protein_allatoms_edge_to_edge_r_matrix_pmd_rates_heatmap.png
```

Heatmap of

$$
\log_{10}(k_{ET})
$$

---

# Electron Transfer Categories

The script classifies ET rates as

| Category   | log₁₀(kET) |
| ---------- | ------------ |
| Ultrafast  | ≥12         |
| Fast       | 10–12       |
| Moderate   | 7–10        |
| Slow       | 4–7         |
| Negligible | <4           |

---

# Installation

```bash
pip install numpy pandas biopython matplotlib seaborn
```

---

# Usage

## Analyze a PDB Structure

```bash
python pmd_rates_simple.py protein.pdb
```

---

## Analyze an Existing Distance Matrix

```bash
python pmd_rates_simple.py distance_matrix.csv
```

---

# Example

Input

```text
1abc.pdb
```

Output

```text
1abc_chromophore_atom_list.csv

1abc_allatoms_edge_to_edge_pairwise_table.csv

1abc_allatoms_edge_to_edge_r_matrix.csv

1abc_allatoms_edge_to_edge_r_matrix_pmd_rates.csv

1abc_allatoms_edge_to_edge_r_matrix_pmd_rates_heatmap.png
```

---

# Project Structure

```text
pmd_rates_simple.py
│
├── PDB Parsing
├── Chromophore Detection
├── Distance Calculations
├── Pairwise Distance Table
├── Distance Matrix Generation
├── PMD Rate Calculation
├── Statistical Summary
├── Top ET Pair Ranking
└── Heatmap Generation
```

---

# Dependencies

- Python ≥ 3.8
- NumPy
- Pandas
- Biopython
- Matplotlib
- Seaborn

---

# Computational Complexity

For **N chromophoric residues**, the algorithm computes all donor–acceptor combinations:

- Distance calculations: **O(N²)**, with each pair evaluating all atom–atom combinations between residues.
- PMD rate calculations: **O(N²)**.
- Memory usage: **O(N²)** for storing the pairwise distance and rate matrices.

---

# Citation

If you use this software in academic work, please cite the original **Page–Moser–Dutton** electron transfer model together with this implementation.
