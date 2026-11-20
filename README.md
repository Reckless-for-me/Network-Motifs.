# Network-Motifs.
Network Motif Analysis using Python and Jupyter Notebook.
# Network Motifs Analysis

## Project Overview

This project analyzes network motifs in the SNAP Email-Eu-core dataset.

A network motif is a small connection pattern that occurs significantly more often in a real network than in comparable randomized networks.

The project studies email communication between people and identifies recurring patterns of communication using directed triads.

## Dataset

**Dataset:** SNAP Email-Eu-core

The dataset contains:
- 1,005 nodes representing users.
- 25,571 directed edges representing email communication.
- 42 departments.

Self-loops and duplicate directed edges are removed during data cleaning.

Dataset source: https://snap.stanford.edu/data/email-Eu-core.html

## Project Objectives

1. Load and clean the email communication dataset.
2. Create department-specific induced subgraphs.
3. Identify the 13 connected directed triad classes.
4. Count connected triads in each department.
5. Generate randomized networks while preserving in-degree and out-degree sequences.
6. Calculate Z-scores for each triad class.
7. Identify motifs and anti-motifs.
8. Compare the results of two selected departments using a grouped bar chart.

## Methodology

### 1. Triad Classification

Each group of three connected nodes is classified into one of 13 connected directed triad classes.

Canonical codes:

`[3, 6, 7, 10, 11, 15, 21, 23, 25, 27, 30, 31, 63]`

### 2. Triad Census

The project counts connected triads using pairs of neighbours rather than checking every possible combination of three nodes.

### 3. Randomized Networks

Random networks are generated using directed edge swaps.

The randomization preserves the in-degree and out-degree of every node.

### 4. Z-Score Calculation

The Z-score measures how much the real triad count differs from the average count in randomized networks.

\[
Z = \frac{N_{\text{real}}-\mu_{\text{random}}}{\sigma_{\text{random}}}
\]

Where:

- \(N_{\text{real}}\) = triad count in the original network.
- \(\mu_{\text{random}}\) = mean triad count in randomized networks.
- \(\sigma_{\text{random}}\) = standard deviation of randomized triad counts.

Classification:

- Z-score greater than 2: Motif.
- Z-score less than -2: Anti-motif.
- Otherwise: No significant pattern.

If the standard deviation is zero, the Z-score is set to zero.

## Assignment Parameters

The following parameters are calculated from the group roll numbers:

| Parameter | Value |
|---|---:|
| R | 2073 |
| Selected pair index | 173 |
| Number of randomized networks | 47 |
| Randomization factor | 8 |
| Random seed | 2073 |

The number of swap attempts for each randomized network is:

`8 × number of internal department edges`

## Technologies Used

- Python
- Jupyter Notebook
- Matplotlib
- NetworkX (optional verification only)

## Requirements and Installation

Install the required libraries:

```bash
pip install -r requirements.txt
```

Launch Jupyter Notebook:

```bash
jupyter notebook
```

Open `Network_Motifs.ipynb` and run the cells in order.

## Project Structure

```text
Network_Motifs/
│
├── Network_Motifs.ipynb
├── README.md
├── requirements.txt
│
├── data/
│
└── results/
```

## Output

The notebook generates:

- Cleaned network data.
- Department-specific network information.
- Real triad counts.
- Randomized triad counts.
- Z-scores for all 13 triad classes.
- Motif and anti-motif classifications.
- Department comparison bar chart.

## Conclusion

This project uses network motif analysis to study patterns in email communication. By comparing real networks with randomized networks that preserve degree sequences, it identifies statistically significant communication patterns between users.
