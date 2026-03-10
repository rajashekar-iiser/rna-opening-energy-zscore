# RNA Opening Energy Z-score Calculator

A Python implementation for calculating **RNA motif opening energy z-scores** using the **ViennaRNA package** and **dinucleotide shuffling**.

This tool was developed as an in-house implementation of the methodology used in the Perl utility:

https://github.com/mtw/plfoldz

The script computes **tetranucleotide opening energies** and evaluates statistical significance using **1000 dinucleotide shuffling events**.

---

## Features

* ViennaRNA-based thermodynamic calculations
* Dinucleotide-preserving sequence shuffling
* Z-score calculation for motif accessibility
* Automated processing of multiple FASTA and BED inputs
* Generates analysis-ready CSV outputs

---

## Requirements

Python ≥ 3.8

Dependencies:

* ViennaRNA
* Biopython
* numpy
* pandas
* tqdm
* matplotlib
* seaborn

Install using:

```
pip install -r requirements.txt
```

---

## Input

### FASTA

Genome or sequence FASTA files:

```
*_ungapped.fasta
```

### BED

Motif coordinates:

```
chrom   start   end
```

---

## Usage

Run in the directory containing the FASTA files and motif directories.

```
python z_score_calculator.py
```

---

## Output

The script generates:

```
z_scores_df.csv
z_scores_df_melt.csv
seq_data_df.csv
```

---

## Method Summary

1. Extract motif sequences from BED coordinates
2. Add ±100 nt flanking regions
3. Calculate RNA opening energy using ViennaRNA
4. Perform 1000 dinucleotide shuffles of flanking sequences
5. Compute z-score relative to shuffled background distribution

---

## Citation

If you use this code, please cite the associated manuscript.

---

## License

MIT License
