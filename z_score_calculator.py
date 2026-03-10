#!/usr/bin/env python3
"""
Z-score calculator for RNA opening energy around genomic motifs.

This script computes tetranucleotide opening energies using ViennaRNA
and evaluates statistical significance using dinucleotide shuffling.

Method summary
--------------
1. Extract motif sequences from BED coordinates
2. Add ±100 nt flanking regions
3. Calculate opening energy using ViennaRNA
4. Perform 1000 dinucleotide shuffles of the flanks
5. Compute z-score relative to shuffled background

Outputs
-------
z_scores_df.csv
z_scores_df_melt.csv
seq_data_df.csv
"""

import os
import math
import random
from collections import defaultdict

import numpy as np
import pandas as pd
import tqdm
import RNA
from Bio import SeqIO

# Thermodynamic constants
R = 0.00198717  # kcal/mol/K
T0 = 273.15
T = 37
tempK = T0 + T
kT = R * tempK


# ------------------------------------------------------------
# File Readers
# ------------------------------------------------------------

def read_fasta(fasta_file):
    """Read FASTA file and return dictionary of sequences."""
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))


def read_bed(bed_file):
    """Read BED file."""
    return pd.read_csv(bed_file, sep="\t", header=None)


# ------------------------------------------------------------
# Sequence Utilities
# ------------------------------------------------------------

def shuffle_sequence(seq, dinuc=2):
    """
    Shuffle sequence while preserving dinucleotide structure.
    """
    chunks = [seq[i:i + dinuc] for i in range(0, len(seq), dinuc)]

    for i in range(len(chunks) - 1, 0, -1):
        j = random.randint(0, i)
        chunks[i], chunks[j] = chunks[j], chunks[i]

    return "".join(chunks)


# ------------------------------------------------------------
# RNA Opening Energy
# ------------------------------------------------------------

def compute_opening_energy(sequence, ulength, window_size, max_bp_span, where):
    """
    Compute RNA opening energy using ViennaRNA.
    """
    up = RNA.pfl_fold_up(sequence, ulength, window_size, max_bp_span)
    oe = [-kT * math.log(up[i][ulength]) for i in range(ulength, len(sequence) + 1)]
    return oe[where]


# ------------------------------------------------------------
# Main Z-score Calculation
# ------------------------------------------------------------

def main(fasta_file, bed_file):

    z_scores = []

    fasta_dict = read_fasta(fasta_file)
    bed_df = read_bed(bed_file)

    seq_count = bed_df.shape[0]

    for _, row in bed_df.iterrows():

        chrom, start, end = row[0], row[1], row[2]

        seq = str(fasta_dict[chrom].seq[start:end])

        seq5 = str(fasta_dict[chrom].seq[max(0, start - 100):start])
        seq3 = str(fasta_dict[chrom].seq[end:min(len(fasta_dict[chrom]), end + 100)])

        seq_ori = seq5 + seq + seq3
        motif_length = len(seq)

        motif_oe = compute_opening_energy(
            seq_ori,
            motif_length,
            100,
            100,
            len(seq5)
        )

        shuf_oe = []

        for _ in range(1000):

            seq5_shuf = shuffle_sequence(seq5, dinuc=2)
            seq3_shuf = shuffle_sequence(seq3, dinuc=2)

            seq_shuf = seq5_shuf + seq + seq3_shuf

            shuf_oe.append(
                compute_opening_energy(
                    seq_shuf,
                    motif_length,
                    100,
                    100,
                    len(seq5_shuf)
                )
            )

        z_score = (motif_oe - np.mean(shuf_oe)) / np.std(shuf_oe)

        z_scores.append(z_score)

    return chrom, seq_count, z_scores


# ------------------------------------------------------------
# Pipeline
# ------------------------------------------------------------

if __name__ == "__main__":

    z_scores = defaultdict(list)
    seq_data = defaultdict(list)

    fasta_files = [
        f for f in os.listdir(".")
        if os.path.isfile(f) and f.endswith("_ungapped.fasta")
    ]

    motifs_dirs = [
        d for d in os.listdir(".")
        if os.path.isdir(d) and d.endswith("_gtags")
    ]

    total_files = len(fasta_files)
    total_motifs = len(motifs_dirs)

    with tqdm.tqdm(total=total_files,
                   desc="Processing fasta files",
                   dynamic_ncols=True) as pbar_files:

        for fasta_file in fasta_files:

            with tqdm.tqdm(total=total_motifs,
                           desc=f"Processing motifs in {fasta_file}",
                           leave=False,
                           dynamic_ncols=True) as pbar_motifs:

                for motif in motifs_dirs:

                    bed_file = (
                        f"{motif}/"
                        f"{fasta_file.split('_ungapped.fasta')[0]}_"
                        f"{motif.replace('_gtags','')}_motifs.bed"
                    )

                    if os.path.exists(bed_file) and os.stat(bed_file).st_size > 0:

                        seq_name, seq_count, z_scores_for_file = main(
                            fasta_file,
                            bed_file
                        )

                        z_scores[motif.replace("_gtags", "")].extend(
                            z_scores_for_file
                        )

                        seq_data[seq_name].append(
                            (motif.replace("_gtags", ""), seq_count, z_scores_for_file)
                        )

                    else:
                        print(f"Bed file {bed_file} does not exist, skipping...")

                    pbar_motifs.update()

            pbar_files.update()


    # ------------------------------------------------------------
    # Output tables
    # ------------------------------------------------------------

    z_scores_df = pd.DataFrame(
        dict([(k, pd.Series(v)) for k, v in z_scores.items()])
    )

    z_scores_df.to_csv("z_scores_df.csv")

    z_scores_df = z_scores_df.melt(
        var_name="motif",
        value_name="z_score"
    )

    z_scores_df.to_csv("z_scores_df_melt.csv")

    seq_data_items = [
        (k, *tup)
        for k, list_of_tuples in seq_data.items()
        for tup in list_of_tuples
    ]

    seq_data_df = pd.DataFrame(
        seq_data_items,
        columns=["sequence_name", "motif", "count", "z_scores"]
    )

    seq_data_df.to_csv("seq_data_df.csv", index=False)