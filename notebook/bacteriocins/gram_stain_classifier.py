#!/usr/bin/env python
"""Estimate Gram stain (positive / negative / indeterminate) for GTDB bacteria.

Gram staining reflects cell-envelope architecture (monoderm vs diderm) and is
strongly conserved at the phylum level, with a handful of well-known class-level
exceptions. There is no Gram column in the GTDB metadata, so we infer it from the
`gtdb_taxonomy` string (phylum + class are enough).

Rules (see PHYLUM_GRAM / CLASS_OVERRIDES below):
  * Monoderm, stain-positive lineages  -> "positive"
      Actinomycetota, Bacillota (Firmicutes, incl. all Bacillota_* splits)
  * Diderm, stain-negative lineages    -> "negative"
      Pseudomonadota and essentially every other characterised diderm phylum
  * Class-level exceptions override the phylum default, e.g. Negativicutes
    (diderm Firmicutes; GTDB phylum Bacillota_C) -> "negative"
  * Anything not covered by curated knowledge -> "indeterminate"
      (candidate/uncultured phyla: Patescibacteria/CPR, alphabet-soup MAGs, and
       lineages whose envelope is genuinely variable, e.g. Chloroflexota,
       Deinococcota)

Limitations: cell-wall-less Mollicutes (Mycoplasma, Phytoplasma, Acholeplasma)
are nested inside Bacillota>Bacilli by GTDB and cannot be separated at the class
level here, so they are (mis)labelled "positive". They are a small minority.

Usage:
    python gram_stain_classifier.py [path/to/gtdb_metadata.csv]
"""

import csv
import os
import sys

# --- Curated Gram classification by GTDB phylum -----------------------------
# "positive" = monoderm / Gram-positive; "negative" = diderm / Gram-negative.
# Phyla not listed here fall through to "indeterminate".
PHYLUM_GRAM = {
    # ---- Gram-positive (monoderm) ----
    "Actinomycetota": "positive",   # Actinobacteria
    "Bacillota":      "positive",   # Firmicutes  (and every Bacillota_* below)
    "Bacillota_A":    "positive",
    "Bacillota_B":    "positive",
    "Bacillota_C":    "positive",   # overridden to negative at class level (Negativicutes)
    "Bacillota_D":    "positive",
    "Bacillota_E":    "positive",
    "Bacillota_F":    "positive",
    "Bacillota_G":    "positive",
    "Bacillota_H":    "positive",
    # ---- Gram-negative (diderm) ----
    "Pseudomonadota":     "negative",  # Proteobacteria
    "Bacteroidota":       "negative",
    "Campylobacterota":   "negative",
    "Spirochaetota":      "negative",
    "Chlamydiota":        "negative",
    "Cyanobacteriota":    "negative",
    "Planctomycetota":    "negative",
    "Verrucomicrobiota":  "negative",
    "Acidobacteriota":    "negative",
    "Desulfobacterota":   "negative",
    "Desulfobacterota_B": "negative",
    "Desulfobacterota_C": "negative",
    "Desulfobacterota_D": "negative",
    "Desulfobacterota_E": "negative",
    "Desulfobacterota_F": "negative",
    "Desulfobacterota_G": "negative",
    "Myxococcota":        "negative",
    "Myxococcota_A":      "negative",
    "Bdellovibrionota":   "negative",
    "Nitrospirota":       "negative",
    "Nitrospirota_A":     "negative",
    "Nitrospinota":       "negative",
    "Fusobacteriota":     "negative",
    "Fibrobacterota":     "negative",
    "Aquificota":         "negative",
    "Thermotogota":       "negative",
    "Deferribacterota":   "negative",
    "Chrysiogenota":      "negative",
    "Synergistota":       "negative",
    "Elusimicrobiota":    "negative",
    "Gemmatimonadota":    "negative",
    "Calditrichota":      "negative",
    "Armatimonadota":     "negative",
    "Methylomirabilota":  "negative",
    "Chlorobiota":        "negative",
    "Ignavibacteriota":   "negative",
    "Coprothermobacterota": "negative",
    "Dictyoglomota":      "negative",
    "Caldisericota":      "negative",
    "Thermodesulfobiota": "negative",
    # ---- Deliberately NOT assigned (left indeterminate) ----
    # Chloroflexota  : monoderm but stains negative -> variable
    # Deinococcota   : diderm but stains positive   -> variable
    # Patescibacteria / CPR and uncultured candidate phyla -> unknown envelope
}

# Class-level exceptions that override the phylum default.
# Keyed by (phylum, class).
CLASS_OVERRIDES = {
    ("Bacillota_C", "Negativicutes"): "negative",  # diderm Firmicutes (Veillonella, Selenomonas, ...)
}


def parse_taxonomy(tax_string):
    """Return a dict of rank-letter -> name, e.g. {'d':'Bacteria','p':'Bacillota',...}."""
    ranks = {}
    for token in tax_string.split(";"):
        token = token.strip()
        if "__" in token:
            rank, name = token.split("__", 1)
            ranks[rank] = name
    return ranks


def classify(tax_ranks):
    """Classify a parsed taxonomy dict as positive / negative / indeterminate."""
    phylum = tax_ranks.get("p", "")
    klass = tax_ranks.get("c", "")
    if (phylum, klass) in CLASS_OVERRIDES:
        return CLASS_OVERRIDES[(phylum, klass)]
    return PHYLUM_GRAM.get(phylum, "indeterminate")


def main():
    default_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "data", "gtdb_metadata.csv",
    )
    path = sys.argv[1] if len(sys.argv) > 1 else default_path

    K = 15
    counts = {"positive": 0, "negative": 0, "indeterminate": 0}
    n_bacteria = 0

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("domain") != "Bacteria":
                continue
            n_bacteria += 1
            ranks = parse_taxonomy(row["gtdb_taxonomy"])
            counts[classify(ranks)] += 1

    # --- Report -------------------------------------------------------------
    print(f"Bacterial genomes analysed: {n_bacteria:,}\n")
    print(f"{'Gram class':<15}{'count':>10}{'fraction':>12}")
    print("-" * 37)
    for label in ("positive", "negative", "indeterminate"):
        n = counts[label]
        print(f"{label:<15}{n:>10,}{n / n_bacteria:>12.4f}")

    p_pos = counts["positive"] / n_bacteria

    print()
    print(f"Drawing K = {K} genomes uniformly at random from the {n_bacteria:,} bacteria,")
    print("the number of Gram-positive draws is Binomial(K, p_pos) (with replacement;")
    print("Hypergeometric without replacement -- same expectation by linearity).")
    print()
    print(f"  P(Gram-positive) = {counts['positive']:,} / {n_bacteria:,} = {p_pos:.4f}")
    print(f"  Expected # Gram-positive = K * p_pos = {K} * {p_pos:.4f} = {K * p_pos:.3f}")

    # Also report the bound treating indeterminate as an upper bound on "positive".
    p_pos_upper = (counts["positive"] + counts["indeterminate"]) / n_bacteria
    print()
    print("  If 'indeterminate' genomes were all counted as Gram-positive (upper bound):")
    print(f"    Expected # Gram-positive <= {K} * {p_pos_upper:.4f} = {K * p_pos_upper:.3f}")

    # --- Significance test: is an observed count of Gram-positive > chance? ---
    # One-sample test of a proportion: null is that the K draws come from the
    # background population (P(Gram-positive) = p_pos). With N = 50,640 and K = 15
    # the draws are effectively independent, so the exact test is Binomial(K, p_pos)
    # (the finite-population hypergeometric is numerically indistinguishable here).
    observed_pos = 11
    from scipy.stats import binomtest, hypergeom

    bt = binomtest(observed_pos, K, p_pos, alternative="greater")
    bt_two = binomtest(observed_pos, K, p_pos, alternative="two-sided")

    # Exact finite-population equivalent (no independence assumption), for comparison.
    # P(X >= observed) drawing K without replacement from n_bacteria with counts['positive'] successes.
    hg_greater = hypergeom.sf(observed_pos - 1, n_bacteria, counts["positive"], K)

    print()
    print(f"Significance of observing {observed_pos}/{K} Gram-positive vs. background p = {p_pos:.4f}:")
    print(f"  Expected under null      : {K * p_pos:.2f}")
    print(f"  Binomial one-sided  P(X >= {observed_pos}) = {bt.pvalue:.3e}")
    print(f"  Binomial two-sided                = {bt_two.pvalue:.3e}")
    print(f"  Hypergeometric one-sided (exact, finite pop) = {hg_greater:.3e}")


if __name__ == "__main__":
    main()
