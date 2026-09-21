# AMP Challenge 2027 — Submission Package

A reproducible, zero-cost generative pipeline for the AMP Challenge 2027
(NeurIPS 2026 Competition Track). Produces a 50,000-sequence antimicrobial
peptide library plus a surrogate-ranked top-100 for each of five categories.

## Method (filter-first, on-manifold, local-only)

**Generator.** An order-3 residue Markov model trained on 40,047 real AMPs
(the MarLys antibacterial subset + QMAP peptides with panel MIC ≤ 16 µM).
Length targets are drawn from the training length distribution and the
end-token is honored only at/after the target, so output length matches real
AMPs (median 18 aa, mean net charge +2.9 — vs MarLys 18 / +2.8). On-manifold
by construction; no GPU, no fine-tuning, no network at generation time.

**Filter.** A surrogate ensemble (ESM-2 650M embeddings + physicochemical
features) scores every candidate: a potency classifier (active = MIC ≤ 16 µM
on any panel strain) and per-strain XGBoost MIC regressors, all calibrated on
QMAP's homology-aware split. Only the four strain heads that clear PCC ≥ 0.5
(E. coli, S. aureus, K. pneumoniae, S. enterica) drive selection; low-data
heads are excluded rather than trusted blindly. An HC50 regressor supplies a
soft safety signal.

**Selection.** Each category applies its own objective to the shared scored
pool, then greedily builds a top-100 under two identity gates —
Levenshtein ratio ≤ 0.8 (the shipped verifier's rule) AND MMseqs2 identity
< 0.8 (the competition PDF's rule) vs MarLys — plus an internal diversity
cap. Because 25 of the top-100 are drawn uniformly at random for wet-lab
assay, selection maximizes the list floor, not a hand-picked frontier.

## Categories

| Category | Objective |
|---|---|
| `broad_spectrum` | potency + coverage across all 4 reliable strains + safety |
| `gram_neg` | activity across the 3 reliable Gram-negative strains |
| `gram_pos` | S. aureus activity (the reliable Gram-positive head) |
| `mdr` | potency + Gram-negative coverage (MDR panel is GN-heavy) |
| `therapeutic` | selectivity: high HC50 with low MIC (safety window) |

## Quality (vs the starter-kit toy K/P baseline)

- **Compliance:** all 5 libraries + top-100 lists pass every gate (unique,
  8–50 aa, canonical alphabet, no exact MarLys match, both identity gates).
- **Diversity (seqme):** 0.86 vs 0.33 toy.
- **Realism (FBD to MarLys, ESM-2 space):** 3.0 vs 94.2 toy — 31× closer to
  the real-AMP distribution.
- **Top-100 predicted activity:** median potency 0.92–0.98; ≥97% of each
  list predicted active (MIC ≤ 16 µM) on the reliable strains, except the
  therapeutic list which deliberately trades S. aureus coverage (55%) for a
  4× higher safety window (median HC50 211 vs ~50 µM).

## Reproducibility

```
uv sync
uv run generate_broad_spectrum   # → generate_broad_spectrum/{library,top}.fasta
uv run generate_gram_pos
uv run generate_gram_neg
uv run generate_mdr
uv run generate_therapeutic
```

Fixed seed (2027), pure Python/NumPy. Re-running reproduces `library.fasta`
byte-for-byte (verified: identical SHA-256 across runs).

## Caveats

Surrogate scores are self-consistent predictions (PCC 0.5–0.64 heads), not
ground truth — the library is only as strong as those heads, which is why the
design also rewards realism and diversity, not activity alone. The Markov
generator is a strong baseline, not a ceiling; a fine-tuned generator is an
optional later upgrade.
