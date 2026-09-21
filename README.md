# AMP Challenge 2027 — Filter-First Generative AMP Design

A fully reproducible package for the [AMP Challenge 2027](https://szczurek-lab.github.io/amp-challenge-website/) (NeurIPS Competition Track): a 50,000-peptide antimicrobial library, a ranked top-100 candidate list, the code that produced them, calibrated surrogate models, an independent three-judge validation panel, and a paper draft.

**Thesis:** *filter-first.* Sequence quality is set by the ranking stage, not by generator complexity. We generate broadly on the manifold of real AMPs with a simple, interpretable model, then invest the effort in a well-calibrated, independently-validated ranking and selection stage. The whole pipeline runs on one workstation at ~$0 compute.

## Headline results

| | |
|---|---|
| Library | 50,000 unique peptides, 8–50 aa, 20 canonical AA, 0 exact matches to MarLys-AMP |
| Distribution realism | Fréchet Biological Distance **3.0** vs 94.2 for a naive baseline (31× closer) |
| Top-100 potency | median predicted *E. coli* MIC **~1.8 µM** (strongest independent judge), ~10× below the 16 µM threshold |
| Independent validation | 3 held-out judges (MBC-Attention PCC 0.71, AMPredictor 0.50, Pandi DeepAMP); 92–100 % of the list beats the real-AMP median |
| Novelty | 0 exact matches vs MarLys / DBAASP / UniProt; worst identity 0.77 (< 0.80 threshold) |
| Reproducibility | library regenerates **byte-for-byte** from a fixed seed; passes the organizers' `verify_submission.py` |

See the [paper draft](paper/AMP_Challenge_paper.md) for the full method and figures.

## Repository map

```
├── paper/                     # paper draft (Markdown) + figures/
├── submission/                # the competition deliverables
│   ├── ClaudeAMP_library_50k.fasta   # 50,000-sequence library
│   ├── ClaudeAMP_top100.fasta        # ranked top-100 (broad-spectrum, submitted)
│   ├── ABSTRACT_1500.txt / METHOD_SUMMARY_10k.txt / SUBMISSION_FORM_ANSWERS.md
│   ├── WHY_LOCK_NOW.md               # the lock/resource-allocation rationale
│   └── category_lists/               # 5 per-category top-100 lists (design exploration)
├── src/
│   ├── amp_challenge_2027_pkg/  # byte-reproducible generator package (uv/pyproject)
│   ├── generate/                # generation + library scoring scripts
│   ├── eval/                    # compliance + seqme evaluation harness
│   ├── surrogate/               # ESM-2 + XGBoost ensemble (weight-fetch script)
│   ├── benchmark/               # the 3 independent judges (AMPredictor, MBC, Pandi)
│   └── optimize/                # the (rejected) DeepAMP optimization pass
├── checkpoints/               # SMALL model files: markov_model.json, ensemble.pkl
├── results/                   # all benchmark figures (*.png) + results/json/*.json
├── docs/                      # analysis & planning notes
├── env/                       # environment.yml + requirements.txt
└── LICENSE · CITATION.cff · .gitignore
```

## Quickstart — regenerate the library

```bash
# 1. environment (Python 3.11)
conda env create -f env/environment.yml   # or: pip install -r env/requirements.txt

# 2. regenerate the 50,000-sequence library byte-for-byte (pure Python/NumPy, seed 2027)
cd src/amp_challenge_2027_pkg
pip install -e .
generate_broad_spectrum                    # writes the library + top.fasta

# 3. verify compliance with the organizers' checker
python ../../src/eval/eval_harness_compliance.py <library.fasta> <top.fasta>
```

The generator needs only `numpy` and the small `checkpoints/markov_model.json`. Reproducing the *ranking* additionally needs the surrogate stack (ESM-2, XGBoost); reproducing the *independent benchmark* needs the three judge repos and their weights, fetched by the scripts in `src/benchmark/` and `src/surrogate/`.

## What is intentionally not vendored

Large model weights are fetched by scripts, not committed (see `.gitignore`):

- **ESM-2 / ESM-1b** — `src/surrogate/fetch_weights.py` (Hugging Face; set `HF_HUB_DISABLE_XET=1`).
- **AMPredictor, MBC-Attention, Pandi Deep_AMP** — clone the upstream repos (URLs in the paper references); wrappers in `src/benchmark/` reuse their shipped weights.
- **DeepAMP (Li 2024) generator** — weights on the authors' Google Drive; used only by the *rejected* optimization pass in `src/optimize/`.

## Method in one paragraph

An order-3 residue Markov model trained on ~40,000 real AMPs (MarLys-AMP + DBAASP peptides with MIC ≤ 16 µM) generates a length- and charge-matched 50,000-peptide library. Candidates are ranked by an ESM-2 + XGBoost surrogate ensemble calibrated on a homology-aware split, keeping only the four per-strain MIC heads that reach PCC ≥ 0.5. The top-100 is selected for potency, multi-strain coverage, safety (HC50), diversity, and hard novelty/compliance gates, then validated by three activity predictors not used in selection. A generative-optimization pass that improved the list under its own acceptance gate was **rejected** after two held-out judges showed the gain did not transfer — see §3.5 of the paper.

## License

MIT (see `LICENSE`). All training data is open (MarLys-AMP CC-0; DBAASP via QMAP). Prepared for co-authorship eligibility: code, small checkpoints, and training-data provenance are all open.
