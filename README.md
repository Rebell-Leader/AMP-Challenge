# AMP Challenge 2027 — oracle-selected Markov AMP library

**Submission v2.** A simple, interpretable generator samples on the manifold of
real antimicrobial peptides; the 50,000-member library is then *selected* on the
published AMP oracles the competition's own evaluation framework ships.

The reason this is v2 is worth stating up front: we reproduced the Phase-1
scorer locally and found our own v1 library scoring **at the shuffled-AMP decoy
floor** on predicted potency and embedding realism, while looking excellent on
diversity and novelty — which were never the binding components. See
[DECISION_RECORD.md](DECISION_RECORD.md).

## Quick start

```bash
uv sync
uv run generate_broad_spectrum                     # library.fasta + top.fasta
uv run generate_broad_spectrum --verify-selection   # also re-derive the selection
```

CPU-only, no network, a few minutes. Output is byte-identical on every run
(`library.fasta` sha256 `8a1535f8e509fcb6635ef9be789f842fee9356c0fa61bd495b4fb9f4d0842793`).

## Entry points

| entry point | category | top-100 ranking |
|---|---|---|
| `generate_broad_spectrum` | Broad spectrum | potency consensus |
| `generate_gram_pos` | Optimal activity, Gram-positive | potency consensus |
| `generate_gram_neg` | Optimal activity, Gram-negative | potency consensus |
| `generate_mdr` | Optimal activity, MDR / WHO-priority | potency consensus |
| `generate_therapeutic` | Optimal selectivity | safety window (HC50/MIC) |

All five emit the same 50,000-sequence library. The four activity categories
share the potency-ranked top-100; `generate_therapeutic` ships a separate list
ranked on predicted selectivity.

## Method in brief

1. **Generate.** Order-3 residue Markov model trained on 40,047 real AMPs
   (MarLys antibacterial set + DBAASP/QMAP peptides with panel MIC ≤ 16 µM),
   seed 2027, 1,000,000 candidates → 998,452 after compliance filtering.
2. **Select the library.** Length-stratified top-50,000 by mean percentile rank
   of **amPEPpy** and **AMPredictor**. **AMPlify is held out of selection** and
   used only as an independent judge.
3. **Rank the top-100.** Four-judge consensus (amPEPpy, AMPredictor, AMPlify,
   and an in-house ESM-2/XGBoost per-strain MIC ensemble) under hard gates:
   predicted HC50 ≥ 16 µM, ≤ 0.70 identity to any already-selected candidate,
   ≤ 0.79 identity to every reference sequence.

Oracles are the three AMP models published in the
[seqme third-party plugin registry](https://github.com/szczurek-lab/seqme-thirdparty),
each run through seqme's own `ThirdPartyModel` wrapper.

## Results vs. the v1 library

Eight of eleven seqme components improved, two unchanged, one regressed
(Diversity 0.855 → 0.816, landing at the real-potent-AMP level of 0.820).
Held-out AMPlify rose **0.577 → 0.907** against a real-potent-AMP reference of
0.929, so the gain is not an artifact of optimising against the selecting
models.

Full tables: [DECISION_RECORD.md](DECISION_RECORD.md),
[METHOD_SUMMARY_10k.txt](METHOD_SUMMARY_10k.txt).

## Reproducibility

Oracle inference needs a GPU and three mutually incompatible legacy
environments, so per-candidate oracle scores ship as a checkpoint
(`checkpoint/selection.npz`, float64) exactly as trained weights would. The
candidate pool itself regenerates from the seed on every run. The scoring code
that produced the checkpoint is [`src/score_pool.py`](src/score_pool.py), and
`--verify-selection` re-derives the library from the shipped scores and asserts
it matches the emitted FASTA — so the selection step is auditable, not merely
asserted.

## Layout

```
pyproject.toml               root uv project (five entry points)
uv.lock                      committed lockfile
src/amp_challenge_2027/      generator + selection
src/score_pool.py            regenerates the oracle-score checkpoint (GPU)
checkpoint/                  markov_model.json, selection.npz, top100_scores.npz
data/antibacterial.fasta     competition reference set
submission/                  final FASTAs + compliance report
DECISION_RECORD.md           why v1 was replaced
METHOD_SUMMARY_10k.txt       full method disclosure
KAGGLE_WRITEUP.md            writeup form content
```

## License

MIT — see [LICENSE](LICENSE). Code, weights and training-data provenance are all
open, so nothing is withheld under the challenge's full-requirements disclosure
rule.
