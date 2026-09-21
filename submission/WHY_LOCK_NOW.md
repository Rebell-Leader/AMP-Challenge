# Why we lock now — and what it would actually take to break through

## The decision
We lock the submission now. This is not "we ran out of time" — it is a
resource-allocation decision grounded in what the benchmarks told us.

## What the benchmarks established
Our top-100 was validated by **four independent activity predictors**: our own
ESM-2/XGBoost ensemble (used for selection) and three surrogates named in the
competition that had **no role** in selection — AMPredictor (ESM-1b+GNN,
PCC 0.50 vs ground truth), MBC-Attention (iFeature+CNN, PCC 0.71), and Pandi
DeepAMP (one-hot CNN/LSTM).

- All four agree the list is **substantially more potent than typical natural
  AMPs**: 83–100% of the top-100 beats the median real AMP under each judge.
- The **best-validated** independent judge (MBC-Attention, PCC 0.71) predicts a
  **median E. coli MIC of ~2 µM** — roughly an order of magnitude below the
  competition's 16 µM potency threshold.
- The list is fully compliant, reproducible byte-for-byte, and passes the
  organizers' own `verify_submission.py` end-to-end.

The convergence of decorrelated judges is the key signal. It means the quality is
**real**, not an artifact of scoring peptides with the same model that selected
them. That is the single most important thing to establish before a one-shot
submission, and we have established it.

## Why pushing further is low expected value
We *tried* to push further, so this is empirical, not a guess. We ran a bounded
DeepAMP (Li 2024) seed-and-optimize pass over the top-100 using the strongest
judge, MBC-Attention, as the acceptance gate. Under MBC-Attention it looked like a
big win (median MIC roughly halved, 71–87% of seeds "improved"). But the two
**held-out** judges told the real story:

| Judge | Role | Verdict on "optimized" list |
|---|---|---|
| MBC-Attention | optimization target | better (as expected — it was the gate) |
| Pandi DeepAMP | held-out | **worse in 8/10** category×head comparisons |
| AMPredictor | held-out | **worse in 5/5** categories |

The gain lived entirely in the model we optimized against. This is textbook
surrogate overfitting, and shipping it would have **degraded** the submission.
We rejected it and kept the original list.

The lesson generalizes: our activity predictors share a **~0.5 PCC ceiling** with
every published AMP MIC model (including the competition's own named surrogates).
That ceiling — not the generator — is the binding limit. No amount of generator
cleverness moves a candidate past a scorer that can only rank at 0.5 correlation.
And the exact Phase-1 ranking weights are withheld by the organizers, so we cannot
tune to the actual objective even if we wanted to.

## What a genuine breakthrough would require (and cost)
To move the needle we would have to raise the **activity-prediction ceiling**, not
the generator. Concretely:

1. **Better MIC data for the hard strains.** Six of ten panel strains
   (P. aeruginosa, A. baumannii, both enterococci, E. cloacae, B. subtilis) are
   label-starved and their heads never reached PCC 0.5. Mining and curating more
   per-strain MIC data (DBAASP full export, dbAMP 3.0, APD, Peptipedia,
   literature) is the highest-value work — but it is a **data** project, not a
   compute project, and its payoff is uncertain because the labels are genuinely
   scarce in the world.

2. **A multi-task / transfer MIC model** (share signal across strains) — a bounded
   training effort, ~1 GPU-day on the local RTX 4080, no cloud needed. This is the
   one item with a plausible positive expected value, and it improves *ranking*,
   not generation.

3. **A LoRA-fine-tuned generator with an ensemble acceptance gate** (all judges
   must agree, over a 10^5–10^6 variant pool). This is the only place cloud/GPU
   would help. Cost estimate: a single A40/A100 spot instance for <24 h (~$5–30).
   But given that even the honest judges already place our list ~10× under
   threshold, and that our one optimization attempt was rejected for overfitting,
   the expected marginal gain does not justify the added overfitting risk on a
   **single** submission.

## Bottom line
- The submission is validated on every axis we can measure locally, by models we
  did not train against.
- The one attempt to improve it made it worse under independent judges.
- The binding limit is predictor reliability (~0.5 PCC) and the hidden ranking
  weights — neither is fixable by more generation or more compute.
- Rational resource use is therefore **~$0 and lock now**. The only work with
  positive expected value (better strain-MIC data + a small multi-task predictor)
  is a data-curation effort for a *future* round, not a reason to delay this
  submission.
