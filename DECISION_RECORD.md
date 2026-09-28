# Re-evaluation decision record — replacing the locked submission

## Decision

**Replace** the locked 50,000-sequence library and top-100. The previous
`WHY_LOCK_NOW.md` argued for freezing on the grounds that three independent
activity predictors agreed the top-100 was strong and that a DeepAMP
optimisation pass had failed its held-out judges. Both statements were true and
remain true. Neither addressed the library, and the library is what the
Aggregation Score is computed on.

## What changed the decision

We reproduced the actual Phase-1 scorer — `seqme` (named in the rules) plus the
three AMP oracles in its third-party plugin registry (amPEPpy, AMPredictor,
AMPlify) — and scored the locked library against the contrast set the rules say
the score was tuned to discriminate: real potent AMPs, weak AMPs, UniProt
non-AMPs, and character-permuted AMP decoys.

The locked library was at or below shuffled-decoy level on three of four
potency/realism components and tied on the fourth:

| component | locked | decoy | UniProt | potent AMPs |
|---|---|---|---|---|
| amPEPpy | 0.528 | 0.564 | 0.232 | 0.758 |
| AMPlify | 0.577 | 0.565 | 0.103 | 0.929 |
| AMPredictor MIC (µM) | 41.59 | 33.95 | 33.56 | 15.58 |
| FBD → potent AMPs | 3.64 | 3.04 | 3.55 | 0.01 |

It was strong on Diversity (0.855), Novelty (1.000), Uniqueness (1.000),
Conformity (0.392) and Authenticity (0.779) — but those were never the binding
components, and the permuted-decoy cohort scores *higher* on diversity (0.857)
than any real cohort, so a high diversity number was not evidence of quality.

Root cause: the library was raw generator output — the first 50,000 compliant
draws from a ~70,000 pool. Only the top-100 had ever been selected, and only
from 5,000 pre-filtered candidates. An order-3 Markov chain over AMPs reproduces
local composition but not global amphipathic architecture, so to a
whole-sequence oracle it resembles the decoy cohort by construction.

## What was done

Same generator, same corpus, same seed 2027; pool enlarged to 1,000,000 draws
(998,452 compliant/unique/novel). Library selected as the length-stratified
top-50,000 by mean percentile rank of amPEPpy and AMPredictor, with AMPlify
excluded from selection.

## Acceptance rule and outcome

Because the aggregation weights are withheld, we required no regression on any
named seqme component. Outcome: **8 improved, 2 unchanged, 1 regressed**
(`no_regression_decision.csv`).

The single regression is Diversity, 0.855 → 0.816. We accepted it because:

- 0.816 is effectively the diversity of real potent AMPs (0.820), i.e. the
  realism target rather than a deficit;
- the frontier sweep (`frontier.csv`) shows the trade is near-unavoidable —
  softening selection until amPEPpy falls back to 0.758 and median MIC rises to
  10.6 µM recovers only +0.011 diversity;
- the components that improved (predicted potency, embedding distance to potent
  AMPs, physicochemical conformity) are the ones the rules state the Aggregation
  Score was explicitly tuned to detect.

## The control that justifies the change

AMPlify had **no role** in library selection and improved 0.577 → 0.907 against
a real-potent-AMP reference of 0.929. This is the same test that caused the
earlier DeepAMP pass to be rejected — there, the held-out judges got worse. Same
test, opposite result, so the opposite decision follows.

Mean net charge moved 3.89 → 6.05 against a real-potent-AMP value of 5.90:
toward the reference distribution, not past it.

## Top-100

Four-judge rank consensus (amPEPpy, AMPredictor, AMPlify, in-house
ESM-2/XGBoost strain-MIC ensemble) over the 5,000 highest-consensus library
members, with hard gates (HC50 ≥ 16 µM; ≤ 0.70 identity to any already-selected
candidate; ≤ 0.79 identity to every reference sequence).

Both lists re-scored through the identical pipeline in this session
(`top100_comparison.csv`), so the comparison is like-for-like:

| | locked top-100 | new top-100 |
|---|---|---|
| amPEPpy (mean) | 0.727 | 0.905 |
| AMPlify (mean) | 0.961 | 1.000 |
| AMPredictor median MIC | 2.75 µM | 1.11 µM |
| median predicted HC50 | 99.4 µM | 44.6 µM |
| median predicted safety window | 35.1× | 45.7× |
| worst identity to reference | — | 0.743 |

Predicted HC50 falls, but the selectivity window — the quantity that matters —
improves, because MIC falls faster. A separate Optimal-Selectivity list ranks on
the window directly: median HC50 199.9 µM, median MIC 0.69 µM, median window
276×.

(Note: a prior-session artifact recorded the locked list's AMPredictor median as
3.66 µM. The 2.75 µM above is the value obtained here under the same batching and
settings used for the new list, and is the one used for comparison.)

Leave-one-judge-out audit (`loo_audit.csv`): dropping any single judge changes
34–58% of membership, so no judge dominates and the four carry independent
information. The four-judge list beats the 5,000-candidate pool average on all
four judges simultaneously (0.905 vs 0.872; 1.11 vs 1.997 µM; 1.000 vs 0.937;
0.971 vs 0.846).

## Compliance issue found and fixed

The previous repository had **no root `pyproject.toml`** — the uv project was
nested one directory down. The organizers' `verify_submission.py` runs `uv sync`
at the repository root, so verification would have failed on the locked
submission regardless of sequence quality. The repository is now a root-level uv
project with all five category entry points and a committed `uv.lock`.

Verified: all five entry points run from a clean clone, emit identical
libraries, and pass the organizers' own `_verify_sequences`,
`_verify_no_overlap`, `_verify_top` and `_veritfy_max_simularity` checks.
`library.fasta` sha256
`8a1535f8e509fcb6635ef9be789f842fee9356c0fa61bd495b4fb9f4d0842793`, identical
across repeated runs. `--verify-selection` re-derives the library from the
shipped oracle scores and matches.

## Not assessed

- **Wet-lab activity.** Everything above is *predicted*. The oracles agree with
  each other and with real-AMP reference distributions, which is the strongest
  computational evidence available, but no peptide here has been synthesised.
- **Aggregation Score itself.** The component weights are withheld, so we
  optimise for monotone improvement across components rather than for the score.
  We cannot state a rank.
- **Oracle calibration on our own sequences.** The oracles were trained on
  natural and literature AMPs; their absolute MIC predictions on de-novo Markov
  sequences carry unquantified extrapolation error. We rely on them
  comparatively (our library vs. real cohorts scored identically), not
  absolutely.
- **Hemolysis beyond the in-house HC50 head** (Pearson r 0.51 on a
  homology-aware split) — a weak predictor used only as a gate and a tiebreak.

## Outstanding before upload

1. Replace the `ClaudeAMP` filename/team-name placeholder with the registered
   Kaggle team name.
2. Push the updated repository and confirm the public clone URL resolves.
3. Register with an institutional email if co-authorship eligibility is wanted.
4. Submit before the deadline (October 1, 2026 AOE).
