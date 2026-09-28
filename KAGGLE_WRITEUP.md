# Kaggle writeup — ready-to-paste form content

---

## Title  (limit 80)

```
Scoring our own library with the real benchmark
```
*(47 characters)*

**Alternatives, same constraint:**
- `Our library looked like a decoy. We measured it, then fixed it.` (63)
- `Reproducing the scorer changed our submission` (45)

---

## Subtitle  (limit 140)

```
Reproducing seqme and its three published AMP oracles showed our 50k library sat at the decoy floor. Reselection fixed it.
```
*(122 characters)*

---

## Card / thumbnail image  (560 × 280)

`writeup_card_560x280.png`

---

## Submission track

Computational track (auto-selected).

---

## Project description

**The short version:** we built the benchmark before we trusted our submission,
and the benchmark told us our 50,000-peptide library was indistinguishable from
shuffled-AMP decoys. Fixing that — not improving the generator — is the entire
contribution.

### 1. We reproduced the scorer instead of guessing at it

The rules specify that full libraries are evaluated with **seqme**, and that the
Aggregation Score combines physicochemical properties, predicted potency from
published oracles, embeddings compared against known potent peptides,
synthesizability, novelty and diversity — tuned to separate potent AMPs from weak
ones and from negatives including UniProt and synthetic decoys.

Rather than optimise a proxy, we stood the real thing up locally: `seqme`, plus
all three antimicrobial-peptide oracles published in its third-party plugin
registry (`szczurek-lab/seqme-thirdparty`) — **amPEPpy**, **AMPredictor** and
**AMPlify** — each in its own uv-locked environment through seqme's own
`ThirdPartyModel` wrapper. We then assembled the exact contrast set the rules
describe: real potent AMPs (DBAASP/QMAP, min panel MIC ≤ 16 µM), weak AMPs
(all MIC > 64 µM), UniProt non-AMP negatives, and character-permuted AMP decoys.

### 2. What it said about our own submission

Our library — the one we had already locked, after a careful multi-judge
validation of its top-100 — landed here:

| component | our library | shuffled decoys | UniProt neg. | real potent AMPs |
|---|---|---|---|---|
| amPEPpy | **0.528** | 0.564 | 0.232 | 0.758 |
| AMPlify | **0.577** | 0.565 | 0.103 | 0.929 |
| AMPredictor MIC (µM) | **41.6** | 33.9 | 33.6 | 15.6 |
| FBD → potent AMPs | **3.64** | 3.04 | 3.55 | 0.01 |

Below the decoy floor on three of four potency/realism components, tied on the
fourth, and below even the UniProt negative control on two. Meanwhile it looked
excellent precisely where we had spent our effort: Diversity 0.855, Novelty
1.000, Uniqueness 1.000, Conformity 0.392 — at or above the real-AMP reference.

The reason is structural, and we think it generalises beyond our entry: **an
order-3 Markov chain trained on AMPs is, to an oracle that reads the whole
sequence, close to a smart shuffle of AMPs** — which is exactly what the decoy
cohort is. It reproduces composition and 3-mer statistics but not the global
amphipathic architecture that makes a peptide kill bacteria. High diversity was
not a virtue here: the permuted decoys are the *most* diverse cohort in the
entire panel (0.857).

We had also mis-targeted our filtering. Our top-100 was genuinely strong because
it *was* selected — but from only 5,000 pre-filtered candidates. The library,
which is what decides advancement, had never been selected at all: it was simply
the first 50,000 compliant draws.

### 3. The fix: select the library, don't just generate it

Same generator, same corpus, same seed — we enlarge the draw to **1,000,000
candidates** (998,452 after compliance filtering) so the library can be *chosen*:

1. amPEPpy scores all 998k; a length-stratified top-4× prefilter keeps 200,000.
2. AMPredictor scores those.
3. The library is the length-stratified top-50,000 by mean percentile rank of the
   two oracles. Stratifying against the previous library's length histogram holds
   the length distribution fixed exactly.
4. **AMPlify is excluded from selection entirely** and used only as a judge.

Because the aggregation weights are withheld, we required **no regression on any
named component**. Locked → rebuilt:

| component | locked | rebuilt | verdict |
|---|---|---|---|
| amPEPpy | 0.528 | 0.797 | improved |
| AMPlify *(held out)* | 0.577 | 0.907 | improved |
| AMPredictor MIC (µM) | 41.59 | 9.28 | improved |
| FBD → potent AMPs | 3.64 | 2.89 | improved |
| FBD → UniProt | 2.34 | 4.00 | improved |
| Conformity score | 0.392 | 0.509 | improved |
| Authenticity | 0.779 | 0.788 | improved |
| FKEA | 478.1 | 505.2 | improved |
| Novelty | 1.000 | 1.000 | unchanged |
| Uniqueness | 1.000 | 1.000 | unchanged |
| **Diversity (5)** | **0.855** | **0.816** | **regressed** |

Eight improved, two held, one regressed. We report the regression plainly:
0.816 is essentially the diversity of real potent AMPs (0.820), and a frontier
sweep shows the trade is close to unavoidable — diluting the selection pressure
until amPEPpy falls back to 0.758 and median MIC rises to 10.6 µM buys only
**+0.011** diversity. We judged the potency and embedding-realism gains, which
are the components the Aggregation Score is explicitly built to detect, worth
0.04 of diversity. Mean net charge moves 3.89 → 6.05, against a real-potent-AMP
value of 5.90 — toward the reference, not past it.

**The held-out control is the point.** AMPlify had no role in selection and still
rose 0.577 → 0.907 (real potent AMPs: 0.929). Our previous version had *rejected*
a DeepAMP optimisation pass precisely because the judges we held out then got
worse — textbook surrogate overfitting. Same test, opposite result, opposite
decision.

### 4. Top-100

The 5,000 highest-consensus library members are re-ranked by four judges —
amPEPpy, AMPredictor, AMPlify, and an in-house ESM-2/XGBoost per-strain MIC
ensemble (only the four strain heads reaching Pearson r ≥ 0.5 on a
homology-aware split are used) — with hard gates: predicted HC50 ≥ 16 µM,
≤ 0.70 identity to any already-selected candidate so the random 25-peptide draw
cannot hit near-duplicates, and ≤ 0.79 identity to every reference sequence.

Result: amPEPpy **0.905**, AMPlify **1.000**, AMPredictor median MIC **1.11 µM**,
predicted-active on all four reliable strains, median HC50 44.6 µM, median
predicted safety window **45.7×**, worst reference identity **0.743**. The
four-judge list beats the 5,000-candidate pool average on all four judges
simultaneously — the consensus is not trading one judge against another.

The Optimal-Selectivity entry point ranks instead on predicted safety window
subject to MIC ≤ 16 µM: median HC50 199.9 µM, median MIC 0.69 µM, median window
**276×**.

### 5. Reproducibility and compliance

- Repository is a **uv project at its root** with all five category entry points,
  a committed `uv.lock` and a pinned Python version, so `uv sync` then
  `uv run generate_broad_spectrum` works from a clean clone. *(Our previous
  layout nested the uv project one level down — the organizers' `uv sync` step
  would have failed on it. Worth checking in your own repo.)*
- Oracle inference needs a GPU and three mutually incompatible legacy
  environments, so per-candidate oracle scores ship as a checkpoint, exactly as
  trained weights would. The 1M pool itself regenerates from seed 2027 on every
  run, **CPU-only, in minutes**. Two runs give byte-identical output.
- `--verify-selection` re-derives the library from the shipped oracle scores and
  asserts it matches the emitted FASTA, so the selection step is auditable rather
  than merely asserted.
- Passes the organizers' own `verify_submission.py` checks for all five entry
  points: 50,000 unique 8–50 aa canonical sequences, zero exact matches to the
  reference set, top-100 a unique library subset, no top-100 sequence above 0.80
  identity to any reference peptide.
- Total compute: one workstation, one consumer GPU, ~3 hours, no cloud spend.
- All data open (MarLys-AMP, DBAASP via QMAP, UniProt); no proprietary inputs, so
  nothing is withheld under the full-requirements disclosure rule.
- No manual sequence curation whatsoever: no peptide was hand-picked,
  hand-edited, or dropped by inspection. Every filter is a scripted threshold in
  version control.

### 6. What we'd tell other teams

1. **Stand up the real scorer first.** It is open, it installs, and it will
   probably surprise you. Ours took a few hours and changed our submission.
2. **Check your library, not just your shortlist.** The library decides who
   advances, and it is the easiest thing to leave unselected.
3. **Include a decoy cohort.** Shuffled versions of your own training AMPs are
   the single most informative baseline we ran; "better than random" is a much
   weaker claim than "better than a shuffle".
4. **Diversity can be a symptom.** Ours was high because our output was close to
   noise. The right target is the diversity of real potent AMPs, not the maximum.
5. **Always hold a judge out.** It is the only thing separating a real gain from
   scorer overfitting — and it told us to reject our last optimisation and accept
   this one.

---

## Attachments / project links

- GitHub repository (code, weights, checkpoints, uv.lock, docs)
- `reevaluation_panel.png` — per-component locked vs rebuilt, with real-AMP and
  decoy reference lines
- `seqme_panel_locked.csv`, `seqme_panel_rebuild.csv` — full metric panels across
  all nine cohorts
- `no_regression_decision.csv` — the accept/reject table
- `loo_audit.csv`, `frontier.csv` — held-out judge audit and the
  potency/diversity frontier

## DOI citation

Opt in.
