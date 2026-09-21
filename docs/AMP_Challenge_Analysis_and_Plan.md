# AMP Challenge 2027 — Context, Doc Audit, and Actionable Plan

*Prepared from the authoritative competition PDF (`AMPChallenge_NeurIPSCompetition.pdf`, Szczurek et al., "under review"), the starter-kit repo (`amp-challenge-2027/`), and the three internal strategy docs. Every load-bearing claim below is traced to the PDF or verified against source; where the strategy docs disagree with the PDF, the PDF wins.*

---

## 1. Ground truth: what the competition actually is

The AMP Challenge is a **NeurIPS 2026 Competition-Track** benchmark run by the Szczurek lab (Helmholtz Munich / Warsaw) for generation + Phase-1 scoring, with wet-lab validation by the **de la Fuente lab (UPenn)**. It is the first standardized, wet-lab-falsified benchmark for generative AMP design. The organizers' own question: *which generative strategies, computational metrics, and selection rules predict experimental efficacy.*

### Hard facts (from the PDF — these override the strategy docs)

| Item | Authoritative value | Notes |
|---|---|---|
| **Platform** | Kaggle Community Hackathons | Institutional email required to register |
| **Deadline** | **October 1, 2026 AOE** | Docs assume "September 1" — **wrong** |
| **Phase 1 result** | Dec 2026 NeurIPS workshop | Qualification announced here |
| **Phase 2 assays** | Mar–Apr 2027 | Results TBA 2027; paper late 2027 |
| **Deliverable** | 50,000-seq library **+** ranked top-100, per category | One team = one library + one top-100 list |
| **Sequence rules** | 20 canonical AA; **8–50 aa**; unique; **linear, free termini**; no non-canonical/modified | Enforced by `verify_submission.py` |
| **Novelty (library)** | No sequence **exactly** matching MarLys-AMP | Repo `data/antibacterial.fasta` = this set |
| **Novelty (top-100)** | ≤ **80% identity** (MMseqs2) to any MarLys entry | Stricter; non-compliant auto-replaced by next valid |
| **Advance to Phase 2** | ≤ **20 teams**; 25 peptides each (≤500 total) | If >20 qualify, Phase-1 rank selects cohort |
| **Phase 2 sampling** | 25 drawn **uniformly at random** from your top-100 | You do not choose which 25 |
| **MIC ceiling** | 64 µM (>64 = inactive); **potency threshold MIC ≤ 16 µM** | Tighter than the 32 µM used in prior studies |
| **HC50 ceiling** | 128 µM | Safety window SW = HC50 / MIC50 |
| **Baselines** | **HydrAMP** + **AMP-Diffusion** (excluded from ranking) | Both from organizers; weights public |
| **Eval library** | **seqme** (`arxiv 2511.04239`, first author = organizer Møller-Larsen) | Deterministic given fixed version |
| **Compute stance** | *"participation does not require substantial compute … single CPU eligible … fine-tuning not required"* | Wet-lab is the funded prize, not compute |
| **Cost to enter** | €0; synthesis + assays funded centrally | The prize IS the wet-lab validation |

### The 20-strain panel (PDF Appendix B) — decisive and under-used by the docs

- **15 Gram-negative**, **5 Gram-positive**. **8 are MDR** (5 GN + 3 GP: CRE *E. coli* ×2, MDR *A. baumannii*, *K. pneumoniae*, MDR *P. aeruginosa*; MRSA, VRE *E. faecalis*, VRE *E. faecium*).
- The panel is **3:1 Gram-negative-weighted**. Broad-spectrum success rate is therefore dominated by Gram-negative activity — the harder target (outer membrane / LPS). This asymmetry is a strategic lever none of the three docs discuss.

### Phase 1 scoring — what is actually measured (PDF §1.5)

Four metric families under seqme, **weights and reference-set composition withheld until Phase 1 closes** (explicitly to prevent optimizing against the ranking function):
1. **Surrogate activity** — aggregated AMP-classifier + MIC-predictor scores (named: AMPredictor, MBC-Attention, DeepAMP), ensembled to reduce single-model bias.
2. **Sequence-level** — uniqueness, internal diversity, novelty (alignment bit-scores), clustering coverage.
3. **Embedding distributional similarity** — Fréchet Biological Distance, MMD, precision/recall vs. two reference sets (known-AMP + generic-peptide), using **ESM2 and ESM-C** embeddings.
4. **Property distribution** — charge & amphiphilicity conformity to known AMPs + synthesizability rate.

**Key correction:** Phase 1 is **not** "diversity/novelty-first." Surrogate *activity* is a co-equal metric family, and distributional similarity rewards looking like real AMPs — i.e., **realism is scored alongside diversity**. A pure-diversity pool that ignores predicted activity is penalized.

---

## 2. Audit of the three internal strategy docs

The docs are well-researched and the literature is real (I verified OmegAMP `2504.17247`, AmpLyze `2507.08162`, ProSpero `2505.22494` on arXiv; QMAP, STAMP, BATTLE-AMP, GenPept, the withdrawn GFlowNet paper on bioRxiv). But there are material errors and misalignments with the PDF.

### 2.1 Factual errors to fix before acting

| # | Claim in docs | Reality | Severity |
|---|---|---|---|
| E1 | Deadline "September 1, 2026" | **October 1, 2026 AOE** | High — 1 month more runway |
| E2 | Phase 1 = "diversity, novelty, physicochemical" only | Adds a **surrogate-activity** family + embedding similarity (ESM2/ESM-C) | High — reshapes the objective |
| E3 | "Novelty threshold = ≥30% identity from known AMPs" | No such rule. Library = **no exact match**; top-100 = **≤80% identity** (MMseqs2) to MarLys | High — the docs over-constrain generation |
| E4 | Baselines to beat = OmegAMP | Baselines are **HydrAMP + AMP-Diffusion**; OmegAMP is *a competing method*, not the yardstick | Med — changes the target |
| E5 | Length windows "6–30", "7–35", "5–50" (inconsistent across docs) | Competition window is **8–50** (verifier rejects <8 or >50) | High — sequences outside 8–50 fail compliance |
| E6 | seqme "confirmed to be Phase 1 eval" with known metrics | seqme is confirmed, but **weights/reference sets are withheld** — you cannot fit to the exact ranking | Med — kills any "optimize the metric" plan |
| E7 | "25 sampled from top-100" framed as curation leverage | 25 are **uniform random** — all 100 must be good; no cherry-picking | Med — argues for a *uniformly* strong top-100 |
| E8 | Reference DB variously "dbAMP/DRAMP/DBAASP" for the identity filter | The competition's reference is specifically **MarLys-AMP (MLAMP)** — it's the repo FASTA | Med — filter against the *right* set |
| E9 | Repo `weights.csv` implies a real model | It is the literal value `0.6` — a toy K/P Bernoulli sampler. The repo is a **format template**, not a model | Low — but the docs never say this |

### 2.2 Strategic misalignments (defensible calls, but worth revisiting)

- **Over-indexing on RL/LoRA fine-tuning.** The PDF says fine-tuning is *not required* and single-CPU methods are eligible. The July doc partly corrects this (ProSpero, frozen-EvoDiff), which is the right direction. RL should be *optional*, not the spine.
- **Gram-selective categories dismissed as "sparse data."** Given the 15:5 GN:GP panel, **Gram-negative and Broad-spectrum are effectively the same objective**, and Gram-positive is a *small, winnable* sub-panel (only 5 strains, 3 of them the well-studied *S. aureus*/enterococci). The docs under-weight Gram-positive as an easier category to place well in.
- **Top-100 curation framed around choosing winners.** Because 25 are random, the real objective is **maximize the worst-case / mean of the top-100**, not build a Pareto front you then hand-pick from. This favors tight, uniformly-high-quality lists (with the ≤80% identity + diversity constraints).
- **Compute budgets are over-engineered.** The whole task is doable on the local box (RTX 4080, 32 CPU, 61 GB) plus free public model weights. The AWS/Vast.ai plans ($95→$29) are fine but the honest number is **≈$0** for a first strong submission.
- **Surrogate-activity family ignored as a design target.** Since Phase 1 explicitly ensembles AMP classifiers + MIC predictors, our *own* generation should be filtered by a similar ensemble so our library scores well on *their* surrogate family — even though exact models/weights are hidden, the *class* of model is named.

### 2.3 What the docs get right (keep)

- "Filter-first" thesis (BATTLE-AMP, OmegAMP): the classifier/filter stage is the quality bottleneck, not the generator. **Correct and central.**
- Use organizer-adjacent tooling (seqme, HydrAMP, OmegAMP, ProSpero) — both to benchmark and to anticipate evaluation.
- MIC regressors beat binary classifiers (BATTLE-AMP); composition carries most signal.
- HC50 prediction is unreliable (QMAP) → treat AmpLyze as a soft/directional filter, keep boundary sequences.
- Reproducibility hygiene from day 1 (uv, fixed seed, public weights) for co-authorship.

---

## 3. Actionable plan — low-cost, staged, local-first

Philosophy: **spend the first two weeks proving we can score a library the way the organizers will, then generate against that.** Do not build infrastructure the competition says we don't need. Every stage produces a checkpoint that is already a valid (if weak) submission, so we are never at zero.

### Stage 0 — Reproduce the evaluation harness (Week 1) — **local, $0**
- Clone/run the starter kit end-to-end (`uv run generate_broad_spectrum`, then `verify_submission.py`) so we own the exact compliance contract.
- Install **seqme** and reproduce Phase-1 metrics on the two **baseline libraries** (HydrAMP, AMP-Diffusion) once the organizers release them; until then, run seqme on our own trial pools.
- Build the **MarLys identity gate** (exact-match for library; MMseqs2 ≤80% for top-100) as a standalone, tested module. This is a hard gate we must pass.
- **Deliverable:** `eval_harness/` that takes any FASTA → seqme report + compliance verdict. *This is the single highest-leverage artifact.*

### Stage 1 — Assemble the surrogate ensemble (Weeks 1–2) — **local, $0**
- Stand up an activity/safety scoring stack mirroring the *class* of models Phase 1 uses: an AMP classifier (OmegAMP-style XGBoost on physicochemical + composition + PLM features), a MIC regressor (per-strain / species-aware — STAMP or an ESM2 head on DBAASP), and AmpLyze for HC50.
- **Calibrate every predictor on QMAP's homology-aware split** before trusting it (this is the July doc's best single recommendation). Drop any strain head with test PCC < 0.5.
- Curate training data from the PDF-endorsed sources: **DBAASP, APD6, dbAMP 3.0, AMPSphere, Peptipedia v2, MarLys** (all open). De-dup, log-transform MIC, geometric-mean duplicates.
- **Deliverable:** `score_library(fasta) →` per-seq {p_AMP, per-strain MIC, HC50, physchem} table + a QMAP calibration report.

### Stage 2 — Baseline generation (Weeks 2–3) — **local GPU, $0**
- Generate the 50,000-library from **public pretrained weights only**: HydrAMP (cVAE, analogue + unconstrained) and/or OmegAMP (diffusion, property + subset conditioning). No training required.
- Tier the pool for the scored diversity Phase 1 rewards: unconditional broad sampling + property-conditioned (charge/length matched to the MarLys distribution: median 18 aa, mean charge +2.8, but *not* over-fit — 23% of real AMPs are ≤0 charge) + analogue diversification.
- Pass everything through Stage-0 gates + Stage-1 ensemble; keep predicted-active, novel, diverse.
- **Deliverable:** first fully-compliant 50k library + seqme score beating both baselines. *This alone is a submittable entry.*

### Stage 3 — Top-100 selection, per category (Weeks 3–4) — **local, $0**
- Because 25/100 are **random**, optimize the **mean and floor** of the list, not a hand-picked frontier.
- Category routing off the panel structure:
  - **Broad-spectrum** (primary): coverage = fraction of 20 strains with predicted MIC ≤ 16 µM, GN-weighted to match the 15:5 panel.
  - **Gram-positive** (opportunistic, only 5 strains): a genuinely winnable, lower-variance category — worth a dedicated list.
  - **Gram-negative / MDR** (stretch): only if QMAP-calibrated MIC heads clear PCC ≥ 0.5 on the ESKAPE strains.
  - **Selectivity**: SW = HC50/MIC50 with AmpLyze as a *soft* ranker; keep HC50-boundary sequences.
- Enforce ≤80% MarLys identity (auto-replace) + cluster cap (≤5/cluster, MMseqs2 50%) so the random-25 draw never lands on redundant peptides.
- **Deliverable:** category-specific top-100 lists + selection documentation (required by rules).

### Stage 4 — (Optional) surrogate-guided optimization (Weeks 4–6) — **local or ≤$30 cloud**
- Only if Stage-3 lists look weak: add **ProSpero** (frozen EvoDiff + SMC, no fine-tune cost) or a small LoRA-RL loop to push the top-100 toward higher predicted activity/selectivity while holding novelty.
- Keep entropy/diversity regularization (GFlowNet lesson) to avoid mode collapse.

### Stage 5 — Reproducibility & submission (Week 6+) — **local, $0**
- Package as the required `uv` project with fixed seed; verify `uv sync` + entry point reproduces the library byte-for-byte (organizers check this).
- Public repo, permissive license, full data disclosure → co-authorship eligibility.
- Submit via Kaggle by **Oct 1, 2026 AOE**.

---

## 4. Cloud, open models, and tooling — what we actually need

### Compute
- **Default: local.** The RTX 4080 (16 GB) + 32 CPU + 61 GB RAM box runs HydrAMP/OmegAMP inference, XGBoost, MMseqs2, seqme, and QMAP comfortably. The organizers explicitly designed the task for single-CPU eligibility.
- **Cloud only if** we (a) fine-tune ESM-2 650M / run AmpLyze's dual 3B-param embedding stack at scale, or (b) do a large RL loop. Then a **single A40/A100 spot instance for < 24 h (≈$5–30)** covers it. No compute provider is currently configured in this workspace — we'd add one (Customize → Compute) only at Stage 4, and only if needed.
- **Honest budget for a strong first submission: ~$0.**

### Open models / weights to pull (all public)
- **Generators:** HydrAMP (cVAE weights), AMP-Diffusion (baseline), OmegAMP (diffusion + XGBoost filter), EvoDiff (ProSpero backbone).
- **Scorers:** AmpLyze (HC50), STAMP / LLAMP (species MIC), an ESM-2 classifier head; ESM2 + ESM-C embeddings (used by the Phase-1 FBD/MMD family — worth mirroring).
- **Eval:** seqme, QMAP (`pip install qmap-benchmark`), MMseqs2.
- **Data:** DBAASP, APD6, dbAMP 3.0, AMPSphere, Peptipedia v2, MarLys (repo already ships the antibacterial MarLys subset).

### Do my current MCP tools suffice for the investigation?
**Yes for the research/analysis phase; partially for data pulls; no live gap that blocks starting.**

- **Literature & method scouting — sufficient.** `literature` (arXiv verified working; OpenAlex needs a free key — add under Customize → Credentials if we want citation-graph queries), `biorxiv` (verified — pulled QMAP/STAMP/BATTLE-AMP/GenPept), `pubmed`, `fetch_article_fulltext` cover the whole SOTA-tracking need.
- **Protein/structure context — sufficient.** `structures-interactions` (PDB/AlphaFold), `protein-annotation` (InterPro/STRING), `uniprot` via `genes-ontologies` for any mechanistic follow-up.
- **Chemistry — available if needed** (`chembl`, `chemistry`, `zinc`) though peptides are mostly out of scope for small-molecule tools.
- **Gaps (not blockers):** No MCP tool pulls DBAASP/dbAMP/APD/MarLys AMP records directly — those come via their web downloads / the repo FASTA / `qmap-benchmark`, all reachable from the analysis kernel. No MCP wraps seqme or the surrogate models — we run those as normal Python packages. Neither is a wall; both are "download + pip" steps.
- **Recommendation:** add an **OpenAlex API key** (free) for citation-graph literature sweeps, and note that everything else runs in-kernel. No new connector is required to begin.

---

## 5. Immediate next actions (this week)

1. Run the starter kit + `verify_submission.py` locally; confirm the exact compliance contract. *(Stage 0)*
2. Build + unit-test the **MarLys identity gate** (exact + MMseqs2 ≤80%). *(Stage 0)*
3. `pip install qmap-benchmark seqme`; reproduce QMAP splits and a seqme run on a trial pool. *(Stages 0–1)*
4. Pull DBAASP + dbAMP + MarLys; build the de-duplicated, log-MIC training table. *(Stage 1)*
5. Stand up one AMP classifier + one MIC head + AmpLyze, **calibrate on QMAP**, drop anything < 0.5 PCC. *(Stage 1)*
6. Generate a first compliant 50k library from HydrAMP/OmegAMP public weights and score it. *(Stage 2)*

*Correct the three docs' deadline (Oct 1), length window (8–50), novelty rule (exact / ≤80% MarLys), and Phase-1 objective (activity + realism + diversity, weights hidden) before anyone executes against them.*
