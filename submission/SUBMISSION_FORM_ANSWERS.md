# AMP Challenge 2027 — completed submission form answers (v2)

> Replace `ClaudeAMP` / `[TEAM NAME]` with the registered Kaggle team name before upload.

## Team / basic details
- Team name: [TEAM NAME]
- Track: Computational track (single track, auto-selected)
- Writeup title (47/80): Scoring our own library with the real benchmark
- Subtitle (122/140): Reproducing seqme and its three published AMP oracles showed our 50k library sat at the decoy floor. Reselection fixed it.
- Card image: writeup_card_560x280.png (560x280)
- DOI citation: opt in

## Categories entered
- Broad spectrum
- Optimal activity, Gram-positive
- Optimal activity, Gram-negative
- Optimal activity, MDR/WHO-priority
- Optimal selectivity (separate safety-window ranking)

All five entry points share the same 50,000-sequence library; the four activity
categories share the potency-ranked top-100, and Optimal Selectivity uses a
safety-window ranking.

## Files
- ClaudeAMP_library_50k.fasta — 50,000 sequences
- ClaudeAMP_top100.fasta — ranked top-100 (activity categories)
- ClaudeAMP_therapeutic_top100.fasta — ranked top-100 (Optimal Selectivity)
- Repository: github.com/Rebell-Leader/AMP-Challenge (root uv project, five entry points)

## Abstract (1,498 / 1,500 characters)
See ABSTRACT_1500.txt

## Method summary (9,999 / 10,000 characters)
See METHOD_SUMMARY_10k.txt

## Project description
See KAGGLE_WRITEUP.md

## Key numbers for quick reference
- Generator: order-3 Markov, 40,047-peptide corpus, seed 2027, 1,000,000-candidate pool
- Library selection: amPEPpy + AMPredictor rank consensus, length-stratified; AMPlify held out
- Library: 8 of 11 seqme components improved vs the previously locked library, 2 unchanged,
  1 regressed (Diversity 0.855 -> 0.816, i.e. to the real-potent-AMP level of 0.820)
- Held-out AMPlify: 0.577 -> 0.907 (real potent AMPs 0.929)
- Top-100: amPEPpy 0.905, AMPlify 1.000, AMPredictor median MIC 1.11 uM,
  median HC50 44.6 uM, median safety window 45.7x, worst reference identity 0.743
- Selectivity list: median HC50 199.9 uM, median MIC 0.69 uM, median window 276x
- library.fasta sha256 8a1535f8e509fcb6635ef9be789f842fee9356c0fa61bd495b4fb9f4d0842793

## Pre-upload checklist
- [ ] Rename files and fill the team-name field with the registered team name
- [ ] Push the updated repository; confirm the public clone URL resolves
- [ ] Register on Kaggle with an institutional email (co-authorship eligibility)
- [ ] Paste title, subtitle, description; upload the card image
- [ ] Submit by October 1, 2026 AOE
