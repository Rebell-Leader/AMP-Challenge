# AMP Challenge 2027 — Completed Submission Form Answers

> Team name placeholder: **ClaudeAMP** — rename the two FASTA files and this field
> to your actual team name before uploading (files must be
> `TeamName_library_50k.fasta` and `TeamName_top100.fasta`).

---

**Please provide an abstract summarizing your method (maximum 1500 characters)** *
→ See `ABSTRACT_1500.txt` (1,341 characters). Paste its contents.

**Is your generative AI model or method previously published?** *
→ **No.**
(The Markov generator and the selection/validation pipeline were built for this
challenge. We *use* previously published models — ESM-2, and the named surrogates
AMPredictor/MBC-Attention/DeepAMP — but only for scoring/validation, not as our
generative method. If the form intends "does your generative method reuse any
published component," the honest note is: the generator itself is novel; the
scorers are published tools used off-the-shelf.)

**If previously published, provide the publication link.**
→ N/A (not previously published).

**In which competition categories are you primarily interested in having your
model assessed? (Multiple selections allowed)** *
→ ☑ Broad Spectrum Activity
→ ☑ Gram-Positive Activity
→ ☑ Gram-Negative Activity
→ ☑ Multi-Drug Resistant (MDR) Activity
→ ☐ Optimal Selectivity  (not selected — our single top-100 is optimized for
   broad activity; a selectivity-tuned list would trade away the activity
   categories where we are strongest)

**Short summary describing training data, external databases used, and any manual
intervention or computational filters applied (maximum 10,000 characters)** *
→ See `METHOD_SUMMARY_10k.txt` (7,513 characters). Paste its contents.

**Upload generated peptide library (50,000 sequences), FASTA, named
TeamName_library_50k.fasta** *
→ `ClaudeAMP_library_50k.fasta` (50,000 unique sequences, 8–50 aa, 20 canonical
   AA, no exact MarLys match, verified by organizers' verify_submission.py).
   Rename to your team name.

**Upload candidate list (top 100 sequences), FASTA, named TeamName_top100.fasta** *
→ `ClaudeAMP_top100.fasta` (100 unique sequences, all members of the library,
   worst Levenshtein identity 0.77 < 0.80 to MarLys, verified). Rename to your
   team name.

**All submitted sequences and experimental results will be made publicly
available under CC-BY 4.0 after the competition. Do you agree?** *
→ **Yes.** (Our sequences are generated from open data and we are committed to
   open release — this is consistent with co-authorship eligibility.)

**Do you intend to meet the Full Requirements for co-authorship (open-source code,
weights, and training data)?**
→ **Recommended: Yes.** Everything we used is already open: the generator is pure
   Python/NumPy with a fixed seed (byte-reproducible), the training corpus is
   MarLys-AMP (CC-0) + DBAASP-via-QMAP (open), and the model "weights" are a small
   JSON Markov table plus public ESM-2. There is no proprietary component to
   withhold, so meeting the full requirements costs us nothing and earns
   co-authorship. **Confirm with the team before selecting.**

---

## Pre-submission checklist (all verified)
- [x] Library: 50,000 unique sequences, 8–50 aa, 20 canonical AA
- [x] Library: no exact match to MarLys-AMP
- [x] Top-100: 100 unique sequences, all subset of the library
- [x] Top-100: worst Levenshtein identity 0.769 < 0.80 to MarLys
- [x] No secret-pattern-corrupted sequences
- [x] Passes organizers' own verify_submission.py end-to-end
- [x] Library reproducible byte-for-byte from fixed seed
- [x] Independently validated by 3 held-out judges (92–100% beat real-AMP median)
- [ ] Rename files + team-name field to your actual team name
- [ ] Register on Kaggle with an institutional email
- [ ] Submit by **October 1, 2026 AOE**
