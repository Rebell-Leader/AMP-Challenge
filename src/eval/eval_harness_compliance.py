"""
AMP Challenge 2027 — Stage 0 compliance gate.

Owns the exact submission contract from the starter kit's verify_submission.py,
plus the MarLys top-100 identity check in BOTH forms:
  (a) Levenshtein.ratio > 0.80   -> what the shipped verifier actually enforces
  (b) MMseqs2 alignment identity  -> what the competition PDF (§1.3) specifies

Both are provided because they are NOT equivalent; a top-100 list must satisfy
whichever the organizers run at scoring time. Treat >0.80 Levenshtein as the
hard local gate and MMseqs2 >=0.80 as an additional flag.

Usage:
    from eval_harness_compliance import check_library, check_top100, load_fasta
    ok, report = check_library("library.fasta")
    ok, report = check_top100("top.fasta", "library.fasta",
                              marlys="data/antibacterial.fasta")
"""
from __future__ import annotations
from pathlib import Path
import subprocess, tempfile, shutil, os

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 8, 50            # PDF + verifier: 8..50 inclusive
LIBRARY_SIZE, TOP_SIZE = 50_000, 100
SIM_THRESHOLD = 0.80               # top-100 vs MarLys


def load_fasta(path) -> tuple[list[str], list[str]]:
    headers, seqs, hdr, parts = [], [], None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if hdr is not None:
                headers.append(hdr); seqs.append("".join(parts))
            hdr, parts = line[1:], []
        else:
            parts.append(line.upper())
    if hdr is not None:
        headers.append(hdr); seqs.append("".join(parts))
    return headers, seqs


def check_library(lib_path, marlys=None) -> tuple[bool, dict]:
    """Full-library gate: count, alphabet, length, uniqueness, no exact MarLys match."""
    headers, seqs = load_fasta(lib_path)
    errs, warns = [], []
    if len(seqs) != LIBRARY_SIZE:
        errs.append(f"library size {len(seqs)} != {LIBRARY_SIZE}")
    seen = set()
    bad_alpha = bad_len = dups = 0
    for h, s in zip(headers, seqs):
        if not h.strip():
            errs.append("missing header")
        if set(s) - STANDARD_AA:
            bad_alpha += 1
        if not (MIN_LEN <= len(s) <= MAX_LEN):
            bad_len += 1
        if s in seen:
            dups += 1
        seen.add(s)
    if bad_alpha: errs.append(f"{bad_alpha} seqs with non-canonical residues")
    if bad_len:   errs.append(f"{bad_len} seqs outside {MIN_LEN}-{MAX_LEN} aa")
    if dups:      errs.append(f"{dups} duplicate sequences")
    exact_overlap = 0
    if marlys is not None:
        ref = set(load_fasta(marlys)[1])
        exact_overlap = len(seen & ref)
        if exact_overlap:
            errs.append(f"{exact_overlap} library seqs are exact MarLys matches")
    return (len(errs) == 0), {
        "n": len(seqs), "unique": len(seen), "bad_alphabet": bad_alpha,
        "bad_length": bad_len, "duplicates": dups,
        "exact_marlys_overlap": exact_overlap, "errors": errs, "warnings": warns,
    }


def _levenshtein_max_ratio(top_seqs, ref_seqs):
    """Return worst (seq, ref, ratio) among pairs exceeding threshold, else best-below."""
    import Levenshtein
    ref_by_len = {}
    for r in ref_seqs:
        ref_by_len.setdefault(len(r), []).append(r)
    lengths = sorted(ref_by_len)
    worst = (None, None, 0.0)
    violations = []
    for s in top_seqs:
        ls = len(s)
        # ratio = 2*M/(|s|+|r|); to exceed T need |r| in [ (1-T)/(1+T)... ] window.
        # bound: ratio<=2*min(ls,lr)/(ls+lr); require > T -> lr in [ls*(1-T)/(1+T), ls*(1+T)/(1-T)]
        lo = ls * (1 - SIM_THRESHOLD) / (1 + SIM_THRESHOLD)
        hi = ls * (1 + SIM_THRESHOLD) / (1 - SIM_THRESHOLD)
        best = 0.0; best_ref = None
        for lr in lengths:
            if lr < lo or lr > hi:
                continue
            for r in ref_by_len[lr]:
                ratio = Levenshtein.ratio(s, r)
                if ratio > best:
                    best, best_ref = ratio, r
                    if best >= 1.0:
                        break
        if best > best if False else best > worst[2]:
            worst = (s, best_ref, best)
        if best > SIM_THRESHOLD:
            violations.append((s, best_ref, round(best, 3)))
    return worst, violations


def check_top100(top_path, lib_path, marlys=None, run_mmseqs=True) -> tuple[bool, dict]:
    """Top-100 gate: count, subset-of-library, dedup, MarLys identity (Levenshtein + optional MMseqs2)."""
    _, top = load_fasta(top_path)
    lib = set(load_fasta(lib_path)[1])
    errs, warns = [], []
    if len(top) != TOP_SIZE:
        errs.append(f"top size {len(top)} != {TOP_SIZE}")
    seen = set(); not_in = 0; dups = 0
    for s in top:
        if s not in lib: not_in += 1
        if s in seen: dups += 1
        seen.add(s)
    if not_in: errs.append(f"{not_in} top seqs not in library")
    if dups:   errs.append(f"{dups} duplicate top seqs")

    lev_violations = mm_violations = None; lev_worst = None
    if marlys is not None:
        ref = load_fasta(marlys)[1]
        lev_worst, lev_violations = _levenshtein_max_ratio(list(seen), ref)
        if lev_violations:
            errs.append(f"{len(lev_violations)} top seqs exceed Levenshtein 0.80 vs MarLys")
        if run_mmseqs and shutil.which("mmseqs"):
            mm_violations = _mmseqs_identity(list(seen), ref, marlys_path=marlys)
            if mm_violations:
                warns.append(f"{len(mm_violations)} top seqs >=80% MMseqs2 identity vs MarLys (PDF rule)")
    return (len(errs) == 0), {
        "n": len(top), "not_in_library": not_in, "duplicates": dups,
        "levenshtein_worst": (lev_worst[0], lev_worst[2]) if lev_worst else None,
        "levenshtein_violations": lev_violations,
        "mmseqs_violations": mm_violations,
        "errors": errs, "warnings": warns,
    }


def _mmseqs_identity(top_seqs, ref_seqs, marlys_path=None):
    """Flag top seqs with >=80% MMseqs2 pairwise identity to any MarLys entry."""
    tmp = Path(tempfile.mkdtemp(prefix="mmseqs_"))
    try:
        q = tmp/"query.fasta"; r = tmp/"ref.fasta"
        q.write_text("".join(f">t{i}\n{s}\n" for i, s in enumerate(top_seqs)))
        r.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(ref_seqs)))
        res = tmp/"res.m8"
        cmd = ["mmseqs","easy-search",str(q),str(r),str(res),str(tmp/"tmp"),
               "--min-seq-id","0.8","-s","7.5","--max-seqs","5",
               "--format-output","query,target,fident","-v","1"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            return {"mmseqs_error": p.stderr[-300:]}
        hits = {}
        if res.exists():
            for line in res.read_text().splitlines():
                qn, tn, fid = line.split("\t")[:3]
                fid = float(fid)
                if fid >= 0.8:
                    i = int(qn[1:])
                    hits[top_seqs[i]] = max(hits.get(top_seqs[i], 0), fid)
        return [(s, round(v, 3)) for s, v in hits.items()]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    import sys, json
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    marlys = "/home/obi/AMP_challenge/amp-challenge-2027/data/antibacterial.fasta"
    ok1, r1 = check_library(base/"library.fasta", marlys=marlys)
    ok2, r2 = check_top100(base/"top.fasta", base/"library.fasta", marlys=marlys)
    print(json.dumps({"library_ok": ok1, "library": r1,
                      "top100_ok": ok2, "top100": r2}, indent=2))
