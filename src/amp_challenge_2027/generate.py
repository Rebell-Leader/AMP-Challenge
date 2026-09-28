"""Deterministic, oracle-selected AMP library generator — AMP Challenge 2027.

Pipeline
--------
1. An order-3 residue Markov model (trained on 40,047 real AMPs: the MarLys
   antibacterial subset + QMAP/DBAASP peptides with panel MIC <= 16 uM) draws a
   1,000,000-candidate pool under fixed seed 2027. Pure Python/NumPy, no GPU,
   no network.
2. Candidates are filtered to the competition constraints (8-50 aa, 20 canonical
   amino acids, unique, no exact match to data/antibacterial.fasta).
3. The 50,000-member library is *selected* from that pool on published AMP
   oracles, taken from the seqme third-party plugin registry
   (https://github.com/szczurek-lab/seqme-thirdparty):
       amPEPpy     - https://github.com/szczurek-lab/seqme-amPEPpy
       AMPredictor - https://github.com/szczurek-lab/seqme-ampredictor
   Selection is by rank consensus within length strata, so the library's length
   distribution matches the training corpus exactly. The third registry oracle,
   AMPlify (https://github.com/szczurek-lab/seqme-amplify), was deliberately
   held out of library selection and used only as an independent judge.
4. The ranked top-100 lists are selected from the library by a four-judge
   consensus (amPEPpy, AMPredictor, AMPlify, and an in-house ESM-2/XGBoost
   strain-MIC ensemble) under hard novelty, hemolysis and internal-diversity
   gates. The generate_therapeutic list is ranked instead on predicted
   selectivity (HC50 / MIC safety window) subject to MIC <= 16 uM.

Reproducibility
---------------
Oracle inference needs a GPU and three separate legacy environments, so the
per-candidate oracle scores are shipped as a checkpoint (checkpoint/selection.npz)
exactly as trained model weights would be. The pool itself is regenerated from
the seed on every run, so `uv run generate_broad_spectrum` reproduces both
library.fasta and top.fasta byte-for-byte on any machine, CPU-only, in minutes.
The code that produced the shipped scores is in src/score_pool.py, and
`--verify-selection` re-derives the selection from the shipped scores and
asserts it matches the emitted library.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MODEL = os.path.join(ROOT, "checkpoint", "markov_model.json")
SELECTION = os.path.join(ROOT, "checkpoint", "selection.npz")
REFERENCE = os.path.join(ROOT, "data", "antibacterial.fasta")

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
SEED = 2027
POOL_N = 1_000_000
LIBRARY_SIZE = 50_000
TOP_K = 100


# --------------------------------------------------------------------------- #
# 1. generator
# --------------------------------------------------------------------------- #
def _load_markov(path: str = MODEL):
    mk = json.load(open(path))
    tc = {
        k: (list(v.keys()), np.array(list(v.values()), float) / sum(v.values()))
        for k, v in mk["trans"].items()
    }
    return mk, tc


def _gen_pool(n: int, seed: int, mk, tc, verbose: bool = True) -> list[str]:
    """Sample n peptides from the order-3 Markov model.

    Target length is drawn from the training-corpus length distribution and the
    end-token is honored only at/after the target length, so output lengths
    match real AMPs (median 18) instead of collapsing short.
    """
    rng = np.random.default_rng(seed)
    order = mk["order"]
    start_ctx = mk["start_ctx"]
    start_p = np.array(mk["start_p"])
    targets = rng.choice(mk["lengths"], size=n)

    out: list[str] = []
    for j, L in enumerate(targets):
        tgt = max(8, int(L))
        s = start_ctx[rng.choice(len(start_ctx), p=start_p)]
        guard = 0
        while len(s) < tgt and guard < 400:
            guard += 1
            key = s[-order:]
            if key not in tc:
                break
            toks, probs = tc[key]
            nxt = toks[rng.choice(len(toks), p=probs)]
            if nxt == "$":
                if len(s) >= tgt:
                    break
                continue
            s += nxt
        if len(s) > 50:
            s = s[:50]
        out.append(s)
        if verbose and j and j % 200_000 == 0:
            print(f"    sampled {j:,}/{n:,}", flush=True)
    return out


def _reference_set() -> set[str]:
    if not os.path.exists(REFERENCE):
        return set()
    return {
        line.strip().upper()
        for line in open(REFERENCE)
        if line.strip() and not line.startswith(">")
    }


def build_pool(verbose: bool = True) -> list[str]:
    """Regenerate the compliant, unique, novel candidate pool from the seed."""
    mk, tc = _load_markov()
    reference = _reference_set()
    if verbose:
        print(f"[1/3] sampling {POOL_N:,} candidates (seed {SEED})", flush=True)
    raw = _gen_pool(POOL_N, SEED, mk, tc, verbose=verbose)

    seen: set[str] = set()
    pool: list[str] = []
    for s in raw:
        if (
            8 <= len(s) <= 50
            and set(s) <= AA_SET
            and s not in reference
            and s not in seen
        ):
            seen.add(s)
            pool.append(s)
    if verbose:
        print(f"[2/3] compliant / unique / novel pool: {len(pool):,}", flush=True)
    return pool


# --------------------------------------------------------------------------- #
# 2. selection
# --------------------------------------------------------------------------- #
def _average_rank_pct(x: np.ndarray, higher_better: bool = True) -> np.ndarray:
    """Percentile rank with ties averaged (matches pandas rank(pct=True))."""
    n = len(x)
    order = np.argsort(x, kind="stable")
    ranks = np.empty(n, dtype=np.float64)
    sorted_x = x[order]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_x[j + 1] == sorted_x[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    ranks /= n
    return ranks if higher_better else 1.0 - ranks


def _stratified_top(score, keep_per_len, length_arr):
    chosen = []
    for L in sorted(keep_per_len):
        cand = np.flatnonzero(length_arr == L)
        if cand.size == 0:
            continue
        k = min(int(keep_per_len[L]), cand.size)
        chosen.append(cand[np.argsort(score[cand], kind="stable")[::-1][:k]])
    return np.concatenate(chosen)


def derive_selection(pool: list[str], ck) -> np.ndarray:
    """Re-derive the library indices from the shipped oracle scores."""
    lengths = np.array([len(s) for s in pool])
    hist = dict(zip(ck["length_hist_len"].tolist(), ck["length_hist_cnt"].tolist()))
    stage_b = ck["stage_b_idx"]
    cons = 0.5 * _average_rank_pct(ck["ampeppy"][stage_b], True) + 0.5 * _average_rank_pct(
        ck["ampredictor"], False
    )
    local = _stratified_top(cons, hist, lengths[stage_b])
    return stage_b[local]


# --------------------------------------------------------------------------- #
# 3. emit
# --------------------------------------------------------------------------- #
def _write_fasta(sequences, path, prefix):
    with open(path, "w") as f:
        for i, s in enumerate(sequences, start=1):
            f.write(f">{prefix}{i}\n{s}\n")


def build(category: str = "generate_broad_spectrum",
          verify_selection: bool = False, verbose: bool = True):
    pool = build_pool(verbose=verbose)
    ck = np.load(SELECTION)
    if int(ck["pool_size"]) != len(pool):
        raise RuntimeError(
            f"pool size mismatch: regenerated {len(pool)}, checkpoint "
            f"{int(ck['pool_size'])} — generator and checkpoint disagree."
        )

    library = [pool[i] for i in ck["library_idx"]]
    # The four activity categories share the potency-ranked top-100; the
    # Optimal-Selectivity category uses a separate safety-window ranking.
    top_key = "therapeutic_idx" if category == "generate_therapeutic" else "top100_idx"
    top = [pool[i] for i in ck[top_key]]

    if verify_selection:
        derived = derive_selection(pool, ck)
        if not np.array_equal(np.sort(derived), np.sort(ck["library_idx"])):
            raise RuntimeError("re-derived selection does not match shipped library_idx")
        print("    selection re-derived from shipped oracle scores: MATCH", flush=True)

    if len(library) != LIBRARY_SIZE or len(set(library)) != LIBRARY_SIZE:
        raise RuntimeError("library is not 50,000 unique sequences")
    if len(top) != TOP_K or not set(top) <= set(library):
        raise RuntimeError("top-100 is not a unique subset of the library")
    return library, top


def _emit(category: str, outdir: str | None = None, **kw):
    outdir = outdir or category
    os.makedirs(outdir, exist_ok=True)
    library, top = build(category=category, **kw)
    _write_fasta(library, os.path.join(outdir, "library.fasta"), "seq")
    _write_fasta(top, os.path.join(outdir, "top.fasta"), "top")
    print(
        f"[3/3] wrote {len(library):,} library + {len(top)} top sequences to {outdir}/",
        flush=True,
    )
    return library


def _main(category: str):
    p = argparse.ArgumentParser(description=f"Generate the {category} submission.")
    p.add_argument("--outdir", default=category)
    p.add_argument("--verify-selection", action="store_true",
                   help="re-derive the selection from the shipped oracle scores")
    p.add_argument("--quiet", action="store_true")
    a = p.parse_args()
    _emit(category, a.outdir, verify_selection=a.verify_selection, verbose=not a.quiet)


def generate_broad_spectrum():
    _main("generate_broad_spectrum")


def generate_gram_pos():
    _main("generate_gram_pos")


def generate_gram_neg():
    _main("generate_gram_neg")


def generate_mdr():
    _main("generate_mdr")


def generate_therapeutic():
    _main("generate_therapeutic")


if __name__ == "__main__":
    _main("generate_broad_spectrum")
