"""Deterministic AMP library generator for the AMP Challenge 2027.

One order-3 residue Markov model (trained on 40,047 real AMPs: the MarLys
antibacterial subset + QMAP peptides with panel MIC <= 16 uM) generates a
50,000-sequence library. Length targets are drawn from the training-corpus
length distribution and the end-token is only honored at/after the target,
so the output length distribution matches real AMPs (median 18, mean charge
+2.9). All five category entry points share this library; each supplies its
own surrogate-ranked top-100 (shipped as data/<category>_top.fasta).

Reproducibility: fixed seed 2027, pure-Python/NumPy, no network, no GPU.
Re-running reproduces library.fasta byte-for-byte.
"""
import json, os, re, sys
import numpy as np

HERE = os.path.dirname(__file__)
PKG = os.path.dirname(os.path.dirname(HERE))
MODEL = os.path.join(PKG, "checkpoint", "markov_model.json")
DATA = os.path.join(PKG, "data")
AAset = set("ACDEFGHIKLMNPQRSTVWY")
SEED = 2027
LIBRARY_SIZE = 50000

# Exclude sequences whose substrings match common cloud-credential regexes
# (e.g. AWS access-key IDs). These are perfectly valid peptides, but they trip
# automated secret scanners in CI/hosting pipelines and get redacted in transit,
# which would corrupt the shipped FASTA. Filtering them keeps the submission
# machine-safe with negligible effect on library quality.
_SECRET_RE = re.compile(r'(?:AKIA|AGPA|AIDA|AROA|AIPA|ANPA|ANVA|ASIA)[ACDEFGHIKLMNPQRSTVWY]{16}')
def _trips_secret(s): return bool(_SECRET_RE.search(s))

def _load():
    mk = json.load(open(MODEL))
    tc = {k: (list(v.keys()), np.array(list(v.values()), float) / sum(v.values()))
          for k, v in mk["trans"].items()}
    return mk, tc

def _gen_pool(n, seed, mk, tc):
    rng = np.random.default_rng(seed)
    order = mk["order"]; sc = mk["start_ctx"]; sp = np.array(mk["start_p"])
    lengths = np.array(mk["lengths"])
    tgts = rng.choice(lengths, size=n)
    out = []
    for L in tgts:
        tgt = max(8, int(L))
        s = sc[rng.choice(len(sc), p=sp)]; g = 0
        while len(s) < tgt and g < 400:
            g += 1; key = s[-order:]
            if key not in tc: break
            toks, probs = tc[key]; nxt = toks[rng.choice(len(toks), p=probs)]
            if nxt == "$":
                if len(s) >= tgt: break
                continue
            s += nxt
        if len(s) > 50: s = s[:50]
        out.append(s)
    return out

def _marlys():
    fa = os.path.join(DATA, "antibacterial.fasta")
    if not os.path.exists(fa): return set()
    return set(l.strip().upper() for l in open(fa) if l.strip() and not l.startswith(">"))

def build_library():
    mk, tc = _load()
    marlys = _marlys()
    raw = _gen_pool(75000, SEED, mk, tc)
    seen = set(); lib = []
    for s in raw:
        if (8 <= len(s) <= 50 and set(s) <= AAset and s not in marlys
                and s not in seen and not _trips_secret(s)):
            seen.add(s); lib.append(s)
        if len(lib) >= LIBRARY_SIZE: break
    assert len(lib) == LIBRARY_SIZE, f"only {len(lib)}"
    return lib

def _emit(category, outdir):
    os.makedirs(outdir, exist_ok=True)
    lib = build_library()
    with open(os.path.join(outdir, "library.fasta"), "w") as f:
        for i, s in enumerate(lib): f.write(f">seq{i+1}\n{s}\n")
    top_src = os.path.join(DATA, f"{category}_top.fasta")
    if os.path.exists(top_src):
        top = [l.strip() for l in open(top_src) if l.strip() and not l.startswith(">")]
        with open(os.path.join(outdir, "top.fasta"), "w") as f:
            for i, s in enumerate(top): f.write(f">top{i+1}\n{s}\n")
    return lib

def generate_broad_spectrum(outdir="generate_broad_spectrum"): return _emit("broad_spectrum", outdir)
def generate_gram_pos(outdir="generate_gram_pos"):           return _emit("gram_pos", outdir)
def generate_gram_neg(outdir="generate_gram_neg"):           return _emit("gram_neg", outdir)
def generate_mdr(outdir="generate_mdr"):                     return _emit("mdr", outdir)
def generate_therapeutic(outdir="generate_therapeutic"):     return _emit("therapeutic", outdir)

if __name__ == "__main__":
    cat = sys.argv[1] if len(sys.argv) > 1 else "broad_spectrum"
    lib = _emit(cat, f"generate_{cat}")
    print(f"generated {len(lib)} sequences for {cat}")
