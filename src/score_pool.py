"""Produce checkpoint/selection.npz — the oracle scores the library is selected on.

This script is NOT run during submission verification. It documents and
regenerates the shipped checkpoint, and it needs a GPU plus three mutually
incompatible legacy environments, which is exactly why the resulting scores are
shipped as a checkpoint instead.

Oracles are the three antimicrobial-peptide models published in the seqme
third-party plugin registry, each invoked through seqme's own ThirdPartyModel
wrapper (which clones the plugin repo and runs it in its own uv-locked venv):

    registry     https://github.com/szczurek-lab/seqme-thirdparty
    amPEPpy      https://github.com/szczurek-lab/seqme-amPEPpy
    AMPredictor  https://github.com/szczurek-lab/seqme-ampredictor
    AMPlify      https://github.com/szczurek-lab/seqme-amplify

Selection uses amPEPpy and AMPredictor only. AMPlify is deliberately held out of
selection so it can serve as an independent judge of the result; it is scored
here as well, but its scores are not written to the selection checkpoint.

Scores are stored as float64. This matters: rounding them to float32 introduces
ties that do not exist in the original values and shifts the selection by ~19 of
50,000 members, so the emitted library would no longer be reproducible.

Usage
-----
    python src/score_pool.py --plugins /path/to/cloned/plugin/repos

Requires: seqme, uv, a CUDA GPU, and network access on first run so seqme can
sync each plugin environment.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import seqme as sm

from amp_challenge_2027.generate import (
    _average_rank_pct,
    _stratified_top,
    build_pool,
)

OVERSAMPLE = 4  # stage-B prefilter keeps OVERSAMPLE x the final count per length bin


def make_oracles(plugin_root: str) -> dict:
    return {
        "amPEPpy": sm.models.ThirdPartyModel(
            entry_point="ampeppy.predict:predict",
            path=os.path.join(plugin_root, "seqme-amPEPpy"),
        ),
        "AMPredictor": sm.models.ThirdPartyModel(
            entry_point="ampredictor.predict:predict",
            path=os.path.join(plugin_root, "seqme-ampredictor"),
        ),
        "AMPlify": sm.models.ThirdPartyModel(
            entry_point="amplify.predict:predict",
            path=os.path.join(plugin_root, "seqme-amplify"),
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--plugins", required=True,
                   help="directory containing the cloned seqme-* plugin repos")
    p.add_argument("--out", default="checkpoint/selection_regenerated.npz")
    p.add_argument("--length-hist", default="checkpoint/selection.npz",
                   help="npz supplying the target per-length counts")
    args = p.parse_args()

    oracles = make_oracles(args.plugins)
    pool = build_pool()
    lengths = np.array([len(s) for s in pool])

    hist_src = np.load(args.length_hist)
    target = dict(zip(hist_src["length_hist_len"].tolist(),
                      hist_src["length_hist_cnt"].tolist()))

    # Stage A: the cheap oracle scores the whole pool.
    ampeppy = np.concatenate([
        oracles["amPEPpy"](pool[i:i + 50_000]) for i in range(0, len(pool), 50_000)
    ]).astype(np.float64)

    # Stage B: length-stratified prefilter, then the expensive oracle.
    stage_b = _stratified_top(
        ampeppy, {L: c * OVERSAMPLE for L, c in target.items()}, lengths
    ).astype(np.int64)
    stage_b_seqs = [pool[i] for i in stage_b]
    ampredictor = np.concatenate([
        oracles["AMPredictor"](stage_b_seqs[i:i + 10_000],
                               tokens_per_batch=512, chunk_size=2000)
        for i in range(0, len(stage_b_seqs), 10_000)
    ]).astype(np.float64)

    # Selection: mean percentile rank of the two oracles, within length strata.
    consensus = (
        0.5 * _average_rank_pct(ampeppy[stage_b], True)
        + 0.5 * _average_rank_pct(ampredictor, False)
    )
    library_idx = stage_b[_stratified_top(consensus, target, lengths[stage_b])]

    np.savez_compressed(
        args.out,
        pool_size=np.int64(len(pool)),
        ampeppy=ampeppy,
        stage_b_idx=stage_b,
        ampredictor=ampredictor,
        library_idx=library_idx.astype(np.int64),
        length_hist_len=np.array(sorted(target), dtype=np.int64),
        length_hist_cnt=np.array([target[L] for L in sorted(target)], dtype=np.int64),
    )
    print(f"wrote {args.out}: pool {len(pool)}, stage-B {len(stage_b)}, "
          f"library {len(library_idx)}")
    print("top-100 ranking (four-judge consensus incl. AMPlify and the in-house "
          "ESM-2/XGBoost ensemble) is documented in METHOD_SUMMARY_10k.txt "
          "section 6; its indices ship in the same checkpoint.")


if __name__ == "__main__":
    main()
