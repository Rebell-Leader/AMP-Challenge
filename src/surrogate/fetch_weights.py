"""
Stage 1 weight fetcher — caches all model weights under
/home/obi/AMP_challenge/model_cache so we load from disk, never re-hit HF.

IMPORTANT (sandbox proxy): HF's Xet backend hangs through the SOCKS proxy.
We force the classic LFS route with HF_HUB_DISABLE_XET=1 and require
httpx[socks]/socksio. This gives ~16 MB/s instead of a 0-byte stall.
"""
import os
CACHE = "/home/obi/AMP_challenge/model_cache"
os.environ["HF_HOME"] = f"{CACHE}/hf"
os.environ["HF_HUB_DISABLE_XET"] = "1"
os.environ["TORCH_HOME"] = f"{CACHE}/torch"
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

from huggingface_hub import snapshot_download

# repo_id -> allow_patterns (None = whole repo)
REPOS = {
    "facebook/esm2_t33_650M_UR50D": ["*.json","*.txt","pytorch_model.bin","*.py"],
    "facebook/esm2_t12_35M_UR50D":  ["*.json","*.txt","pytorch_model.bin","*.py"],
}

def fetch_all():
    import time
    for repo, pats in REPOS.items():
        t = time.time()
        p = snapshot_download(repo, allow_patterns=pats)
        print(f"[ok] {repo}  ({time.time()-t:.1f}s) -> {p}", flush=True)

if __name__ == "__main__":
    fetch_all()
    print("ALL WEIGHTS CACHED", flush=True)
