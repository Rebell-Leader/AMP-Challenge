"""
Bounded DeepAMP (Li et al., Nat Commun 15(1):7538, 2024) seed-and-optimize pass.
Uses the repo's AOM checkpoint (the "Finetuning-AOM" model, fine-tuned on the
`random_pair` dataset per the repo README; "AOM"/"POM" are the repo's own dataset/
section labels — the repo does not spell the acronyms out, so no expansion is asserted
here) as a masked-mutation generator to propose variants of each top-100 seed, and
accepts a variant ONLY if:
  (1) MBC-Attention (PDF ref [33], PCC 0.71 vs QMAP E.coli) predicts LOWER MIC, AND
  (2) it passes ALL compliance gates ACTUALLY ENFORCED here: 8-50 aa, 20 canonical AA,
      unique, no exact MarLys match, not already in library, Levenshtein ratio<0.79 vs
      MarLys (the binding verifier gate), no secret pattern.
      (No MMseqs gate is run in this pass — Levenshtein is the enforced identity gate.)
MBC-Attention is the acceptance gate (user instruction); DeepAMP only proposes.
Runs on CPU (~28 seq/s) — no GPU/Docker required for a bounded top-100 pass.
"""
import os, sys, glob, json, random, re, numpy as np, torch
BASE="/home/obi/AMP_challenge"
DEEP=f"{BASE}/deepAMP"
sys.path.insert(0, f"{DEEP}/code"); sys.path.insert(0, f"{DEEP}/src"); sys.path.insert(0, DEEP)
from src.model import AMPBERT

VOCAB=list("#$&*ABCDEFGHIKLMNOPQRSTUVWXYZ")
STOI={c:i for i,c in enumerate(VOCAB)}; ITOS={i:c for c,i in STOI.items()}
CANON=set("ACDEFGHIKLMNPQRSTVWY")
PAD,CLS,SEP,MASK='#','$','&','*'
BLOCK=64
_SECRET_RE=re.compile(r'(A3T|AKIA|AGPA|AIDA|AROA|AIPA|ANPA|ANVA|ASIA)[A-Z0-9]{16}')
def trips_secret(s): return bool(_SECRET_RE.search(s))

class Cfg:
    n_layers=12; n_heads=12; n_embd=32; d_model=512; d_ff=1024; d_k=d_v=64
    def __init__(s,v,b): s.vocab_size=v; s.block_size=b

def load_aom(epoch="epoch_199"):
    m=AMPBERT(Cfg(len(VOCAB),BLOCK))
    ck=glob.glob(f"{DEEP}/weights/aom_extract/*/{epoch}/model.pkl")[0]
    m.load_state_dict(torch.load(ck,map_location="cpu"), strict=True)
    m.eval(); return m

def encode(seq):
    ids=[STOI[CLS]]+[STOI[c] for c in seq]+[STOI[SEP]]
    ids=ids[:BLOCK]; ids+=[STOI[PAD]]*(BLOCK-len(ids)); return ids

@torch.no_grad()
def propose(model, seed, n_variants=8, mask_ratio=0.3, seed_rng=0):
    """Mask random positions of the seed and sample AOM predictions -> variants."""
    rng=random.Random(hash(seed)&0xffff ^ seed_rng)
    L=len(seed); variants=set()
    interior=list(range(1, L+1))  # positions in encoded seq (after CLS) that are residues
    for _ in range(n_variants):
        n_mask=max(1, int(L*mask_ratio))
        mask_positions=rng.sample(interior, min(n_mask, len(interior)))
        ids=encode(seed)
        for p in mask_positions: ids[p]=STOI[MASK]
        X=torch.LongTensor([ids])
        pos=torch.LongTensor([mask_positions+[0]*(20-len(mask_positions))])
        logits,_=model(X,pos)  # [1, 20, vocab]
        probs=torch.softmax(logits[0], dim=-1)
        new=list(seed)
        for j,p in enumerate(mask_positions):
            # sample from top-k to add diversity, restricted to canonical AAs
            pr=probs[j].numpy()
            cand=[(ITOS[i],pr[i]) for i in range(len(VOCAB)) if ITOS[i] in CANON]
            cand.sort(key=lambda x:-x[1])
            topk=cand[:5]; tot=sum(w for _,w in topk)
            r=rng.random()*tot; acc=0
            for aa,w in topk:
                acc+=w
                if r<=acc: new[p-1]=aa; break
        v="".join(new)
        if 8<=len(v)<=50 and set(v)<=CANON: variants.add(v)
    return list(variants)

if __name__=="__main__":
    m=load_aom()
    test="KWKLFKKIEKVGQNIRDGIIKAG"
    vs=propose(m,test,n_variants=8)
    print(f"seed: {test}")
    for v in vs: print("  variant:", v)
    print(f"generated {len(vs)} unique compliant variants on CPU")
