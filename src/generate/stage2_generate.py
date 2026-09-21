"""
Stage 2 — compliant 50k library generation + surrogate scoring.

Route decision (documented). Of the two organizer baselines named in the
competition PDF (HydrAMP + AMP-Diffusion), only HydrAMP was locatable and
evaluated here; a public AMP-Diffusion weight release was not found (the
guessed repo 404'd) and was not pursued further. A third, unrelated model
(LAMP GRU-VAE, from a competition co-author) was also trialed as a
convenient modern-PyTorch generator. Findings:
  - HydrAMP (szczurek-lab/HydrAMP): TF 2.2.1 / Keras 2.3.1 / Py3.8 legacy
    stack + custom wheels + 466 MB Google-Drive checkpoints — conflicts with
    the torch/numpy-2 env.
  - AMP-Diffusion: no public weight release located; not evaluated.
  - LAMP GRU-VAE (pszmk/lamp-grugru-vae): NOT AMP-Diffusion — a separate
    co-author model. Loads cleanly but is a low-fidelity demo export (0/10
    exact reconstruction; non-autoregressive GRU decoder loses long-range
    coherence). Off-manifold under prior sampling (negative net charge vs
    AMPs' +2.8).

Adopted route (filter-first + zero-cost): an order-3 residue Markov model
trained on 40k real AMPs (MarLys antibacterial subset + QMAP potent peptides,
MIC<=16uM). On-manifold by construction (median length 18, mean charge +2.9),
100% compliant. The Stage-1 surrogate ensemble does the quality lifting.
This is a strong, honest first library; a fine-tuned generator is a later,
optional upgrade (Stage 4) if the ensemble-filtered pool proves insufficient.
"""
import os, json, numpy as np, pickle, torch
os.environ.setdefault("HF_HOME","/home/obi/AMP_challenge/model_cache/hf")
os.environ.setdefault("HF_HUB_OFFLINE","1")
os.environ.setdefault("HF_HUB_DISABLE_XET","1")

BASE="/home/obi/AMP_challenge"
AA="ACDEFGHIKLMNPQRSTVWY"; AAset=set(AA)

# ---------- generator ----------
def load_markov(path=f"{BASE}/markov_model.json"):
    mk=json.load(open(path))
    tc={k:(list(v.keys()), np.array(list(v.values()),float)/sum(v.values()))
        for k,v in mk["trans"].items()}
    return mk, tc

def gen_pool(n, seed, mk, tc):
    """Sample n peptides. Target length is drawn from the training length
    distribution; the end-token is only honored AT/AFTER the target length,
    so the output length distribution matches the corpus (median ~18) instead
    of collapsing short (an early-$ bug shifted the median to 13)."""
    rng=np.random.default_rng(seed)
    order=mk["order"]; sc=mk["start_ctx"]; sp=np.array(mk["start_p"])
    lens=rng.choice(mk["lengths"], size=n)
    out=[]
    for L in lens:
        tgt=max(8,int(L))
        s=sc[rng.choice(len(sc),p=sp)]; g=0
        while len(s)<tgt and g<400:
            g+=1; key=s[-order:]
            if key not in tc: break
            toks,probs=tc[key]; nxt=toks[rng.choice(len(toks),p=probs)]
            if nxt=="$":
                if len(s)>=tgt: break   # only end at/after target length
                continue                 # early end: resample this step
            s+=nxt
        if len(s)>50: s=s[:50]
        out.append(s)
    return out

# ---------- scorer ----------
_kd={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
     'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def physchem(s):
    n=len(s); c=np.array([s.count(a) for a in AA],float)/n
    charge=s.count('K')+s.count('R')-s.count('D')-s.count('E')
    return np.concatenate([c,[n/50,charge/10,charge/n,np.mean([_kd[a] for a in s])/5,
                              sum(s.count(a) for a in "AILMFVWC")/n,
                              (s.count('K')+s.count('R')+s.count('H'))/n]])

def get_scorer():
    from transformers import AutoTokenizer, AutoModel
    tok=AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    esm=AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").eval()
    torch.set_num_threads(16)
    bundle=pickle.load(open(f"{BASE}/surrogate_models/ensemble.pkl","rb"))
    keep=[sp for sp in bundle["calibration"] if not sp.startswith("_")
          and bundle["calibration"][sp]["keep"]]
    @torch.no_grad()
    def embed(seqs,bs=64):
        out=[]
        for i in range(0,len(seqs),bs):
            enc=tok(seqs[i:i+bs],return_tensors="pt",padding=True)
            rep=esm(**enc).last_hidden_state
            m=enc["attention_mask"].clone().unsqueeze(-1).float(); m[:,0]=0
            for j,l in enumerate(enc["attention_mask"].sum(1)): m[j,l-1]=0
            out.append(((rep*m).sum(1)/m.sum(1).clamp(min=1)).cpu().numpy().astype(np.float32))
            if i%3200==0: print(f"    embed {i}/{len(seqs)}",flush=True)
        return np.concatenate(out)
    def score(seqs):
        X=embed(seqs); P=np.stack([physchem(s) for s in seqs]); Xf=np.concatenate([X,P],1)
        res={"p_potent":bundle["potency_clf"].predict_proba(Xf)[:,1],
             "hc50":10**bundle["hc50"].predict(Xf)}
        for sp in keep: res[f"mic_{sp}"]=10**bundle["mic_heads"][sp].predict(Xf)
        return res, keep
    return score

if __name__=="__main__":
    mk,tc=load_markov()
    marlys=set(l.strip().upper() for l in open(f"{BASE}/amp-challenge-2027/data/antibacterial.fasta")
               if l.strip() and not l.startswith(">"))
    # over-generate, filter to compliant + novel + unique, take 50k
    print("generating...",flush=True)
    raw=gen_pool(70000, seed=2027, mk=mk, tc=tc)
    seen=set(); lib=[]
    for s in raw:
        if 8<=len(s)<=50 and set(s)<=AAset and s not in marlys and s not in seen:
            seen.add(s); lib.append(s)
        if len(lib)>=50000: break
    print(f"compliant unique novel library: {len(lib)}",flush=True)
    assert len(lib)==50000, f"only {len(lib)} — increase over-generation"
    # score all 50k
    print("scoring 50k through surrogate ensemble...",flush=True)
    score=get_scorer()
    res,keep=score(lib)
    # composite selection score: potency prob + safety + mean strain activity
    pp=res["p_potent"]; hc=res["hc50"]
    mic_mat=np.stack([res[f"mic_{sp}"] for sp in keep],1)   # [N, n_strain]
    frac_active=(mic_mat<=16).mean(1)                       # coverage across reliable strains
    sel = pp + 0.3*frac_active + 0.1*np.clip(np.log10(hc)/np.log10(128),0,1)
    np.save(f"{BASE}/library_scores.npy", np.column_stack([pp,hc,frac_active,sel]))
    # write library.fasta + top.fasta (top-100 by sel), enforce diversity later
    order=np.argsort(sel)[::-1]
    with open(f"{BASE}/broad_spectrum_library.fasta","w") as f:
        for i,s in enumerate(lib): f.write(f">seq{i+1}\n{s}\n")
    top=order[:100]
    with open(f"{BASE}/broad_spectrum_top.fasta","w") as f:
        for r,i in enumerate(top): f.write(f">top{r+1}\n{lib[i]}\n")
    json.dump({"n":len(lib),"potency_med":float(np.median(pp)),
               "frac_potent_gt0.5":float(np.mean(pp>0.5)),
               "top100_potency_med":float(np.median(pp[top])),
               "top100_frac_active_med":float(np.median(frac_active[top])),
               "keep_strains":keep},
              open(f"{BASE}/stage2_report.json","w"), indent=2)
    print("DONE — wrote library, top, scores, report",flush=True)
