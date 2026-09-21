import os, json, numpy as np, pickle, torch
os.environ.setdefault("HF_HOME","/home/obi/AMP_challenge/model_cache/hf")
os.environ.setdefault("HF_HUB_OFFLINE","1"); os.environ.setdefault("HF_HUB_DISABLE_XET","1")
BASE="/home/obi/AMP_challenge"; AA="ACDEFGHIKLMNPQRSTVWY"
from transformers import AutoTokenizer, AutoModel
tok=AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
esm=AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").eval(); torch.set_num_threads(16)
bundle=pickle.load(open(f"{BASE}/surrogate_models/ensemble.pkl","rb"))
keep=[sp for sp in bundle["calibration"] if not sp.startswith("_") and bundle["calibration"][sp]["keep"]]
_kd={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,
     'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def physchem(s):
    n=len(s); c=np.array([s.count(a) for a in AA],float)/n
    charge=s.count('K')+s.count('R')-s.count('D')-s.count('E')
    return np.concatenate([c,[n/50,charge/10,charge/n,np.mean([_kd[a] for a in s])/5,
                              sum(s.count(a) for a in "AILMFVWC")/n,(s.count('K')+s.count('R')+s.count('H'))/n]])
@torch.no_grad()
def embed(seqs,bs=64):
    out=[]
    for i in range(0,len(seqs),bs):
        enc=tok(seqs[i:i+bs],return_tensors="pt",padding=True)
        rep=esm(**enc).last_hidden_state
        m=enc["attention_mask"].clone().unsqueeze(-1).float(); m[:,0]=0
        for j,l in enumerate(enc["attention_mask"].sum(1)): m[j,l-1]=0
        out.append(((rep*m).sum(1)/m.sum(1).clamp(min=1)).cpu().numpy().astype(np.float32))
        if i%1280==0: print(f"  embed {i}/{len(seqs)}",flush=True)
    return np.concatenate(out)

cand=json.load(open(f"{BASE}/candidates.json"))
cand_idx=np.load(f"{BASE}/cand_idx.npy")
print(f"scoring {len(cand)} candidates",flush=True)
X=embed(cand); P=np.stack([physchem(s) for s in cand]); Xf=np.concatenate([X,P],1)
pp=bundle["potency_clf"].predict_proba(Xf)[:,1]; hc=10**bundle["hc50"].predict(Xf)
mic={sp:10**bundle["mic_heads"][sp].predict(Xf) for sp in keep}
mic_mat=np.stack([mic[sp] for sp in keep],1); frac_active=(mic_mat<=16).mean(1)
sel=pp+0.3*frac_active+0.1*np.clip(np.log10(hc)/np.log10(128),0,1)
np.savez(f"{BASE}/cand_scores.npz", idx=cand_idx, pp=pp, hc=hc, frac_active=frac_active,
         sel=sel, mic=mic_mat, keep=np.array(keep))
# top-100 with a light diversity cap: greedy, skip if >0.8 identity to a chosen one
import Levenshtein
order=np.argsort(sel)[::-1]; chosen=[]; chosen_seq=[]
for i in order:
    s=cand[i]
    if any(Levenshtein.ratio(s,c)>0.8 for c in chosen_seq): continue
    chosen.append(i); chosen_seq.append(s)
    if len(chosen)==100: break
with open(f"{BASE}/broad_spectrum_top.fasta","w") as f:
    for r,i in enumerate(chosen): f.write(f">top{r+1}\n{cand[i]}\n")
json.dump({"n_cand":len(cand),"potency_med":float(np.median(pp)),"frac_potent_gt0.5":float(np.mean(pp>0.5)),
           "top100_potency_med":float(np.median(pp[chosen])),
           "top100_frac_active_med":float(np.median(frac_active[chosen])),
           "top100_hc50_med":float(np.median(hc[chosen])),"keep_strains":keep},
          open(f"{BASE}/stage2_report.json","w"),indent=2)
print("DONE",flush=True)
