"""
Faithful AMPredictor inference wrapper.
Citation (competition PDF ref [4], verified against PDF bibliography):
  Ruihan Dong, Rongrong Liu, Ziyu Liu, Yangang Liu, et al.
  "Exploring the repository of de novo-designed bifunctional antimicrobial
   peptides through deep learning." eLife 13:RP97330 (2025).
  Code: github.com/ruihan-dong/AMPredictor
Held-out MIC benchmark: uses ESM-1b node features + ESM-1b contact-map graph +
per-residue Morgan fingerprint, exactly as the published model was trained.
Output: predicted logMIC (lower = more potent). NOT used anywhere in our selection.
"""
import os, sys, json, numpy as np, torch
os.environ.setdefault("TORCH_HOME", "/home/obi/AMP_challenge/model_cache/torch")
AMPDIR = "/home/obi/AMP_challenge/AMPredictor"
sys.path.insert(0, AMPDIR)
import esm
from rdkit.Chem import AllChem as Chem
from torch_geometric.data import Data as GData
from torch_geometric.data import Batch
from AMPredictor import GNNPredictor

MAXLEN = 65
_esm = _alph = _bc = None
_fp_lookup = None

def _load_esm():
    global _esm, _alph, _bc
    if _esm is None:
        _esm, _alph = esm.pretrained.esm1b_t33_650M_UR50S()
        _esm.eval()
        _bc = _alph.get_batch_converter()
    return _esm, _alph, _bc

def _load_fp_lookup():
    global _fp_lookup
    if _fp_lookup is None:
        smiles = json.load(open(os.path.join(AMPDIR, "data/smiles_file.json")))
        lut = {}
        for k, v in smiles.items():
            mol = Chem.MolFromSmiles(v)
            lut[k] = np.array(Chem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048))
        lut[' '] = np.zeros(2048, dtype=int); lut['0'] = np.zeros(2048, dtype=int)
        _fp_lookup = lut
    return _fp_lookup

def _fp_feature(seq):
    lut = _load_fp_lookup()
    fp = np.asarray([lut[a] for a in seq])                 # (L,2048)
    pad = np.zeros((MAXLEN - len(seq), 2048))
    fp = np.concatenate((fp, pad), axis=0)                 # (65,2048)
    return np.mean(fp, axis=1)                             # (65,) mean-pooled

@torch.no_grad()
def _esm_feats(seqs, device):
    """Return per-seq (node_features[65,1280], edge_index[2,E]) using ESM-1b."""
    model, alph, bc = _load_esm(); model = model.to(device)
    out_feats = []
    B = 8
    for i in range(0, len(seqs), B):
        chunk = seqs[i:i+B]
        labels, strs, toks = bc([(str(j), s) for j, s in enumerate(chunk)])
        toks = toks.to(device)
        out = model(toks, repr_layers=[33], return_contacts=True)
        reps = out["representations"][33].cpu()
        cons = out["contacts"].cpu()
        for k, s in enumerate(chunk):
            L = len(s)
            emb = reps[k, 1:L+1]                            # (L,1280) strip BOS/EOS
            if L < MAXLEN:
                emb = torch.cat((emb, torch.zeros(MAXLEN - L, 1280)), 0)
            contact = np.zeros((MAXLEN, MAXLEN))
            cm = cons[k, :L, :L].numpy()
            contact[:L, :L] = cm
            contact += np.eye(MAXLEN)
            r, c = np.where(contact >= 0.5)
            edge = np.array(list(zip(r, c)))
            out_feats.append((emb.float(), edge))
    return out_feats

def score(seqs, device=None, batch=128):
    """seqs: list[str] -> np.array logMIC predictions (lower=more potent)."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model = GNNPredictor().to(device)
    model.load_state_dict(torch.load(os.path.join(AMPDIR, "models/model_GNNPredictor_.model"),
                                     map_location=device))
    model.eval()
    feats = _esm_feats(seqs, device)
    fps = [_fp_feature(s) for s in seqs]
    preds = []
    with torch.no_grad():
        for i in range(0, len(seqs), batch):
            js = list(range(i, min(i+batch, len(seqs))))
            data_pro = []
            for j in js:
                emb, edge = feats[j]
                d = GData(x=emb,
                          edge_index=torch.LongTensor(edge).transpose(1, 0),
                          y=torch.FloatTensor([0.0]))
                data_pro.append(d)
            bpro = Batch.from_data_list(data_pro).to(device)
            bfp = torch.Tensor(np.array([fps[j] for j in js])).to(device)
            out = model(bpro, bfp)
            preds.append(out.cpu().numpy().flatten())
    return np.concatenate(preds)

if __name__ == "__main__":
    # smoke test on 3 known AMPs
    test = ["GIGKFLHSAKKFGKAFVGEIMNS",  # magainin-2, potent
            "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK",  # cecropin-ish
            "AAAAAAAAAA"]  # inert poly-A, should be weak
    p = score(test)
    for s, v in zip(test, p):
        print(f"logMIC={v:.3f}  {s}")
