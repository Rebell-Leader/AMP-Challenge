"""
Faithful MBC-Attention inference wrapper (PDF ref [33]).
Citation (verified against competition PDF bibliography, ref [33]):
  Jielu Yan, Bob Zhang, Mingliang Zhou, Francois-Xavier Campbell-Valois, Shirley W.I. Siu.
  "A deep learning method for predicting the minimum inhibitory concentration of
   antimicrobial peptides against Escherichia coli using Multi-Branch-CNN and Attention."
   mSystems 8(4):e00345-23 (2023).  Code: github.com/jieluyan/MBC-Attention
NOTE: 'MBC' = Multi-Branch CNN (architecture), predicts E. coli MIC (uM).
Output: predicted E. coli MIC in micromolar (lower = more potent).
"""
import os, sys, tempfile, numpy as np, pandas as pd
MBCDIR="/home/obi/AMP_challenge/MBC-Attention"
os.chdir(MBCDIR); sys.path.insert(0, MBCDIR)
from tools.MultiBranchCNN import (geneFastasFromFastaFile, CNNimportFtsDataSetsPoNe,
                                  CNNstandardInputOutput, def_fts, def_scale, def_bias)
from tools.base import readFastaYan
from tensorflow import keras

MDL=os.path.join(MBCDIR,"model/whole_train.mdl")

def score(seqs):
    """seqs: list[str] -> np.array of predicted E. coli MIC (uM)."""
    tmp=tempfile.mkdtemp()
    fpath=os.path.join(tmp,"q.fasta")
    with open(fpath,"w") as f:
        for i,s in enumerate(seqs): f.write(f">q{i}\n{s}\n")
    _, fastas_file = geneFastasFromFastaFile(fpath)
    fastas = pd.read_csv(fastas_file)          # re-read so SEQUENCE is str, not Bio.Seq (matches repo test flow)
    sets = CNNimportFtsDataSetsPoNe(fastas, ft_list=def_fts, target=None)
    X, Y = CNNstandardInputOutput(sets)
    mdl = keras.models.load_model(MDL)
    pred = mdl.predict(X, verbose=0)
    pred = pred/def_scale - def_bias          # de-scale to pMIC-space (same as repo)
    mic_uM = 10 ** (-pred)                     # linear uM (repo's own conversion)
    return np.asarray(mic_uM).flatten()

if __name__=="__main__":
    test=["GIGKFLHSAKKFGKAFVGEIMNS","KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK","AAAAAAAAAA"]
    for s,v in zip(test, score(test)): print(f"E.coli MIC={v:.2f} uM  {s}")
