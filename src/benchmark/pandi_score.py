"""
Faithful Pandi DeepAMP regressor wrapper (competition PDF ref [18]).
Citation (verified against competition PDF bibliography, ref [18]):
  Amir Pandi, David Adam, Amir Zare, Van Tuan Trinh, Stefan L. Schaefer, et al.
  "Cell-free biosynthesis combined with deep learning accelerates de novo-development
   of antimicrobial peptides." Nature Communications 14(1):7197 (2023).
   Code: github.com/amirpandi/Deep_AMP
Ships Gram-split MIC regressors: CNN_gr_neg, LSTM_gr_neg, CNN_gr_pos, LSTM_gr_pos.
Output: predicted MIC (model units; lower = more potent). One-hot encoding, chop<=48.
"""
import os, sys, importlib.util, numpy as np
PANDIR="/home/obi/AMP_challenge/Pandi_DeepAMP"
os.chdir(PANDIR)
# 'code' collides with the stdlib module depending on cwd/sys.path order — load utils.py by path
_spec=importlib.util.spec_from_file_location("pandi_utils", os.path.join(PANDIR,"code","utils.py"))
_pu=importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_pu)
prepare_CNN=_pu.prepare_CNN
from tensorflow import keras

MODELS={
 "cnn_gr_neg":"saved_models/paper/regressor/CNN_gr_neg",
 "lstm_gr_neg":"saved_models/paper/regressor/LSTM_gr_neg",
 "cnn_gr_pos":"saved_models/paper/regressor/CNN_gr_pos",
 "lstm_gr_pos":"saved_models/paper/regressor/LSTM_gr_pos",
}
_cache={}
def _get(name):
    if name not in _cache: _cache[name]=keras.models.load_model(MODELS[name])
    return _cache[name]

def score(seqs, model="cnn_gr_neg"):
    """seqs: list[str] -> np.array MIC pred. Sequences >48aa are chopped out by the repo;
    we pre-truncate to 48 to keep alignment with input order."""
    trunc=[s[:48] for s in seqs]       # keep all seqs, truncate the few >48 (max is 50)
    X=prepare_CNN(trunc)
    mdl=_get(model)
    pred=mdl.predict(X, verbose=0).flatten()
    return pred

if __name__=="__main__":
    test=["GIGKFLHSAKKFGKAFVGEIMNS","KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK","AAAAAAAAAA"]
    for m in ["cnn_gr_neg","cnn_gr_pos"]:
        print(f"--- {m} ---")
        for s,v in zip(test, score(test,m)): print(f"  MIC_pred={v:.3f}  {s}")
