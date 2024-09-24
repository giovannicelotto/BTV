# %%
import uproot
import numpy as np

fileName = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/tight/pt0p8/TTToH_265.root"
file = uproot.open(fileName)
tree = file['tree']
branches = tree.arrays()
maxEntries  = tree.num_entries
# %%
nSV = branches['nSV']
event = branches['event']
distance = branches['distance']

# %%
matched = np.sum(distance>0)
newSV = []
newEV = []
ev = -1
for nsv_, ev_ in zip(nSV, event):
    if ev_ == ev:
        continue
    else:
        ev_ = ev
        newSV.append(nsv_)



# %%
