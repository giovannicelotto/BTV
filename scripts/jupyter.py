#%%
import uproot
import numpy as np
import matplotlib.pyplot as plt
from tuplizer.utilsForScript import getPdgMask
import glob
import os
import mplhep as hep
hep.style.use("CMS")
filePath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/tight/pt0p8"
fileNames = glob.glob(filePath + "/*.root", recursive=True)[:10]
def flatten(xss):
    return [x for xs in xss for x in xs]
# %%

disp = []
dist = []
event = []
for fileName in fileNames:
    print("Opening ", fileName)
    file = uproot.open(fileName)
    tree=file['tree']

    branches = tree.arrays()
    print(branches)
    #disp.append(branches["displacement"])
    #dist.append(branches["distance"])
    event.append(branches["event"])

# %%
print(disp)

# %%
disp=np.array(flatten(disp))
dist=np.array(flatten(dist))
event=np.array(flatten(event))

# %%
bins=np.linspace(0, .3, 31)

fig, ax = plt.subplots(1, 1)
matchedMask = dist>-1
cM = np.histogram(disp[matchedMask], bins=bins)[0]
cM=cM/np.sum(cM)

cNM = np.histogram(disp[~matchedMask], bins=bins)[0]
cNM=cNM/np.sum(cNM)

ax.hist(bins[:-1], bins=bins, weights=cM, alpha=0.4, label='Matched', linewidth=3)
ax.hist(bins[:-1], bins=bins, weights=cNM, alpha=0.4, label='Non Matched', linewidth=3)
ax.hist(bins[:-1], bins=bins, weights=cM, histtype=u'step', linewidth=3, color='C0')
ax.hist(bins[:-1], bins=bins, weights=cNM, histtype=u'step', linewidth=3, color='C1')
ax.legend()
hep.cms.label()
ax.set_xlabel("Distance of flight [cm]")
fig.savefig("/t3home/gcelotto/BTV/plots/displacementM_vs_NM.png", bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 1)
counts = ax.hist(event, bins=np.arange(500))[0]
np.mean(counts)
# %%
