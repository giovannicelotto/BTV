#%%
import uproot
import numpy as np
import matplotlib.pyplot as plt
# %%
filePath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/240506_153739/0000/TTToHadronic_Run2_mc_2024May06_398.root"
file = uproot.open(filePath)
#%%
tree = file['Events']
branches = tree.arrays()
# %%
bins=np.arange(-1, 10)
cTot = np.zeros(len(bins)-1)

numTot=0
denTot=0
for ev in range(10000):#tree.num_entries:
    ProbeTracks_matchedToSV = branches["ProbeTracks_matchedToSV"][ev]
    nSV = branches["nSV"][ev]
    numEv = len(set(ProbeTracks_matchedToSV[ProbeTracks_matchedToSV>-1]))
    #print(ProbeTracks_matchedToSV[ProbeTracks_matchedToSV>-1])
    #print(numEv, nSV)
    denEv = nSV
    numTot = numTot + numEv
    denTot = denTot + denEv
print(numTot, denTot, numTot/denTot )

