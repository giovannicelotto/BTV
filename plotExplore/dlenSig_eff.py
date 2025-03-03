# %%
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import mplhep as hep
hep.style.use("CMS")
# %%
nFiles = 100
rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/noCut/flatGV"
#rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/TTToHadronic_withCutAtNano"
fileNames = glob.glob(rootDir+"/*.root")
if len(fileNames)<nFiles:
    nFiles = len(fileNames)
else:
    fileNames=fileNames[:nFiles]
print("Taking %d files"%len(fileNames))
for fileIdx, fileName in enumerate(fileNames):
    print("Opening %s"%fileName)
    
    f = uproot.open(fileName)
    try:
        tree=f['tree']
    except:
        print("Skipping file. Broken")
        os.remove(fileName)
        continue
    print("File %d"%fileIdx)
    if fileIdx == 0:
        SV_dlenSig = np.array(tree['SV_dlenSig'])
        distance = np.array(tree['distance'])
    else:
        distance        = np.concatenate((distance, np.array(tree['distance'])))
        SV_dlenSig        = np.concatenate((SV_dlenSig, np.array(tree['SV_dlenSig'])))

# %%
        
matched = np.where(distance>0, True, False)
thresholds = np.linspace(0, 20, 101)
roc = []
for t in thresholds:
    roc.append(np.sum(matched[SV_dlenSig>t])/len(matched))
roc = np.array(roc)





# %%
nFiles=100
rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/GVrow"
fileNames = glob.glob(rootDir+"/TTToH_*.root")
if len(fileNames)<nFiles:
    nFiles = len(fileNames)
else:
    fileNames=fileNames[:nFiles]
print("Taking %d files"%len(fileNames))
for fileIdx, fileName in enumerate(fileNames):
    print("Opening %s"%fileName)
    
    f = uproot.open(fileName)
    try:
        tree=f['tree']
    except:
        print("Skipping file. Broken")
        os.remove(fileName)
        continue
    print("File %d"%fileIdx)
    if fileIdx == 0:
        SV_dlenSig = np.array(tree['SV_dlenSig'])
        distance = np.array(tree['distance'])
    else:
        distance        = np.concatenate((distance, np.array(tree['distance'])))
        SV_dlenSig        = np.concatenate((SV_dlenSig, np.array(tree['SV_dlenSig'])))

# %%
        
matched = np.where(distance>0, True, False)
thresholds = np.linspace(0, 20, 101)
roc2 = []
for t in thresholds:
    roc2.append(np.sum(matched[SV_dlenSig>t])/len(matched))
roc2 = np.array(roc2)
# %%
fig, ax = plt.subplots(1, 1)
ax.plot(thresholds, roc, color='blue')
ax.plot(thresholds, roc2, color='red')
ax.set_xlabel("Cut on SV_dlenSig after matching")
ax.set_ylabel("Fraction of GV matched")
ax.set_ylim(0.1, 0.6)
ax.set_xlim(0., 20)
# %%
fig, ax = plt.subplots(1, 1)
ax.hist(SV_dlenSig, bins=np.linspace(0, 10, 30))

# %%
#
#   Now do the same for SV reco
#
# %%
nFiles = 100
rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/TTToHadronic_cut3SV_SVrow"
fileNames = glob.glob(rootDir+"/*.root")
if len(fileNames)<nFiles:
    nFiles = len(fileNames)
else:
    fileNames=fileNames[:nFiles]
print("Taking %d files"%len(fileNames))
for fileIdx, fileName in enumerate(fileNames):
    print("Opening %s"%fileName)
    
    f = uproot.open(fileName)
    try:
        tree=f['tree']
    except:
        print("Skipping file. Broken")
        os.remove(fileName)
        continue
    print("File %d"%fileIdx)
    if fileIdx == 0:
        SV_dlenSig = np.array(tree['SV_dlenSig'])
        SV_matched = np.array(tree['SV_matched'])
    else:
        SV_matched        = np.concatenate((SV_matched, np.array(tree['SV_matched'])))
        SV_dlenSig        = np.concatenate((SV_dlenSig, np.array(tree['SV_dlenSig'])))

# %%
thresholds = np.linspace(0, 20, 101)
roc = []
for t in thresholds:
    roc.append(np.sum(SV_matched[SV_dlenSig>t])/np.sum(SV_dlenSig>t))
roc = np.array(roc)
# %%
fig, ax = plt.subplots(1, 1)
ax.plot(thresholds, roc)
ax.set_xlabel("Cut on SV_dlenSig after matching")
ax.set_ylabel("SV matched / SV")
# %%
fig, ax = plt.subplots(1, 1)
ax.plot(SV_matched, SV_dlenSig, marker='o', linestyle='none', markersize=1)
ax.set_ylim(0, 100)
# %%
fig,ax = plt.subplots(1, 1)
bins = np.linspace(0, 20, 101)
cM =ax.hist(SV_dlenSig[SV_matched==1], bins=bins, label='Matched')[0]
cNM =ax.hist(SV_dlenSig[SV_matched==0], bins=bins, bottom=cM, label='not Matched')[0]
ax.legend()

fig,ax = plt.subplots(1, 1)
ax.hist(bins[:-1], bins=bins, weights=cM/(cM+cNM), label='Matched')[0]
ax.hist(bins[:-1], bins=bins, weights=cNM/(cM+cNM), bottom=cM/(cM+cNM), label='not Matched')[0]
ax.legend()
# %%
