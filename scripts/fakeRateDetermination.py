# %%
import uproot
import numpy as np
import glob
import os
# %%
fileName = "/work/gcelotto/BTV/scripts/tuplizer/TTToH_100.root"
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
    if ev == ev_:
        continue
    else:
        ev = ev_
        newSV.append(nsv_)


# newSV is the list of nSV per event 
# only for those events where there is at least one gv
# the number of matching can be taken from distance>0.
# the ratio gives the matchedSV/totalSV ( a bit optimistic because for those events where there are no GV there might be SV so the denominator can be a bit smaller (order1%))
# %%
print("Fake rate : ", 1-matched/np.sum(newSV))

# %%
#/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/tight/pt0p8
# /pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/tight/pt0p8
rootDir="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/GVrow" # change here
#rootDit="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/TTToHadronic2024Dec13Tuple_noCutSV_dlenSig"
fileNames = glob.glob(rootDir+"/*.root")
nFiles = 100
matched = 0
newSV = []
if nFiles == -1:
    pass
elif len(fileNames)<nFiles:
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
        
        nSV = tree['nSV']
        event = tree['event']
        distance = tree['distance']
    else:
        
        distance        = np.concatenate((distance, np.array(tree['distance'])))
        nSV             = np.concatenate((nSV, tree['nSV']))
        event             = np.concatenate((event, tree['event']))


    ev = -1
    for nsv_, ev_ in zip(nSV, event):
        if ev == ev_:
            continue
        else:
            ev = ev_
            newSV.append(nsv_)

    matched = matched + np.sum(np.array(distance)>0)
        
    #matched = np.array(tree['matched'])

# %%
print("Fake rate : ", 1-matched/np.sum(newSV))
# %%
