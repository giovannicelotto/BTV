# %%
import uproot 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import awkward as ak
import mplhep as hep
hep.style.use("CMS")

def map_to_groups_letter(value):
    if abs(value) in [511, 521, 531, 541]:
        return 'B'
    elif abs(value) in [411, 421, 431]:
        return 'D'
    elif ((abs(value) > 3000) & (abs(value) < 4000)):
        return 'SB'
    elif ((abs(value) > 4000) & (abs(value) < 5000)):
        return 'CB'
    elif ((abs(value) > 5000) & (abs(value) < 6000)):
        return 'BB'
    # Add more conditions as needed
    else:
        return -1  # or any default value for unmatched cases
    
path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul10Tuple"
fileNames = glob.glob(path+"/*.root")

print(len(fileNames), " found")
# %%

events_list = []
pdgIDs_list = []
totalEvents = 0
for id, fileName in enumerate(fileNames[:4]):

#path = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_04_Run2_mc_124X.root"
    f = uproot.open(fileName)
    tree = f['tree']
    branches = tree.arrays()


    event = branches['event']
    pdgID = branches['pdgID']
    totalEvents = totalEvents + np.max(event)
    event = event+500*id
    print(event)
    events_list.append(event)
    pdgIDs_list.append(pdgID)

event = ak.concatenate(events_list, axis=0)
pdgID = ak.concatenate(pdgIDs_list, axis=0)

# %%
df = pd.DataFrame(pdgID)
df.columns=['one']
pdg_mapped = np.array(df.one.apply(map_to_groups_letter))


bins = np.arange(np.max(event))
c=np.histogram(event[(pdg_mapped=='B') | (pdg_mapped=='BB')], bins=bins)[0]
cmed = np.sum(c)/totalEvents
fig, ax = plt.subplots(1, 1)
ax.hist(c, bins=np.arange(7))
ax.set_title("Bottom vertices")
ax.text(x=0.9, y=0.9, s="Mean : %.3f"%cmed, transform=ax.transAxes, ha='right')
print("B gen Vertices", cmed)



# %%
c=np.histogram(event[(pdg_mapped=='D') | (pdg_mapped=='CB')], bins=bins)[0]
cmed = np.sum(c)/totalEvents
fig, ax = plt.subplots(1, 1)
ax.hist(c, bins=np.arange(7))
ax.set_title("Charmed vertices")
ax.text(x=0.9, y=0.9, s="Mean : %.3f"%cmed, transform=ax.transAxes, ha='right')
print("D gen Vertices", cmed)

# %%
c=np.histogram(event[(pdg_mapped=='SB')], bins=bins)[0]
cmed = np.sum(c)/totalEvents
fig, ax = plt.subplots(1, 1)
ax.hist(c, bins=np.arange(7))
ax.set_title("Other vertices")
ax.text(x=0.9, y=0.9, s="Mean : %.3f"%cmed, transform=ax.transAxes, ha='right')
print("Other gen Vertices", cmed)
# %%
