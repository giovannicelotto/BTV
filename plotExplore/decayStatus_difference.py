# %%
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import awkward as ak
# %%
path = "/work/gcelotto/BTV/scripts/tuplizer/TTToH_selPart_100.root"
f = uproot.open(path)
tree = f['tree']
branches=tree.arrays()
# %%
genTracks_pdgId = branches["genTracks_pdgId"]
decayStatus = branches["decayStatus"]
nGenTracks = branches["nGenTracks"]
# %%
arr_0 = ak.ravel(genTracks_pdgId[(decayStatus == 0) & (nGenTracks==2)])
arr_1 = ak.ravel(genTracks_pdgId[(decayStatus == 1) & (nGenTracks==2)])
arr_2 = ak.ravel(genTracks_pdgId[(decayStatus == 2) & (nGenTracks==2)])
arr_3 = ak.ravel(genTracks_pdgId[(decayStatus == 3) & (nGenTracks==2)])

# Count occurrences for each array
values_0, counts_0 = np.unique(arr_0, return_counts=True)
values_1, counts_1 = np.unique(arr_1, return_counts=True)
values_2, counts_2 = np.unique(arr_2, return_counts=True)
values_3, counts_3 = np.unique(arr_3, return_counts=True)
counts_0 = counts_0/np.sum(counts_0)
counts_1 = counts_1/np.sum(counts_1)
counts_2 = counts_2/np.sum(counts_2)
counts_3 = counts_3/np.sum(counts_3)
# Convert values to integers
values_0 = np.array(values_0, dtype=int)
values_1 = np.array(values_1, dtype=int)
values_2 = np.array(values_2, dtype=int)
values_3 = np.array(values_3, dtype=int)

# Create a common x-axis ensuring both sets share the same columns
all_values = np.union1d(values_1, values_2)
all_values = np.union1d(all_values, values_3)
all_values = np.union1d(all_values, values_0)

# Map occurrences to the common x-axis
counts_dict_0 = dict(zip(values_0, counts_0))
counts_dict_1 = dict(zip(values_1, counts_1))
counts_dict_2 = dict(zip(values_2, counts_2))
counts_dict_3 = dict(zip(values_3, counts_3))

# Fill missing values with zeros
counts_0_aligned = [counts_dict_0.get(v, 0) for v in all_values]
counts_1_aligned = [counts_dict_1.get(v, 0) for v in all_values]
counts_2_aligned = [counts_dict_2.get(v, 0) for v in all_values]
counts_3_aligned = [counts_dict_3.get(v, 0) for v in all_values]

# Create the plot
plt.figure(figsize=(10, 5))

# Hollow bars for visibility
plt.bar(range(len(all_values)), counts_0_aligned, edgecolor='red', facecolor='none',   label='B To Other')
plt.bar(range(len(all_values)), counts_1_aligned, edgecolor='blue', facecolor='none',  label='B To D'      )
plt.bar(range(len(all_values)), counts_2_aligned, edgecolor='green', facecolor='none', label='B To X To D' )
plt.bar(range(len(all_values)), counts_3_aligned, edgecolor='purple', facecolor='none', label='B To X To Y To D' )

# Set x-ticks to correspond to the unique values
plt.xticks(range(len(all_values)), all_values, rotation=90)

plt.ylabel("Frequency")
plt.legend()
plt.xlabel("PDG ID of GenTracks from GV")
#plt.title("Comparison of genTracks_pdgId for Different Decay Status")

# Show plot
plt.show()

# %%


nGenTracks = branches["nGenTracks"]
decayStatus = branches["decayStatus"]
# %%
fig, ax = plt.subplots(1, 1)
bins=np.arange(10)
c0 =np.histogram(nGenTracks[decayStatus==0], bins=bins)[0]
c1 =np.histogram(nGenTracks[decayStatus==1], bins=bins)[0]
c2 =np.histogram(nGenTracks[decayStatus==2], bins=bins)[0]
c3 =np.histogram(nGenTracks[decayStatus==3], bins=bins)[0]
c0=c0/np.sum(c0)
c1=c1/np.sum(c1)
c2=c2/np.sum(c2)
c3=c3/np.sum(c3)
ax.hist(bins[:-1], bins=bins, weights=c0, label='B To Other',   histtype='step', color='red', linewidth=2)
ax.hist(bins[:-1], bins=bins, weights=c1, label='B To D',       histtype='step', color='blue', linewidth=2)
ax.hist(bins[:-1], bins=bins, weights=c2, label='B To X To D',  histtype='step', color='green', linewidth=2)
ax.hist(bins[:-1], bins=bins, weights=c3, label='B To X To D',  histtype='step', color='purple', linewidth=2)
ax.set_xlabel("Number of genTracks")
ax.legend()
# %%
