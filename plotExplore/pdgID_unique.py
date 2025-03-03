# %%
import uproot
import awkward as ak
import numpy as np
#path = "/work/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root"
path = "/work/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Pions.root"

# %%
file = uproot.open(path)
tree = file["Events"]
# %%
branches = tree.arrays()
# %%
svDaughters_genPartIdx = branches["svDaughters_genPartIdx"]
GenPart_pdgId = branches["GenPart_pdgId"]
# %%
from collections import Counter
import matplotlib.pyplot as plt

# Flatten and filter PDG IDs
pdg_ids = ak.flatten(GenPart_pdgId[svDaughters_genPartIdx][svDaughters_genPartIdx != -1])
pdg_counts = Counter(pdg_ids)  # Count occurrences of each PDG ID

# Sort PDG IDs for better grouping
sorted_pdg_ids = sorted(pdg_counts.keys())  
occurrences = np.array([pdg_counts[pdg] for pdg in sorted_pdg_ids])
total_count = occurrences.sum()
percentages = (occurrences / total_count) * 100  

# Create figure and axis
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
bars = ax.bar(range(len(sorted_pdg_ids)), occurrences, tick_label=sorted_pdg_ids)

# Annotate bars with normalized percentages
for bar, percent in zip(bars, percentages):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),  
            f"{percent:.2f}%", ha='center', va='bottom', fontsize=9, color='black')

# Improve label readability
ax.set_xticks(range(len(sorted_pdg_ids)))  # Ensure labels are aligned correctly
ax.set_xticklabels(sorted_pdg_ids, rotation=90, fontsize=10)  # Rotate for readability
ax.set_xlabel("PDG ID")
ax.set_ylabel("Occurrences")
#ax.set_title("Histogram of PDG IDs")
ax.set_yscale('log')
# Show the plot
plt.show()

# %%
