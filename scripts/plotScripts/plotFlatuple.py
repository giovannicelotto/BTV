# %%
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.colors import LogNorm
hep.style.use("CMS")
import awkward as ak
import numpy as np
# %%
path = "/work/gcelotto/BTV/scripts/tuplizer/TTToH_selPart_100.root"
f = uproot.open(path)
tree = f['tree']
branches = tree.arrays()

binning ={
#'event'             : [0,1e4,1e4]    ,
'nSV'               : [0,10,11]    ,
'pdgID'             : [-600,600,10001]    ,
'motherPdgID'       : [-600,600,10001]    ,
'pt'                : [0,20,51]    ,
'eta'               : [-2.7,2.7,51]    ,
'phi'               : [-3.5,3.5,51]    ,
'vx'                : [-0.1,0.1,51]    ,
'vy'                : [-0.1,0.1,51]    ,
'vz'                : [-20,20,51]    ,
'displacement'      : [0,1,51]            ,
'promptStatus'      : [-1,5,51]            ,
'decayStatus'       : [-1,5,51]            ,
#'B_eta'             : [3,3,51]    ,
#'B_phi'             : [4,4,51]    ,
'nGenTracks'        : [0,10,11]            ,
'genTracks_pt'      : [0,20,51]            ,
'genTracks_eta'     : [-3,3,51]            ,
'genTracks_phi'     : [-4,4,51]            ,
'genTracks_pdgId'   : [-600,600,51]                ,
'genTracks_charge'  : [-1,1,51]                ,
'nRecoTracks'       : [0,10,11]            ,
'recoTracks_pt'     : [0,20,51]            ,
'recoTracks_eta'    : [-3,3,51]                ,
'recoTracks_phi'    : [-4,4,51]                ,
'recoTracks_pdgID'  : [-600,600,51]                ,
'distance'          : [-1,1,51]        ,
'delta_x'           : [-0.1,0.1,51]        ,
'delta_y'           : [-0.1,0.1,51]        ,
'delta_z'           : [-0.1,0.1,51]        ,
'SV_dlenSig'        : [2,7,51]            ,

}
branches_name = binning.keys()
# %%
# Determine the number of branches
num_branches = len(branches_name)
cols = 6  # Number of columns for subplots
rows = (num_branches // cols) + (num_branches % cols > 0)  # Calculate rows needed
# %%
# Create subplots
fig, axes = plt.subplots(rows, cols, figsize=(30, 4 * rows))
axes = axes.flatten()  # Flatten in case we have 2D array of subplots

for i, branch_name in enumerate(branches_name):
    values = branches[branch_name]
    ax = axes[i]
    
    try:
        ax.hist(values, bins=np.linspace(binning[branch_name][0], binning[branch_name][1], int(binning[branch_name][2])), alpha=0.7, color='blue', edgecolor='black')
    except:
        ax.hist(ak.flatten(values), bins=np.linspace(binning[branch_name][0], binning[branch_name][1], int(binning[branch_name][2])), alpha=0.7, color='blue', edgecolor='black')
    ax.set_title(branch_name, fontsize=10)
    ax.set_xlabel("Value")
    ax.set_ylabel("Frequency")

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()

# %%
