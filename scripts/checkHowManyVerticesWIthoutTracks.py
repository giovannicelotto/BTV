# %%
import uproot
import numpy as np
# %%
f1="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/240723_073316/0000/TTToHadronic_Run2_mc_2024Jul23_537.root"
f1=uproot.open(f1)

# %%
stepSVtoTracks = 0
stepTrackstoGenPart =0
branches = f1['Events'].arrays()
for ev in range(100):

    ProbeTracks_matchedToSV = branches['ProbeTracks_matchedToSV'][ev]
    nSV                     = branches["nSV"][ev]
    nProbeTracks            = branches['nProbeTracks'][ev]
    ProbeTracks_pt          = branches["ProbeTracks_pt"][ev]
    ProbeTracks_eta         = branches["ProbeTracks_eta"][ev]
    ProbeTracks_genPartIdx  = branches["ProbeTracks_genPartIdx"][ev]

    for sv in range(nSV):
        if (len(set(ProbeTracks_matchedToSV[((ProbeTracks_matchedToSV==sv) & (ProbeTracks_pt>1) & (abs(ProbeTracks_eta)<2.5))]))>=1):
            stepSVtoTracks = stepSVtoTracks +1
        
        if (len(set(ProbeTracks_matchedToSV[((ProbeTracks_matchedToSV==sv) & (ProbeTracks_pt>1) & (abs(ProbeTracks_eta)<2.5) & (ProbeTracks_genPartIdx>-1))]))>=1):
            stepTrackstoGenPart = stepTrackstoGenPart +1

print(stepSVtoTracks/(np.sum(branches["nSV"])))
print(stepTrackstoGenPart/(np.sum(branches["nSV"])))
# %%
