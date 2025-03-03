#%%
import uproot
import numpy as np
import sys
import matplotlib.pyplot as plt
filePath = "/work/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root"
file = uproot.open(filePath)
sys.path.append('/work/gcelotto/BTV/scripts/tuplizer/helpers')
from getGV import getMesons_old
from getParams import getParams
#%%
tree = file['Events']
branches = tree.arrays()
# %%
bins=np.arange(-1, 10)
cTot = np.zeros(len(bins)-1)
num0    = 0
num1    = 0
num2    = 0
num3    = 0
nSV_tot = 0
nGV_tot = 0
nSV1_tot = 0
nSV2_tot = 0

tracksPerSV = np.array([])

for ev in range(900):#tree.num_entries:
    if ev % 100 == 0:
        print(ev)
    ProbeTracks_pt              = branches["ProbeTracks_pt"][ev]
    ProbeTracks_eta             = branches["ProbeTracks_eta"][ev]
    ProbeTracks_genPartIdx      = branches["ProbeTracks_genPartIdx"][ev]
    GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
    GenPart_pdgId               = branches["GenPart_pdgId"][ev]
    GenPart_eta                 = branches["GenPart_eta"][ev]
    GenPart_pt                  = branches["GenPart_pt"][ev]
    ProbeTracks_matchedToSV     = branches["ProbeTracks_matchedToSV"][ev]
    GenPart_pdgId               = branches["GenPart_pdgId"][ev]
    nSV                         = branches["nSV"][ev]
    nGenPart                    = branches["nGenPart"][ev]
    GenPart_status              = branches["GenPart_status"][ev]
    GenPart_charge              = branches["GenPart_charge"][ev]

    params = getParams()
    mesons = getMesons_old(params, GenPart_pdgId, GenPart_eta, GenPart_pt, nGenPart, GenPart_status, GenPart_charge, GenPart_genPartIdxMother)     
    nGV_event = len(mesons)
    nGV_tot = nGV_tot + nGV_event
    
    # SV con almeno una traccia da SV con pt>0.8 e eta <2.5
    # set elimina i doppioni
    nSV1track_event = len(set(ProbeTracks_matchedToSV[(ProbeTracks_matchedToSV>-1) & (ProbeTracks_pt>0.8) & (abs(ProbeTracks_eta)<2.5)]))

    nSV1_tot = nSV1_tot + nSV1track_event
    countsPerSV = np.histogram(ProbeTracks_matchedToSV[(ProbeTracks_matchedToSV>-1) & (ProbeTracks_pt>0.8) & (abs(ProbeTracks_eta)<2.5)], bins=np.arange(10))[0]
    tracksPerSV = np.concatenate((tracksPerSV, countsPerSV))

    num1 = num1 + np.min((nSV1track_event, nGV_event))

    # SV con almeno una traccia da SV con pt>0.8 e eta <2.5 & genMatchata ad una genPart
    nSV1trackMatched_event =            len(set(ProbeTracks_matchedToSV[(ProbeTracks_matchedToSV>-1)& (ProbeTracks_genPartIdx>-1) & (ProbeTracks_pt>0.8) & (abs(ProbeTracks_eta)<2.5)]))
    num2 = num2 + np.min((nSV1trackMatched_event, nGV_event))
    nSV2_tot = nSV2_tot + nSV1trackMatched_event

    maskGreen = np.array((ProbeTracks_matchedToSV>-1)& (ProbeTracks_genPartIdx>-1) & (ProbeTracks_pt>0.8) & (abs(ProbeTracks_eta)<2.5))
    maskGreen2 = (np.isin(GenPart_genPartIdxMother[ProbeTracks_genPartIdx],mesons)) | (np.isin(GenPart_genPartIdxMother[GenPart_genPartIdxMother[ProbeTracks_genPartIdx]],mesons)) | (np.isin(GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[ProbeTracks_genPartIdx]]],mesons))
    maskGreen = (maskGreen) & (maskGreen2)
    nSV1trackMatchedToDaughter_event =  len(set(maskGreen))
    num3 = num3 + np.min((nSV1trackMatchedToDaughter_event, nGV_event))


    # den Tot is the total number of SV

    num0 = num0 + np.min((nSV,nGV_event))
    nSV_tot = nSV_tot + nSV

# %%
print("Max reachable efficiency : %.1f%%"%(num0/nGV_tot*100))
print("Blue eff : %.1f%%"%(num1/nGV_tot*100))
print("Blue*Red eff : %.1f%%"%(num2/nGV_tot*100))
print("Blue*Red*Green eff : %.1f%%"%(num3/nGV_tot*100))

print("SV with >=1 track from SV    : %.2f%%"%(nSV1_tot*100/nSV_tot))
print("SV with >=1 track genMatched : %.2f%%"%(nSV2_tot*100/nSV_tot))
# %%
