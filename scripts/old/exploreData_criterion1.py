import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math, sys
import awkward as ak
from scipy.optimize import linear_sum_assignment
import mplhep as hep
hep.style.use("CMS")
from utilsForScript import distance_3d, getPdgMask
def main():
    params = {
        'threshold' : 1,
        'eta_threshold': 2.5
    }

    print("Threshold for row /column : %.2f"%params['threshold'])
    f = uproot.open("/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root")
    tree = f['Events']
    branches = tree.arrays()
    GenPart_pdgId               = branches["GenPart_pdgId"]
    
    # mask to filter events where there are the mesons of interets (aka if an event has no genparticles with these numbers, skip it)
    mask = ak.any(
        (abs(GenPart_pdgId) == 521) | (abs(GenPart_pdgId) == 511) | (abs(GenPart_pdgId) == 531) | (abs(GenPart_pdgId) == 541) |
        (abs(GenPart_pdgId) == 411) | (abs(GenPart_pdgId) == 421) | (abs(GenPart_pdgId) == 431) |
        ((abs(GenPart_pdgId) > 3100) & (abs(GenPart_pdgId) < 3340)) |        #  strange baryons
        ((abs(GenPart_pdgId) > 4100) & (abs(GenPart_pdgId) < 4500)) |        #  charmed baryons
        ((abs(GenPart_pdgId) > 5100) & (abs(GenPart_pdgId) < 5600))          #  bottom baryons
                   , axis=-1)
    print("Events with a meson in the event : ", np.sum(mask1)/tree.num_entries*100, " %")
    
    genVtxCouter = 0
    matchedDistances = []
    
    for ev in np.arange(tree.num_entries)[mask]:
        #if ev>1000:
        #    break
        print("Event ", ev)
        #sys.exit("exit")
        
        nSV                         = branches["nSV"][ev]
        nGenPart                    = branches["nGenPart"][ev]
        GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId               = branches["GenPart_pdgId"][ev]
        GenPart_eta                 = branches["GenPart_eta"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]
        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]
        pdgMask = getPdgMask()
        etaMask = abs(GenPart_eta)<2.5
        mesons = np.arange(nGenPart)[pdgMask & etaMask]
        
        
        daughters=[] 
        for mes in mesons:            
            for gp in range(nGenPart):
                if (GenPart_genPartIdxMother[gp] == mes):
                    daughters.append(gp)                        # save one daughter for each b mesons
                    break
        genVtxCouter=genVtxCouter+len(daughters)
        genVs   =   []
        for genIdx in daughters:
            genPartVert = (GenPart_vx[genIdx], GenPart_vy[genIdx], GenPart_vz[genIdx])
            genVs.append(genPartVert)
              
        SVs     =   []
        for svIdx in range(nSV):
            sv = (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx])
            SVs.append(sv)  

                
        # matrix of distances
        distances = np.zeros((nSV, len(daughters)))
        for i in range(nSV):
            for j in range(len(daughters)):
                distances[i][j] = distance_3d(SVs[i], genVs[j])

        
        # threshold = 2.5 cm of distance
        # if one reco/gen vertex has all the possible matching > 2.5 exlude that row or column
        row_mask = np.all(distances > params['threshold'], axis=1)
        col_mask = np.all(distances> params['threshold'], axis=0)
        row_mask = np.logical_not(row_mask)
        col_mask = np.logical_not(col_mask)
        distances = distances[row_mask][:, col_mask]

        #print(distances)
        
        row_ind, col_ind = linear_sum_assignment(distances)
        for i in range(len(row_ind)):
            #print(f"SV {row_ind[i]} is assigned to Gen Decay Vertex {col_ind[i]} with a distance {distances[row_ind[i], col_ind[i]]}")
            matchedDistances.append(distances[row_ind[i], col_ind[i]])

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    bins= np.linspace(0, 5, 100)
    hep.cms.label(ax=ax)
    ax.hist(np.clip(matchedDistances, bins[0], bins[-1]), bins=bins)
    ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(matchedDistances)), ha='right', transform=ax.transAxes)
    if params['threshold']<999:
        print("Adding text")
        ax.text(x=0.98, y=.75, s="Threshold : %.1f cm"%(params['threshold']), ha='right', transform=ax.transAxes)
    ax.set_xlabel("GenMatched 3d distance [cm]", fontsize=15)
    ax.set_yscale('log')
    fig.savefig("/t3home/gcelotto/BTV/plots/genMatched_distance_cr1.png", bbox_inches='tight')
    plt.close('all')

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    thresholds = bins
    y = []
    hep.cms.label(ax=ax)
    for i in thresholds:
        y.append(np.sum(matchedDistances<i)/len(matchedDistances)*100)
    
    ax.plot(thresholds, y)
    ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(matchedDistances)), ha='right', transform=ax.transAxes)
    ax.set_xlabel("Threshold of matching [cm]")
    ax.set_ylabel("Efficiency [%]")
    ax.set_xlim(bins[0], bins[-1])#thresholds[-1])
    ax.set_ylim(60, 100)
    ax.grid(True)
    #ax.set_ylim(0.1, 1)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig("/t3home/gcelotto/BTV/plots/Eff_vs_threshold_cr1_%d.png"%bins[-1], bbox_inches='tight')


    print(len(matchedDistances)/genVtxCouter) 
    return 

if __name__ =="__main__":
    main()