import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math, sys
import awkward as ak
from scipy.optimize import linear_sum_assignment
import mplhep as hep
hep.style.use("CMS")
def distance_3d(point1, point2):
    """
    Calculates the distance between two points in three-dimensional space.
    
    Args:
        point1 (tuple): The coordinates of the first point (x, y, z).
        point2 (tuple): The coordinates of the second point (x, y, z).
        
    Returns:
        float: The distance between the two points.
    """
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

def main():
    f = uproot.open("/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTTo2L2Nu_Run2_mc_124X.root")
    tree = f['Events']
    branches = tree.arrays()
    GenPart_pdgId               = branches["GenPart_pdgId"]
    Jet_pt                      = branches["Jet_pt"]
    
    mask1 = ak.any((abs(GenPart_pdgId) == 521) |
                   (abs(GenPart_pdgId) == 511) |
                   (abs(GenPart_pdgId) == 411) |
                   (abs(GenPart_pdgId) == 421) |
                   ((abs(GenPart_pdgId) > 4100) & (abs(GenPart_pdgId) < 4500)) |        # charmed baryons, sigma_c, xi_c, lambda_c 
                   ((abs(GenPart_pdgId) > 3100) & (abs(GenPart_pdgId) < 3340))          # strange baryons, lambda, csi
                   , axis=-1)
    print("Events with a meson in the event : ", np.sum(mask1)/tree.num_entries*100, " %")
    genVtxCouter = 0
    matchedDistances = []
    mask = (mask1)
    for ev in np.arange(tree.num_entries)[mask]:
        #if ev>1000:
        #    break
        print("Event ", ev)
        #sys.exit("exit")
        
        nSV                         = branches["nSV"][ev]
        nGenPart                    = branches["nGenPart"][ev]
        GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId               = branches["GenPart_pdgId"][ev]
        SV_pt                      = branches["SV_pt"][ev]
        SV_dxy                      = branches["SV_dxy"][ev]
        SV_dxySig                   = branches["SV_dxySig"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]
        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]

        mesons = np.arange(nGenPart)[np.isin(GenPart_pdgId, [521, -521, 511, -511, 411, 421, -411, -421])]
        #mesons = np.arange(nGenPart)[np.isin(GenPart_pdgId, [ 411, 421, -411, -421])]
        #mesons = np.arange(nGenPart)[np.isin(GenPart_pdgId, [ 521, -521, 511, -511])]
        daughters=[]            # all the daughters of 5 degrees

        for mes in mesons:            
            for gp in range(nGenPart):
                if (GenPart_genPartIdxMother[gp] == mes):
                    daughters.append(gp)                        # save one daughter for each b mesons
                    break
              
        SVs     =   []
        
        genVs   =   []
        for svIdx in range(nSV):
            sv = (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx])
            SVs.append(sv)  
            
        
        genVtxCouter=genVtxCouter+len(daughters)
        for genIdx in daughters:
            genPartVert = (GenPart_vx[genIdx], GenPart_vy[genIdx], GenPart_vz[genIdx])
            genVs.append(genPartVert)
                
        #if nSV == len(daughters):
        distances = np.zeros((nSV, len(daughters)))
        for i in range(nSV):
            for j in range(len(daughters)):
                distances[i][j] = distance_3d(SVs[i], genVs[j])
        while (distances.shape[0]>0) & (distances.shape[1]>0):
            minIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
            matchedDistances.append(distances[minIdx[0], minIdx[1]])
            #print(distances)
            #print(distances[minIdx[0], minIdx[1]])
            #delete the row and column:
            distances = np.delete(distances, minIdx[0], axis=0)
            distances = np.delete(distances, minIdx[1], axis=1)
            #print(distances)
            #input("Next")

        #print("\n\n\n")
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    bins= np.linspace(0, 5, 100)
    hep.cms.label(ax=ax)
    ax.hist(np.clip(matchedDistances, bins[0], bins[-1]), bins=bins)
    ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(matchedDistances)), ha='right', transform=ax.transAxes)
    
    ax.set_xlabel("GenMatched 3d distance [cm]", fontsize=15)
    ax.set_yscale('log')
    fig.savefig("/t3home/gcelotto/BTV/plots/genMatched_distance_cr2.png", bbox_inches='tight')
    plt.close('all')

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    thresholds = np.linspace(0, 5, 100)
    y = []
    hep.cms.label(ax=ax)
    for i in thresholds:
        y.append(np.sum(matchedDistances<i)/len(matchedDistances)*100)
    
    ax.plot(thresholds, y, color='blue')
    ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(matchedDistances)), ha='right', transform=ax.transAxes)
    ax.set_xlabel("Threshold of matching [cm]")
    ax.set_ylabel("Efficiency [%]")
    ax.set_xlim(bins[0], bins[-1])#thresholds[-1])
    ax.set_ylim(60, 100)
    ax.grid(True)
    #ax.set_ylim(0.1, 1)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig("/t3home/gcelotto/BTV/plots/Eff_vs_threshold_cr2.png", bbox_inches='tight')
    print(len(matchedDistances)/genVtxCouter) 
    return 

if __name__ =="__main__":
    main()