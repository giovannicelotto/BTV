import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math, sys
import awkward as ak
import random
from scipy.optimize import linear_sum_assignment
import mplhep as hep
hep.style.use("CMS")
from BTV.scripts.tuplizer.utilsForScript import distance_3d, getPdgMask
from helpers import getTreeAndBranches, criterion0, criterion1, eventDisplay, getTreeAndBranches


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
def main(nEvents, criterion, threshold):
    '''
    19/04: use 0 for iterative and -1 for threshold
    nEvents = number of lines to read from the ROOT file
    criterion = 0: iterative reduction
                1: global minimum
    threshold = (float) if all the entries of a row (column) are above the threshold delete it.
    '''
    params = {
        'threshold' : threshold,
        'eta_threshold': 2.5,
        'pt_threshold': 10,
    }

    print("Threshold for row /column : %.2f"%params['threshold'])
    print("Eta_threshold for mesons : %.2f"%params['eta_threshold'])
    print("Pt_threshold for mesons : %.2f"%params['pt_threshold'])
    
    tree, branches = getTreeAndBranches(fileName = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_2_Run2_mc_124X.root")
    
    GenPart_pdgId               = branches["GenPart_pdgId"]
    genVtxCouter = 0

    df = pd.DataFrame({
            'idx'     :[], 'pdgID'     :[], 'pt'        :[],
            'eta'       :[], 'phi'       :[], 'vx'        :[],
            'vy'        :[], 'vz'        :[], 'displacement':[],
            'matched'   :[], 'distance'  :[], 'normDistance':[]
        })
    nEvents = tree.num_entries if nEvents is None else nEvents
    for ev in np.arange(nEvents):
        #if ev>1000:
        #    break
        if ev%100==0:
            print("Event ", ev, "\r")
        #sys.exit("exit")
        
        nSV                         = branches["nSV"][ev]
        nGenPart                    = branches["nGenPart"][ev]
        GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId               = branches["GenPart_pdgId"][ev]
        GenPart_pt                  = branches["GenPart_pt"][ev]
        GenPart_eta                 = branches["GenPart_eta"][ev]
        GenPart_phi                 = branches["GenPart_phi"][ev]
        PV_x                        = branches["PV_x"][ev]
        PV_y                        = branches["PV_y"][ev]
        PV_z                        = branches["PV_z"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]
        SV_chi2                     = branches["SV_chi2"][ev]
        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]
        ProbeTracks_matchedToSV     = branches["ProbeTracks_matchedToSV"][ev]

        # filter the gen part of interest
        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId)
        etaMask = abs(GenPart_eta)<params['eta_threshold']
        ptMask = GenPart_pt>params['pt_threshold']
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask]
        daughterMask = GenPart_pdgId==GenPart_pdgId[GenPart_genPartIdxMother]     # particles w/ same pdgID of the mother
        # but we need to exclude mothers with the same pdgID of the daughters
        # GenPart_genPartIdxMother[daughterMask] is array of motherIdx with same pdg of daughter
        print(np.intersect1d(mesons, GenPart_genPartIdxMother[daughterMask]))
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask ]
        mesons = np.setdiff1d(mesons, GenPart_genPartIdxMother[daughterMask])
        
        df_event = pd.DataFrame({
            'idx':[],     'pdgID':[],     'pt':[],
            'eta':[],       'phi':[],       'vx':[],
            'vy':[],        'vz':[],        'displacement':[],
            'genDaughters':[],
            'matched':[],   'distance':[],  'normDistance':[],
            'probeTracksFromSV':[]
        })
        
         
        for mes in mesons:            
            nDaughters = sum((GenPart_genPartIdxMother == mes) & (abs(GenPart_pdgId)!=12) & (abs(GenPart_pdgId)!=14) & ((abs(GenPart_pdgId)!=16) ))
            #nDaughters = sum((GenPart_genPartIdxMother == mes))

            for gp in range(nGenPart):
                if (GenPart_genPartIdxMother[gp] == mes):
                    
                    row = [gp,
                           GenPart_pdgId[mes],
                           GenPart_pt[mes],
                           GenPart_eta[mes],
                           GenPart_phi[mes],
                           GenPart_vx[gp],
                           GenPart_vy[gp],
                           GenPart_vz[gp],
                           distance_3d( #displacement
                                        (GenPart_vx[gp], GenPart_vy[gp], GenPart_vz[gp]),
                                        (GenPart_vx[mes], GenPart_vy[mes], GenPart_vz[mes])),
                           nDaughters,
                           False,
                           999.,
                           999.,
                           0,# probetracks from SV
                           ]
                    
                    df_event.loc[len(df_event)]=row
                    
                    break
        genVtxCouter=genVtxCouter+len(df_event.index)
        df_event['pdgClass'] = df_event['pdgID'].apply(map_to_groups_letter)


    # Matching
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])
        tracksCounters = np.array([sum(ProbeTracks_matchedToSV == idx) for idx in range(nSV)])
        
        input("next")
        # row by row plots
        
        
        # matrix of distances
        distances = np.zeros((nSV, len(df_event.index)))
        for i in range(nSV):
            for j in range(len(df_event.index)):
                
                distances[i][j] = distance_3d(SVs[i], (df_event.vx[j], df_event.vy[j], df_event.vz[j]) )
        
        distances_normalized = distances.copy()
        #print(distances)
        for i in range(nSV):
            dist = distance_3d(SVs[i], (PV_x, PV_y, PV_z) )
            distances_normalized[i, :] = distances_normalized[i,:]/dist
        
        # remove row and columns with normalized distances > 1
        if params['threshold']==-1:
            params['threshold']=np.inf
        row_mask = np.all(distances_normalized > params['threshold'], axis=1)
        col_mask = np.all(distances_normalized > params['threshold'], axis=0)
        row_mask = np.logical_not(row_mask)
        col_mask = np.logical_not(col_mask)
        #print("Shape before ", distances.shape)
        distances = distances[row_mask][:, col_mask]
        #print("Shape after ", distances.shape)
        distances_normalized = distances_normalized[row_mask][:, col_mask]
        #print(distances)

        
        if criterion:
            df_event = criterion1(distances=distances, distances_normalized=distances_normalized, df_event=df_event, disaply=True, SVs=SVs,
                                  PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  GenPart_vy=GenPart_vy, col_mask=col_mask, row_mask=row_mask)
        else:
            df_event = criterion0(distances=distances, distances_normalized=distances_normalized, df_event=df_event, display=True,
                                  SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  tracksCounters=tracksCounters,
                                  GenPart_vy=GenPart_vy, col_mask=col_mask, row_mask=row_mask)

        if len(df_event)!=0:
            df = pd.concat((df, df_event))

    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    threshold = float(sys.argv[3])
    print(criterion)
    main(nEvents, criterion, threshold)