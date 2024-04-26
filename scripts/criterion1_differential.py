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
from helpers import getTreeAndBranches, criterion0, criterion1, getTreeAndBranches

def getSuffix(criterion, threshold):
    import math
    print(threshold, math.modf(threshold))
    floatPart, intPart = math.modf(threshold)[0], int(math.modf(threshold)[1])
    floatPart = int(floatPart*10)
    if floatPart == 0:
        suffix = "cr%d_t%d"%(criterion, intPart)
    else:
        suffix = "cr%d_t%dp%d"%(criterion, intPart, floatPart)
    
    return suffix


def main(nEvents, criterion, threshold):
    '''
    nEvents = number of lines to read from the ROOT file
    criterion = 0: iterative reduction
                1: global minimum
                2: itearative reduction with no cut
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
            'index'     :[], 'pdgID'     :[], 'pt'        :[],
            'eta'       :[], 'phi'       :[], 'vx'        :[],
            'vy'        :[], 'vz'        :[], 'displacement':[],
            'matched'   :[], 'distance'  :[], 'normDistance':[]
        })
    nEvents = tree.num_entries if nEvents is None else nEvents
    suffix = getSuffix(criterion, params['threshold'])
    for ev in np.arange(nEvents):
        if ev%100==0:
            print("Event ", ev, "\r")
        
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
        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]
        ProbeTracks_matchedToSV     = branches["ProbeTracks_matchedToSV"][ev]
        ProbeTracks_pt              = branches["ProbeTracks_pt"][ev]
        ProbeTracks_eta             = branches["ProbeTracks_eta"][ev]

        # filter the gen part of interest
        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId)   # list of pdgId consider as genVertices
        etaMask = abs(GenPart_eta)<params['eta_threshold']  # limit the genvertices to GenHadrons within tracker acceptance
        ptMask = GenPart_pt>params['pt_threshold']          # require pT of GenHadrons larger than 10 GeV
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask]
        
        df_event = pd.DataFrame({
            'index':[],     'pdgID':[],     'pt':[],
            'eta':[],       'phi':[],       'vx':[],
            'vy':[],        'vz':[],        'displacement':[],
            'genDaughters':[],
            'matched':[],   'distance':[],  'normDistance':[],
            'probeTracksFromSV':[]
        })
        
         
        for mes in mesons:            
            nDaughters = sum((GenPart_genPartIdxMother == mes) & (abs(GenPart_pdgId)!=12) & (abs(GenPart_pdgId)!=14) & (abs(GenPart_pdgId)!=16) & (GenPart_pt>10) & (abs(GenPart_eta)<2.5) )
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

    # Matching
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])
        tracksCounters = np.array([sum((ProbeTracks_matchedToSV == idx) & (ProbeTracks_pt > 10) & (abs(ProbeTracks_eta)<2.5)) for idx in range(nSV)])

        
        # matrix of distances
        #distances = np.zeros((nSV, len(df_event.index)))
        #for i in range(nSV):
        #    for j in range(len(df_event.index)):
        #        distances[i][j] = distance_3d(SVs[i], (df_event.vx[j], df_event.vy[j], df_event.vz[j]) )
        distances = np.array([[distance_3d(sv, (vx, vy, vz)) for vx, vy, vz in zip(df_event.vx, df_event.vy, df_event.vz)] for sv in SVs ])
        # distances_normalized = distances / np.linalg.norm(SVs - np.array([[PV_x, PV_y, PV_z]]), axis=1).reshape(-1, 1)
        #print(distances)
        distances_normalized = distances.copy()

        for i in range(nSV):
            dist = distance_3d(SVs[i], (PV_x, PV_y, PV_z) )
            distances_normalized[i, :] = distances_normalized[i,:]/dist
        
        # remove row and columns with normalized distances > 1
        if params['threshold']==-1:
            params['threshold']=np.inf
        row_mask = np.all(distances_normalized > params['threshold'], axis=1) if nSV> 0 else None
        col_mask = np.all(distances_normalized > params['threshold'], axis=0) if nSV> 0 else None
        row_mask = np.logical_not(row_mask) if nSV> 0 else None
        col_mask = np.logical_not(col_mask) if nSV>0 else None
        distances = distances[row_mask][:, col_mask]
        distances_normalized = distances_normalized[row_mask][:, col_mask]

        
        if criterion==1:
        # if one reco/gen vertex has all the possible matching > 2.5 exlude that row or column
            df_event = criterion1(distances=distances, distances_normalized=distances_normalized, df_event=df_event, disaply=False, SVs=SVs,
                                  PV_x=None, PV_y=None, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  GenPart_vy=GenPart_vy, col_mask=col_mask, row_mask=row_mask)
        else:
            df_event = criterion0(distances=distances, distances_normalized=distances_normalized, df_event=df_event, display=False,
                                  SVs=SVs, SV_chi2=None, PV_x=None, PV_y=None, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  GenPart_vy=GenPart_vy, tracksCounters=tracksCounters, col_mask=col_mask, row_mask=row_mask)
            
        if not df_event.empty:
            df = pd.concat([df, df_event], ignore_index=True)
    
    outName = "/t3home/gcelotto/BTV/output/df_"+suffix+".parquet"
    df.to_parquet(outName)
    print("Saved in ", outName)
    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    threshold = float(sys.argv[3])
    main(nEvents, criterion, threshold)