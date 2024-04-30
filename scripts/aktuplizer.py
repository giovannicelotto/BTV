import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import mplhep as hep
hep.style.use("CMS")
from utilsForScript import distance_3d, getPdgMask
from helpers import getTreeAndBranches, criterion0_new, criterion1, getTreeAndBranches
import awkward as ak

def getSuffix(criterion, threshold):
    import math
    floatPart, intPart = math.modf(threshold)[0], int(math.modf(threshold)[1])
    floatPart = int(floatPart*10)
    if floatPart == 0:
        suffix = "cr%d_t%d"%(criterion, intPart)
    else:
        suffix = "cr%d_t%dp%d"%(criterion, intPart, floatPart)
    
    return suffix


def main(nEvents, criterion):
    '''
    nEvents = number of lines to read from the ROOT file
    criterion = 0: iterative reduction  
                1: global minimum       #not implemented
                2: itearative reduction with no cut #not implemented
    '''
    params = {
        'threshold':-1,
        'eta_tracker': 2.5,
        'minPt_meson': 10,
        'minPt_daughter':5
    }
    df = None
    for key in params.keys():
        print(key, " : ", params[key])

    tree, branches = getTreeAndBranches(fileName = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_2_Run2_mc_124X.root")
    nEvents = tree.num_entries if nEvents is None else nEvents
    suffix = getSuffix(criterion, params['threshold'])
    
    GenPart_pdgId               = branches["GenPart_pdgId"]

    
    
    for ev in np.arange(nEvents):
        if ev%1==0:
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
        etaMask = abs(GenPart_eta)<params['eta_tracker']  # limit the genvertices to GenHadrons within tracker acceptance
        ptMask = GenPart_pt>params['minPt_meson']          # require pT of GenHadrons larger than 10 GeV
        #daughterMask = GenPart_pdgId==GenPart_pdgId[GenPart_genPartIdxMother]     # particles w/ same pdgID of the mother
        # but we need to exclude mothers with the same pdgID of the daughters
        # GenPart_genPartIdxMother[daughterMask] is array of motherIdx with same pdg of daughter
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask ]
        #mesons = np.setdiff1d(mesons, GenPart_genPartIdxMother[daughterMask])
        
        #df_event = getEmptyDf()
        #assert len(df_event)==0
        rows=[]
        for mes in mesons:        
            notNeutrinoMask = (abs(GenPart_pdgId)!=12) & (abs(GenPart_pdgId)!=14) & (abs(GenPart_pdgId)!=16)    
            #nDaughters = sum((GenPart_genPartIdxMother == mes) & (notNeutrinoMask) & (GenPart_pt>daughterPt) & (abs(GenPart_eta)<2.5) )
            daughters_pt = GenPart_pt[(GenPart_genPartIdxMother == mes) &  (notNeutrinoMask) & (abs(GenPart_eta)<2.5)]
            daughters_pdgId = GenPart_pdgId[(GenPart_genPartIdxMother == mes) &  (notNeutrinoMask) & (abs(GenPart_eta)<2.5)]
            #nDaughters = sum((GenPart_genPartIdxMother == mes))

            dau = np.arange(nGenPart)[GenPart_genPartIdxMother == mes][0] if np.sum(GenPart_genPartIdxMother == mes)>0 else -1
            if dau==-1:
                continue
            row={'idx': mes, 'pdgID': GenPart_pdgId[mes], 'pt': GenPart_pt[mes], 'eta': GenPart_eta[mes], 'phi': GenPart_phi[mes], 'vx': GenPart_vx[dau], 'vy': GenPart_vy[dau], 'vz': GenPart_vz[dau],
                'displacement': distance_3d((GenPart_vx[dau], GenPart_vy[dau], GenPart_vz[dau]), (GenPart_vx[mes], GenPart_vy[mes], GenPart_vz[mes])),
                'daughters_pt': daughters_pt,'daughters_pdgId': daughters_pdgId, 'matched': False, 'distance': 999.1, 'probeTracksFromSV':0}
            rows.append(row)


    # Matching
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])
        tracksCounters = np.array([sum((ProbeTracks_matchedToSV == idx) & (ProbeTracks_pt > params['minPt_daughter']) & (abs(ProbeTracks_eta)<2.5)) for idx in range(nSV)])
        
    # Matrix of distances
        distances = np.zeros((nSV, len(rows)))
        for i in range(nSV):
            for j in range(len(rows)):
                distances[i][j] = distance_3d(SVs[i], (rows[j]['vx'], rows[j]['vy'], rows[j]['vz']) )
        
        if criterion==1:
        # if one reco/gen vertex has all the possible matching > 2.5 exlude that row or column
            pass
        else:
            rows = criterion0_new(distances, rows, tracksCounters, ProbeTracks_pt)
            
        for row in rows:
            array = ak.Array(row)
            if df is None:
                df=array
                print(type(df))
            else:
                print(array)
                df = ak.concatenate([df, array])
        #if not df_event.empty:
        #    pass
            #df = pd.concat([df, df_event], ignore_index=True)
    
    outName = "/t3home/gcelotto/BTV/output/df_"+suffix+".parquet"
    #df['daughters_pt'] = df['daughters_pt'].apply(lambda x: ','.join(map(str, x)))
    #df['daughters_pdgId'] = df['daughters_pdgId'].apply(lambda x: ','.join(map(str, x)))
    #df['probeTracks_pt'] = df['probeTracks_pt'].apply(lambda x: ','.join(map(str, x)))
    #df.to_parquet(outName)
    print("Saved in ", outName)
    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    main(nEvents, criterion)