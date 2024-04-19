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
def getTreeAndBranches(fileName):
    f = uproot.open(fileName)
    tree = f['Events']
    branches = tree.arrays()
    return tree, branches

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
    tree, branches = getTreeAndBranches(fileName = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root")
    
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

        # filter the gen part of interest
        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId)   # list of pdgId consider as genVertices
        etaMask = abs(GenPart_eta)<params['eta_threshold']  # limit the genvertices to GenHadrons within tracker acceptance
        ptMask = GenPart_pt>params['pt_threshold']          # require pT of GenHadrons larger than 10 GeV
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask]
        
        df_event = pd.DataFrame({
            'index':[],     'pdgID':[],     'pt':[],
            'eta':[],       'phi':[],       'vx':[],
            'vy':[],        'vz':[],        'displacement':[],
            'matched':[],   'distance':[],  'normDistance':[]
        })
        
         
        for mes in mesons:            
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
                           False,
                           999.,
                           999.]
                    
                    df_event.loc[len(df_event)]=row
                    
                    break
        genVtxCouter=genVtxCouter+len(df_event.index)



    # Matching
        SVs     =   []
        for svIdx in range(nSV):
            sv = (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx])
            SVs.append(sv)  
        SVs = np.array(SVs)

        
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
        distances = distances[row_mask][:, col_mask]
        distances_normalized = distances_normalized[row_mask][:, col_mask]

        
        if criterion==1:
        # if one reco/gen vertex has all the possible matching > 2.5 exlude that row or column
            row_ind, col_ind = linear_sum_assignment(distances)
            for i in range(len(row_ind)):
                #matchedDistances.append(distances[row_ind[i], col_ind[i]])
                #normMatchedDistances.append(distances_normalized[row_ind[i], col_ind[i]])
                df_event.loc[col_ind[i], 'matched']=True
                df_event.loc[col_ind[i], 'distance']=distances[row_ind[i], col_ind[i]]
                df_event.loc[col_ind[i], 'normDistance']=distances_normalized[row_ind[i], col_ind[i]]
                
        
        else:
            while (np.any(distances < 997)):
                minIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
                if criterion ==0:
                    if distances[minIdx[0], minIdx[1]]*5 - 2*np.array(df_event.displacement[col_mask])[minIdx[1]] > 1:
                        # dont match the gen part
                        distances[:, minIdx[1]]=[998]*distances.shape[0]
                        continue

                df_event.loc[minIdx[1], 'matched']=True
                df_event.loc[minIdx[1], 'distance']=distances[minIdx[0], minIdx[1]]
                df_event.loc[minIdx[1], 'normDistance']=distances_normalized[minIdx[0], minIdx[1]]
                
                # instead of deleting, replace with high values
                distances[minIdx[0], :] = [998]*distances.shape[1]
                distances[:, minIdx[1]] = [998]*distances.shape[0]

                distances_normalized[minIdx[0], :] = [998]*distances.shape[1]
                distances_normalized[:, minIdx[1]] = [998]*distances.shape[0]
            
        if len(df_event)!=0:
            df = pd.concat((df, df_event))
    
    outName = "/t3home/gcelotto/BTV/output/df_"+suffix+".parquet"
    df.to_parquet(outName)
    
    #bins= np.linspace(0, 2, 100)
    #fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    #y = []
    #hep.cms.label(ax=ax)
    #for i in bins:
    #    y.append((df[df.matched==True].normDistance<i).sum()/(df.matched==True).sum()*100*(df.matched==True).sum()/genVtxCouter)
    #ax.plot(bins, y)
    #ax.text(x=0.98, y=0.85, s="Entries : %d"%((df.matched==True).sum()), ha='right', transform=ax.transAxes)
    #ax.text(x=0.98, y=0.75, s="Eff : %.1f%%"%((df.matched==True).sum()/genVtxCouter*100.), ha='right', transform=ax.transAxes)
    #ax.set_xlabel("dist(SV,genVertex)/dist(SV,PV)")
    #ax.set_ylabel(r"N$_\text{genmatched}$/N$_\text{genVertex}$ [%]")
    #ax.set_xlim(bins[0], bins[-1])
    #ax.set_ylim(20, 100)
    #ax.grid(True)
    #fig.savefig("/t3home/gcelotto/BTV/plots/EffVsThreshold_"+suffix+".png", bbox_inches='tight')
    #print("Saved /t3home/gcelotto/BTV/plots/EffVsThreshold_"+suffix+".png")

    #print((df.matched==True).sum()/genVtxCouter) 
    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    threshold = float(sys.argv[3])
    main(nEvents, criterion, threshold)