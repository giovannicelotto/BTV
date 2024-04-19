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
from utilsForScript import distance_3d, getPdgMask
#from efficiency_differential import map_to_groups
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
    
    f = uproot.open("/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root")
    tree = f['Events']
    branches = tree.arrays()
    GenPart_pdgId               = branches["GenPart_pdgId"]
    
    
    genVtxCouter = 0
    #matchedDistances = []
    #normMatchedDistances = []
    df = pd.DataFrame({
            'index'     :[],
            'pdgID'     :[],
            'pt'        :[],
            'eta'       :[],
            'phi'       :[],
            'vx'        :[],
            'vy'        :[],
            'vz'        :[],
            'displacement':[],
            'matched'   :[],
            'distance'  :[]
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
        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]

        # filter the gen part of interest
        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId)
        etaMask = abs(GenPart_eta)<params['eta_threshold']
        ptMask = GenPart_pt>params['pt_threshold']
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask]
        
        df_event = pd.DataFrame({
            'idx':[],
            'pdgID':[],
            'pt':[],
            'eta':[],
            'phi':[],
            'vx':[],
            'vy':[],
            'vz':[],
            'displacement':[],
            'matched':[],
            'distance':[],
            'normDistance':[]
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
                           distance_3d(
                                        (GenPart_vx[gp], GenPart_vy[gp], GenPart_vz[gp]),
                                        (GenPart_vx[mes], GenPart_vy[mes], GenPart_vz[mes])),
                           False,
                           999.,
                           999.]
                    
                    df_event.loc[len(df_event)]=row
                    
                    break
        genVtxCouter=genVtxCouter+len(df_event.idx)
        df_event['pdgClass'] = df_event['pdgID'].apply(map_to_groups_letter)



    # Matching
        SVs     =   []
        for svIdx in range(nSV):
            sv = (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx])
            SVs.append(sv)  
        SVs = np.array(SVs)
        
        input("next")
        # row by row plots
        if ev==0:
            fig, ax = plt.subplots(1, 1)
        else:
            ax.clear()
        ax.scatter(df_event.vx, df_event.vy, label='GenVertices')
        
        for row in range(len(df_event)):
            ax.text(x=df_event.vx.iloc[row], y=df_event.vy.iloc[row], s=df_event.pdgClass.iloc[row], fontsize=12)
            if GenPart_genPartIdxMother[df_event.idx.iloc[row]]!=-1:
                ax.plot([df_event.vx.iloc[row], GenPart_vx[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]], [df_event.vy.iloc[row], GenPart_vy[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]],linestyle='dotted',alpha=0.8, marker='none')
        if SVs.shape[0]>0:
            ax.scatter(SVs[:,0], SVs[:,1], label='Reco SV', marker='x')
        else:
            pass
            assert len(SVs)==0
        ax.set_xlabel("X [cm]")
        ax.set_ylabel("Y [cm]")
        
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
        print("Shape before ", distances.shape)
        distances = distances[row_mask][:, col_mask]
        print("Shape after ", distances.shape)
        distances_normalized = distances_normalized[row_mask][:, col_mask]
        #print(distances_normalized)
        print(distances)

        
        if criterion:
        # if one reco/gen vertex has all the possible matching > 2.5 exlude that row or column
            row_ind, col_ind = linear_sum_assignment(distances)
            for i in range(len(row_ind)):
                #matchedDistances.append(distances[row_ind[i], col_ind[i]])
                #normMatchedDistances.append(distances_normalized[row_ind[i], col_ind[i]])
                df_event.loc[col_ind[i], 'matched']=True
                df_event.loc[col_ind[i], 'distance']=distances[row_ind[i], col_ind[i]]
                df_event.loc[col_ind[i], 'normDistance']=distances_normalized[row_ind[i], col_ind[i]]
                x_values = [SVs[row_mask][row_ind[i]][0], np.array(df_event.vx)[col_mask][col_ind[i]]]
                y_values = [SVs[row_mask][row_ind[i]][1], np.array(df_event.vy)[col_mask][col_ind[i]]]
                ax.plot(x_values, y_values, color='black', marker='none')
        
        else:
            istances = 0
            while (np.any(distances < 997)):
                minIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
                #try:
                    
                if distances[minIdx[0], minIdx[1]]*10 - np.array(df_event.displacement[col_mask])[minIdx[1]] > 5:
                    
                        # dont match the gen part
                    distances[minIdx[0], minIdx[1]]=998
                    if istances==0:
                        ax.scatter(np.array(df_event.vx[col_mask])[minIdx[1]], np.array(df_event.vy[col_mask])[minIdx[1]], label='GenVertices removed', color='pink')
                    else:
                        ax.scatter(np.array(df_event.vx[col_mask])[minIdx[1]], np.array(df_event.vy[col_mask])[minIdx[1]],  color='pink')
                    istances=istances+1
                    continue
                #except:
                #    print("Broken")
                #    sys.exit("Exit")
                    
                    
                #matchedDistances.append(distances[minIdx[0], minIdx[1]])
                #normMatchedDistances.append(distances_normalized[minIdx[0], minIdx[1]])
                df_event.loc[minIdx[1], 'matched']=True
                df_event.loc[minIdx[1], 'distance']=distances[minIdx[0], minIdx[1]]
                df_event.loc[minIdx[1], 'normDistance']=distances_normalized[minIdx[0], minIdx[1]]
                x_values = [SVs[row_mask][minIdx[0]][0], np.array(df_event.vx)[col_mask][minIdx[1]]]
                y_values = [SVs[row_mask][minIdx[0]][1], np.array(df_event.vy)[col_mask][minIdx[1]]]
                ax.plot(x_values, y_values, color='black', marker='none')
                
                # instead of deleting, replace with high values
                distances[minIdx[0], :]=[998]*distances.shape[1]
                distances[:, minIdx[1]]=[998]*distances.shape[0]
                

                
                distances_normalized[minIdx[0], :] = [998]*distances.shape[1]
                distances_normalized[:, minIdx[1]] = [998]*distances.shape[0]


        ax.scatter(PV_x, PV_y, label='PV', color='red')
        theta = np.linspace(0, 2*np.pi, 100)
        x = PV_x + 0.1 * np.cos(theta)
        y = PV_y + 0.1 * np.sin(theta)
        ax.plot(x, y, color='C0')

        ax.legend()
        #print(distances_normalized)
        fig.savefig("/t3home/gcelotto/BTV/plots/EventGenReco.png", bbox_inches='tight')
        #print(df_event.distance)
        if len(df_event)!=0:
            df = pd.concat((df, df_event))

        
    import math
    floatPart, intPart = math.modf(params['threshold'])[0], int(math.modf(params['threshold'])[1])
    floatPart = int(floatPart*10)
    if floatPart == 0:
        suffix = "cr%d_t%d"%(criterion, intPart)
    else:
        suffix = "cr%d_t%dp%d"%(criterion, intPart, floatPart)
    outName = "/t3home/gcelotto/BTV/output/df_"+suffix+".parquet"
    df.to_parquet(outName)
    

# Plots
    #fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    bins= np.linspace(0, 2, 100)
    #hep.cms.label(ax=ax)
    #ax.hist(np.clip(matchedDistances, bins[0], bins[-1]), bins=bins)
    #ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(matchedDistances)), ha='right', transform=ax.transAxes)
    #if params['threshold']<999:
    #    print("Adding text")
    #    ax.text(x=0.98, y=.75, s="Threshold : %.1f"%(params['threshold']), ha='right', transform=ax.transAxes)
    #ax.set_xlabel("GenMatched 3d distance [cm]", fontsize=15)
    #ax.set_yscale('log')
    #fig.savefig("/t3home/gcelotto/BTV/plots/genMatched_distance_cr1.png", bbox_inches='tight')
    #plt.close('all')


    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    thresholds = bins
    y = []
    hep.cms.label(ax=ax)
    for i in thresholds:
        y.append((df[df.matched==True].normDistance<i).sum()/(df.matched==True).sum()*100*(df.matched==True).sum()/genVtxCouter)
    #if params['threshold']==1:
    ax.plot(thresholds, y)
    ax.text(x=0.98, y=0.85, s="Entries : %d"%((df.matched==True).sum()), ha='right', transform=ax.transAxes)
    ax.text(x=0.98, y=0.75, s="Eff : %.1f%%"%((df.matched==True).sum()/genVtxCouter*100.), ha='right', transform=ax.transAxes)

    #else:
    #    ax.plot(thresholds, y)
    #    ax.text(x=0.98, y=0.85, s="Entries : %d"%(len(normMatchedDistances)), ha='right', transform=ax.transAxes)
    #    ax.text(x=0.98, y=0.75, s="Eff : %.1f %%"%(len(normMatchedDistances)/genVtxCouter*100), ha='right', transform=ax.transAxes)
    
    ax.set_xlabel("dist(SV,genVertex)/dist(SV,PV)")
    ax.set_ylabel(r"N$_\text{genmatched}$/N$_\text{genVertex}$ [%]")
    ax.set_xlim(bins[0], bins[-1])#thresholds[-1])
    ax.set_ylim(20, 70)
    ax.grid(True)

    fig.savefig("/t3home/gcelotto/BTV/plots/EffVsThreshold_"+suffix+".png", bbox_inches='tight')
    print("Saved /t3home/gcelotto/BTV/plots/EffVsThreshold_"+suffix+".png")



    print((df.matched==True).sum()/genVtxCouter) 
    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    threshold = float(sys.argv[3])
    print(criterion)
    main(nEvents, criterion, threshold)