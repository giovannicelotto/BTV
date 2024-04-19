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

def eventDisplay(df_event, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy):
    fig, ax = plt.subplots(1, 1)
    
    ax.scatter(df_event.vx, df_event.vy, label='GenVertices')

    ax.text(x=1.02, y=0.45, s="x", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.08, y=0.45, s="y", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.14, y=0.45, s="z", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.2, y=0.45, s="PDG", ha='center', fontsize=14, transform=ax.transAxes)
    y=0.4
    for row in range(len(df_event)):
        ax.text(x=df_event.vx.iloc[row], y=df_event.vy.iloc[row], s=df_event.pdgClass.iloc[row], fontsize=12)
        ax.text(x=1.02, y=y, s="%.1f"%(df_event.vx.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.08, y=y, s="%.1f"%(df_event.vy.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.14, y=y, s="%.1f"%(df_event.vz.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.2, y=y, s="%s"%(df_event.pdgClass[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        y=y-0.05
        if GenPart_genPartIdxMother[df_event.idx.iloc[row]]!=-1:
            ax.plot([df_event.vx.iloc[row], GenPart_vx[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]], [df_event.vy.iloc[row], GenPart_vy[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]],linestyle='dotted',alpha=0.8, marker='none')
    if SVs.shape[0]>0:
        ax.text(x=1.02, y=0.95, s="x", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.08, y=0.95, s="y", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.14, y=0.95, s="z", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.2, y=0.95, s="chi2", ha='center', fontsize=14, transform=ax.transAxes)
        y=0.9
        for row in range(len(SVs)):
            ax.text(x=1.02, y=y, s="%.1f"%(SVs[row][0]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.08, y=y, s="%.1f"%(SVs[row][1]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.14, y=y, s="%.1f"%(SVs[row][2]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.20, y=y, s="%.1f"%(SV_chi2[row]), ha='center', fontsize=14, transform=ax.transAxes)
            y=y-0.05

        ax.scatter(SVs[:,0], SVs[:,1], label='Reco SV', marker='^', color='green')
    else:
        pass
        assert len(SVs)==0
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")

    ax.scatter(PV_x, PV_y, label='PV', color='red')


    theta = np.linspace(0, 2*np.pi, 100)
    x = PV_x + 0.1 * np.cos(theta)
    y = PV_y + 0.1 * np.sin(theta)
    ax.plot(x, y, color='C0')
    # mathcing
    ax.legend()
    return fig, ax

def criterion0(distances, distances_normalized, df_event, display, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy, col_mask, row_mask):
    if display:
        fig, ax = eventDisplay(df_event=df_event, SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx, GenPart_vy=GenPart_vy)
    while (np.any(distances < 997)):
        minIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
        #try:

        if distances[minIdx[0], minIdx[1]]*10 - np.array(df_event.displacement[col_mask])[minIdx[1]] > 5:
            # dont match the gen part
            distances[minIdx[0], minIdx[1]]=998
            continue

        df_event.loc[minIdx[1], 'matched']=True
        df_event.loc[minIdx[1], 'distance']=distances[minIdx[0], minIdx[1]]
        df_event.loc[minIdx[1], 'normDistance']=distances_normalized[minIdx[0], minIdx[1]]
        x_values = [SVs[row_mask][minIdx[0]][0], np.array(df_event.vx)[col_mask][minIdx[1]]]
        y_values = [SVs[row_mask][minIdx[0]][1], np.array(df_event.vy)[col_mask][minIdx[1]]]

        #plot matching
        ax.plot(x_values, y_values, color='black', marker='none')

        # instead of deleting matched vertices, replace them with high values
        distances[minIdx[0], :]=[998]*distances.shape[1]
        distances[:, minIdx[1]]=[998]*distances.shape[0]

        distances_normalized[minIdx[0], :] = [998]*distances.shape[1]
        distances_normalized[:, minIdx[1]] = [998]*distances.shape[0]
    if display:
        fig.savefig("/t3home/gcelotto/BTV/plots/EventGenReco.png", bbox_inches='tight')
        plt.close('all')
    return df_event

def criterion1(distances, distances_normalized, df_event, display, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy, col_mask, row_mask):
    if display:
        fig, ax = eventDisplay(df_event=df_event, SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx, GenPart_vy=GenPart_vy)
    row_ind, col_ind = linear_sum_assignment(distances)
    for i in range(len(row_ind)):
        df_event.loc[col_ind[i], 'matched']=True
        df_event.loc[col_ind[i], 'distance']=distances[row_ind[i], col_ind[i]]
        df_event.loc[col_ind[i], 'normDistance']=distances_normalized[row_ind[i], col_ind[i]]
        x_values = [SVs[row_mask][row_ind[i]][0], np.array(df_event.vx)[col_mask][col_ind[i]]]
        y_values = [SVs[row_mask][row_ind[i]][1], np.array(df_event.vy)[col_mask][col_ind[i]]]
        ax.plot(x_values, y_values, color='black', marker='none')
    return df_event

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
        SV_chi2                     = branches["SV_chi2"][ev]
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
        print(distances)

        
        if criterion:
            df_event = criterion1(distances=distances, distances_normalized=distances_normalized, df_event=df_event, disaply=True, SVs=SVs,
                                  PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  GenPart_vy=GenPart_vy, col_mask=col_mask, row_mask=row_mask)
        else:
            df_event = criterion0(distances=distances, distances_normalized=distances_normalized, df_event=df_event, display=True,
                                  SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx,
                                  GenPart_vy=GenPart_vy, col_mask=col_mask, row_mask=row_mask)

        
        
        
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


    print((df.matched==True).sum()/genVtxCouter) 
    return 

if __name__ =="__main__":
    nEvents = int(sys.argv[1]) if len(sys.argv)>1 else None
    criterion = int(sys.argv[2])
    threshold = float(sys.argv[3])
    print(criterion)
    main(nEvents, criterion, threshold)