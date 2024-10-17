import ROOT
from array import array
import uproot
import numpy as np
from utilsForScript import distance_3d, getPdgMask
import sys
def getSecondMatching():
    secondMatching = False
    return secondMatching
secondMatching = getSecondMatching()
def getMesons(params, GenPart_pdgId_, GenPart_eta_, GenPart_pt_, nGenPart_, GenPart_status_, GenPart_charge, GenPart_genPartIdxMother_):
    pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId_)   # list of pdgId consider as genVertices
    etaMask = abs(GenPart_eta_)<params['max_abs_eta']  # limit the genvertices to GenHadrons within tracker acceptance
    ptMask = GenPart_pt_>params['min_pt_mother']          # require pT of GenHadrons larger than 10 GeV
    mesons = np.arange(nGenPart_)[(pdgMask) & (etaMask) & (ptMask) ] #index of the meosns in nGenPart
    nStableDaughters = np.zeros(len(mesons))
    for idx, mesIdx in enumerate(mesons):
        for genPartIdx in np.arange(nGenPart_)[(GenPart_status_==1) & (abs(GenPart_charge)==1) & (GenPart_pt_>0.8) & (abs(GenPart_eta_)<2.5)]:
            gp = genPartIdx
            while (GenPart_genPartIdxMother_[gp]!=-1):
                if (GenPart_genPartIdxMother_[gp]==mesIdx):
                    nStableDaughters[idx]=nStableDaughters[idx]+1
                    #break the while and look for the next genPart
                    break

                gp = GenPart_genPartIdxMother_[gp]
    mesons = mesons[nStableDaughters>1]
    return mesons

def getParams():
    params = {
            'max_abs_eta': 2.5,
            'min_pt_mother': 10,
        }
    return params

def getOneDaughter(mesons, GenPart_genPartIdxMother_, nGenPart_):
    oneDaughter = []      # one index of a daughter of the meson aligned
    for mes in mesons:
        foundDaughter = -1
        for gp in range(nGenPart_):
            if (GenPart_genPartIdxMother_[gp] == mes):
                oneDaughter.append(gp)
                foundDaughter=gp
                # fill one daughter per meson
                break
        if foundDaughter==-1:
            oneDaughter.append(-1)
    assert len(mesons)==len(oneDaughter)
    return oneDaughter


def matchingEvent(distances, ProbeTracks_matchedToSV_, ProbeTracks_pt_, ProbeTracks_eta_, GenPart_genPartIdxMother_, oneDaughter, ProbeTracks_genPartIdx, nGenPart_):
    distances_filled = distances.copy()
    matchingKey ={}
    while (np.any(distances_filled < 997)):
        #print(distances_filled)
        recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
        #print(recoIdx, genIdx, " possibly matched")
        recoIdx = int(recoIdx)
        genIdx = int(genIdx)

        #tracks from sv
        recoTracksMask = (ProbeTracks_matchedToSV_==recoIdx) & (ProbeTracks_pt_>0.8) & (abs(ProbeTracks_eta_)<2.5)

        # 1. red inefficiency 
#            if len(recoTracksMask)==0:
#                red = red + 1

        # 2. blue inefficiency
#           if len((recoTracksMask) & (ProbeTracks_genPartIdx>-1))==0:
#               blue = blue + 1

        #tracks from gv
        #first mother, grandmother, grandGrandMother
        genTracksMask = (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[oneDaughter[genIdx]])
        genTracksMask = (genTracksMask) | (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]])
        genTracksMask = (genTracksMask) | (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]])
        #genTracksMask = genTracksMask | (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]]])
        #check if genTracks linked to recoTracks are in genVertex

        commonTracks = np.sum(np.in1d(ProbeTracks_genPartIdx[recoTracksMask], np.arange(nGenPart_)[genTracksMask]))
        if commonTracks>=1:
            distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
            matchingKey[genIdx]=recoIdx
        else:
            distances_filled[recoIdx, genIdx]=998
    #print("Distance", distances[recoIdx, genIdx])
        #GenPart_eta_[GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[allDaughters[genIdx]]]
        # if one track reco is similar to a gentrack accept matching

        
    # second matching
    if secondMatching:
        distances_filled = distances.copy()
        while (np.any(distances_filled < 997)):
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)
            skip=False
            if (genIdx in matchingKey.keys()):
                #gen idx already matched so fill the whole column
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                skip=True
            if (recoIdx in matchingKey.values()):
                # same for the recoIdx and row
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
                skip=True
            if skip:
                continue
            if ((distances_filled[recoIdx, genIdx]<0.1)):
                matchingKey[genIdx]=recoIdx
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            else:
                distances_filled[recoIdx, genIdx]=998
    return matchingKey
def main(fileName, fileNumber, prova):
    # inefficiencies
    # red = events where reco SV have no Tracks
    # blue = events where reco SV have no Tracks with genMathcing
    # green = events where reco SV have no Tracks genMatched to daughters of one of the listed hadrons
    # red = 0
    # blue = 0
    # green = 0
    if prova==1:
        fN = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/TTToH_%d.root"%fileNumber
        file = ROOT.TFile(fN, "RECREATE")    
        print(fN)
    elif prova==0:
        file = ROOT.TFile("/scratch/TTToH_%d.root"%fileNumber, "RECREATE")
        print("Created root file in /scratch/TTToH_%d.root"%fileNumber)
    else:
        print("Invalid value for prova")
    tree = ROOT.TTree("tree", "Tree")
    params = getParams()

    event = array('l', [0])
    nSV = array('l', [0])
    #nGenPart = array('d', [0])
    #idx               = array('d', [0])    
    pdgID               = array('l', [0])    
    motherPdgID               = array('l', [0])    
    pt                  = array('d', [0])
    eta                 = array('d', [0])    
    phi                 = array('d', [0])    
    vx                  = array('d', [0])
    vy                  = array('d', [0])
    vz                  = array('d', [0])
    displacement        = array('d', [0])            
    promptStatus        = array('l', [0])   
        # status of D mesons and Charmed Baryons: 
    # 0 prompt (from W)
    # 1,(2),[3],... if B [grand(grand)]mother mother of D
    # -2 neither from W nor from B (from hard scattering e.g.)
    # -1 for all the particles different from D mesons and Charmed Baryons
    decayStatus        = array('l', [0])   
        # status of B mesons and Bottom baryons: 
    # -1 not B nor Bottom baryons
    # 0 does not have D below in the chain
    # 1 mother of D or Charmed baryons
    # 2 grandmother of D or Charmed baryons
    # ...


    B_eta               = array('d', [0])   # eta of the b meson if the D comes from a B.
    B_phi               = array('d', [0])   # phi of the b meson if the D comes from a B.
    nGenTracks          = array('l', [0] )
    genTracks_pt        = array('d', [0.]*15)
    genTracks_eta       = array('d', [0.]*15)
    genTracks_phi       = array('d', [0.]*15)
    genTracks_pdgId     = array('d', [0.]*15)     
    genTracks_charge    = array('l', [0]*15)     
    nRecoTracks         = array('l', [ 0 ] )
    recoTracks_pt       = array('d', [0.]*30)
    recoTracks_eta      = array('d', [0.]*30)                
    recoTracks_phi      = array('d', [0.]*30)                
    distance            = array('d', [0])        
    delta_x            = array('d', [0])        
    delta_y            = array('d', [0])        
    delta_z            = array('d', [0])        
    

    #probeTracksFromSV = array('d', [0])                

    tree.Branch("event",            event,            "event/I")                     
    tree.Branch("nSV",              nSV,                "nSV/I")  
    tree.Branch("pdgID",            pdgID,              "pdgID/I")    
    tree.Branch("motherPdgID",            motherPdgID,              "motherPdgID/I")    
    tree.Branch("pt",               pt,                 "pt/D")
    tree.Branch("eta",              eta,                "eta/D")    
    tree.Branch("phi",              phi,                "phi/D")    
    tree.Branch("vx",               vx,                 "vx/D")
    tree.Branch("vy",               vy,                 "vy/D")
    tree.Branch("vz",               vz,                 "vz/D")
    tree.Branch("displacement",     displacement,       "displacement/D")            
    tree.Branch("promptStatus",     promptStatus,       "promptStatus/I")
    tree.Branch("decayStatus",     decayStatus,       "decayStatus/I")
    tree.Branch("B_eta",     B_eta,       "B_eta/D")
    tree.Branch("B_phi",     B_phi,       "B_phi/D")
    tree.Branch("nGenTracks",  nGenTracks,    "nGenTracks/I")                
    tree.Branch("genTracks_pt",     genTracks_pt,       "genTracks_pt[nGenTracks]/D")            
    tree.Branch("genTracks_eta",     genTracks_eta,       "genTracks_eta[nGenTracks]/D")            
    tree.Branch("genTracks_phi",     genTracks_phi,       "genTracks_phi[nGenTracks]/D")            
    tree.Branch("genTracks_pdgId",  genTracks_pdgId,    "genTracks_pdgId[nGenTracks]/D")                
    tree.Branch("genTracks_charge",  genTracks_charge,    "genTracks_charge[nGenTracks]/I")                

    tree.Branch("nRecoTracks",  nRecoTracks,    "nRecoTracks/I")                
    tree.Branch("recoTracks_pt",     recoTracks_pt,       "recoTracks_pt[nRecoTracks]/D")            
    tree.Branch("recoTracks_eta",  recoTracks_eta,    "recoTracks_eta[nRecoTracks]/D")                
    tree.Branch("recoTracks_phi",  recoTracks_phi,    "recoTracks_phi[nRecoTracks]/D")                
    tree.Branch("distance",         distance,           "distance/D")        
    tree.Branch("delta_x",         delta_x,           "delta_x/D")        
    tree.Branch("delta_y",         delta_y,           "delta_y/D")        
    tree.Branch("delta_z",         delta_z,           "delta_z/D")        

    #tree.Branch("probeTracksFromSV",probeTracksFromSV,  "probeTracksFromSV/D")    


    f = uproot.open(fileName)
    tree_ = f['Events']
    branches = tree_.arrays()
    maxEntries = 100 if prova == 1 else tree_.num_entries
    maxEntries = 500 if maxEntries>10000 else maxEntries

    print("%d entries "%maxEntries)
    for ev in range(maxEntries):
        if ev%100==0:
            print(100 * ev/maxEntries)
        #input("\n\nNext\n\n")
        #if ev == 20:
        #    sys.exit("Exit")
        # input("\n\n next ev %d"%ev)
        # print("event", ev)
        # info extrcted per event
        nSV_                        = branches["nSV"][ev]
        nGenPart_                   = branches["nGenPart"][ev]
        GenPart_genPartIdxMother_   = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId_              = branches["GenPart_pdgId"][ev]
        GenPart_pt_                 = branches["GenPart_pt"][ev]
        GenPart_eta_                = branches["GenPart_eta"][ev]
        GenPart_phi_                = branches["GenPart_phi"][ev]
        GenPart_status_             = branches["GenPart_status"][ev]
        GenPart_charge              = branches["GenPart_charge"][ev]
        PV_x_                       = branches["PV_x"][ev]
        PV_y_                       = branches["PV_y"][ev]
        PV_z_                       = branches["PV_z"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]
        GenPart_vx_                 = branches["GenPart_vx"][ev]
        GenPart_vy_                 = branches["GenPart_vy"][ev]
        GenPart_vz_                 = branches["GenPart_vz"][ev]
        ProbeTracks_matchedToSV_    = branches["ProbeTracks_matchedToSV"][ev]
        ProbeTracks_pt_             = branches["ProbeTracks_pt"][ev]
        ProbeTracks_eta_            = branches["ProbeTracks_eta"][ev]
        ProbeTracks_phi_            = branches["ProbeTracks_phi"][ev]
        ProbeTracks_genPartIdx      = branches["ProbeTracks_genPartIdx"][ev]

        # SVs as reconstructed by IVF
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])

# ***********************************************************************
#                      Gen Vertices definition                         
# ***********************************************************************
        
        mesons = getMesons(params, GenPart_pdgId_, GenPart_eta_, GenPart_pt_, nGenPart_, GenPart_status_, GenPart_charge, GenPart_genPartIdxMother_)     
        oneDaughter = getOneDaughter(mesons, GenPart_genPartIdxMother_, nGenPart_)
        genVertices = np.array([(x, y, z) for x, y, z in zip(GenPart_vx_[oneDaughter], GenPart_vy_[oneDaughter], GenPart_vz_[oneDaughter])])

        #allDaughters = []      # one index of a daughter of the meson aligned
        #for mes in mesons:
        #    for gp in range(nGenPart_):
        #        if (GenPart_genPartIdxMother_[gp] == mes):
        #            allDaughters.append(gp)

# ***********************************************************************
#                            Mathcing                         
# ***********************************************************************
        distances = np.array([[distance_3d(sv, (vx, vy, vz)) for vx, vy, vz in genVertices] for sv in SVs ])
        matchingKey = matchingEvent(distances, ProbeTracks_matchedToSV_, ProbeTracks_pt_, ProbeTracks_eta_, GenPart_genPartIdxMother_, oneDaughter, ProbeTracks_genPartIdx, nGenPart_)
        #input("Next")
        #print(matchingKey)
         

        for genIdx, mes in enumerate(mesons):        
            notNeutrinoMask = (abs(GenPart_pdgId_)!=12) & (abs(GenPart_pdgId_)!=14) & (abs(GenPart_pdgId_)!=16)    
            genDaughtersMask = (GenPart_genPartIdxMother_ == mes) &  (notNeutrinoMask) & (abs(GenPart_eta_)<2.5) & (GenPart_pt_>0.8)
            genTracks_pt_ = GenPart_pt_[genDaughtersMask]
            genTracks_eta_ = GenPart_eta_[genDaughtersMask]
            genTracks_phi_ = GenPart_phi_[genDaughtersMask]
            genTracks_pdgId_ = GenPart_pdgId_[genDaughtersMask]
            genTracks_charge_ = GenPart_charge[genDaughtersMask]
            assert len(genDaughtersMask)==len(GenPart_pdgId_)
            #if np.sum(genDaughtersMask)>4:
                #print(GenPart_pdgId_[genDaughtersMask])
            #print(genTracks_pt_)
            if genIdx in matchingKey:
                recoTracksMask = (abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx]) & (ProbeTracks_pt_>0.8)
                nRecoTracks[0]=np.sum(recoTracksMask)
                recoTracks_pt_ = ProbeTracks_pt_[recoTracksMask]
                recoTracks_eta_ = ProbeTracks_eta_[recoTracksMask]
                recoTracks_phi_ = ProbeTracks_phi_[recoTracksMask]
            else:
                nRecoTracks[0]=0
            for i in range(nRecoTracks[0]):
                recoTracks_pt[i] = recoTracks_pt_[i]
                recoTracks_eta[i] = recoTracks_eta_[i]
                recoTracks_phi[i] = recoTracks_phi_[i]


            if genIdx in matchingKey:
                distance[0]=distances[matchingKey[genIdx], genIdx]
                delta_x[0] = SV_x[matchingKey[genIdx]] - GenPart_vx_[oneDaughter[genIdx]]
                delta_y[0] = SV_y[matchingKey[genIdx]] - GenPart_vy_[oneDaughter[genIdx]]
                delta_z[0] = SV_z[matchingKey[genIdx]] - GenPart_vz_[oneDaughter[genIdx]]
            else:
                distance[0]=-1
                delta_x[0] = -1
                delta_y[0] = -1
                delta_z[0] = -1
            for gp in range(nGenPart_):
                if (GenPart_genPartIdxMother_[gp] == mes):
                    # found the decay of the meson
                    displacement[0]=distance_3d((GenPart_vx_[gp], GenPart_vy_[gp], GenPart_vz_[gp]), (GenPart_vx_[mes], GenPart_vy_[mes], GenPart_vz_[mes]))
                    break

            # meson variables
            event[0]    = ev
            nSV[0]      = nSV_
            pdgID[0]    = GenPart_pdgId_[mes]
            motherPdgID[0]    = GenPart_pdgId_[GenPart_genPartIdxMother_[mes]]
            promptStatus[0] = -1 # not a charmed
            decayStatus[0] = -1 # not a bottom
            B_eta[0]        = -999
            B_phi[0]        = -999
            bottomList = [511, 521, 531,
                            5122, 5132, 5232, 5332,
                            5142, 5242, 5342, 5512, 5532, 5542, 5554,]
            charmedList = [411, 421, 431,   
                           4122,4232,4132,4332,4412,4422,4432,4444,]
            
        # Check if in the tree (on top) there is a B that makes the D or CB not prompt
            if (abs(GenPart_pdgId_[mes]) in charmedList):
                #input("charmed found %d"%(GenPart_pdgId_[mes]))
                promptStatus[0] = -2 # neither W nor B in the upper chain
                daughterToBeChecked = mes
                counterOfDegree = 0
                while ((promptStatus[0] == -2) & (GenPart_genPartIdxMother_[daughterToBeChecked]!=-1)):  # continue from here
                    counterOfDegree = counterOfDegree + 1
                    #print("Mother : ", GenPart_pdgId_[GenPart_genPartIdxMother_[daughterToBeChecked]])
                    if abs(GenPart_pdgId_[GenPart_genPartIdxMother_[daughterToBeChecked]]) == 24:
                        promptStatus[0] = 0
                    elif abs(GenPart_pdgId_[GenPart_genPartIdxMother_[daughterToBeChecked]]) in bottomList:
                        promptStatus[0] = counterOfDegree
                        B_eta[0]=GenPart_eta_[GenPart_genPartIdxMother_[daughterToBeChecked]]
                        B_phi[0]=GenPart_phi_[GenPart_genPartIdxMother_[daughterToBeChecked]]
                        break
                    else:
                        daughterToBeChecked = GenPart_genPartIdxMother_[daughterToBeChecked]
            
        # Check if in the tree (down) there is a D 
            if (abs(GenPart_pdgId_[mes]) in bottomList):
                decayStatus[0] = 0 # does not have a D in the chain

                for charmed in np.arange(nGenPart_)[(np.isin(abs(GenPart_pdgId_),charmedList))]:
                    # you should also consider the case in which you have 2 charmed from the B. then which one you keep?
                    # in this case it is the first one in the list
                    if charmed not in mesons:
                        continue
                    daughterToBeChecked = charmed
                    counterOfDegree = 0
                    while ((decayStatus[0]==0) & (GenPart_genPartIdxMother_[daughterToBeChecked]!=-1)):  # continue from here
                        counterOfDegree = counterOfDegree + 1
                        if GenPart_genPartIdxMother_[daughterToBeChecked] == mes:
                            decayStatus[0] = counterOfDegree
                            break
                        else:
                            daughterToBeChecked = GenPart_genPartIdxMother_[daughterToBeChecked]        

            pt[0]       = GenPart_pt_[mes]
            eta[0]      = GenPart_eta_[mes]
            phi[0]      = GenPart_phi_[mes]
            vx[0]       = GenPart_vx_[gp]
            vy[0]       = GenPart_vy_[gp]
            vz[0]       = GenPart_vz_[gp]
            nGenTracks[0]=len(genTracks_pt_)
            assert len(genTracks_pt_)==np.sum(genDaughtersMask)
            for dau in range(nGenTracks[0]):
                genTracks_pt[dau] = (genTracks_pt_[dau])
                genTracks_eta[dau] = (genTracks_eta_[dau])
                genTracks_phi[dau] = (genTracks_phi_[dau])
                genTracks_pdgId[dau] = float(genTracks_pdgId_[dau])
                genTracks_charge[dau] = int(genTracks_charge_[dau])

            tree.Fill()
    tree.Write()
    file.Close()
    #print("red", red)
    #print("blue", blue)


if __name__ == "__main__":
    fileName = sys.argv[1] if len(sys.argv)>1 else '/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/240723_073316/0000/TTToHadronic_Run2_mc_2024Jul23_100.root'
    fileNumber = int(sys.argv[2]) if len(sys.argv)>1 else 100
    prova = 0 if len(sys.argv)>1 else 1
    main(fileName, fileNumber, prova)