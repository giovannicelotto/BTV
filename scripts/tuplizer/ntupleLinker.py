import ROOT
from array import array
import uproot
import numpy as np
from utilsForScript import distance_3d, getPdgMask
import sys
def main(fileName, fileNumber, prova):
    # inefficiencies
    # red = events where reco SV have no Tracks
    # blue = events where reco SV have no Tracks with genMathcing
    # green = events where reco SV have no Tracks genMatched to daughters of one of the listed hadrons
    # red = 0
    # blue = 0
    # green = 0
    if prova==1:
        file = ROOT.TFile("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul10Tuple/TTToH_%d.root"%fileNumber, "RECREATE")    
        print("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06Tuple/TTToH_%d.root"%fileNumber)
    elif prova==0:
        file = ROOT.TFile("/scratch/TTToH_%d.root"%fileNumber, "RECREATE")
        print("Created root file in /scratch/TTToH_%d.root"%fileNumber)
    else:
        print("Invalid value for prova")
    tree = ROOT.TTree("tree", "Tree")
    params = {
            'max_abs_eta': 2.5,
            'min_pt_mother': 10,
        }

    event = array('l', [0])
    nSV = array('l', [0])
    #nGenPart = array('d', [0])
    #idx               = array('d', [0])    
    pdgID               = array('l', [0])    
    pt                  = array('d', [0])
    eta                 = array('d', [0])    
    phi                 = array('d', [0])    
    vx                  = array('d', [0])
    vy                  = array('d', [0])
    vz                  = array('d', [0])
    displacement        = array('d', [0])            
    nGenTracks          = array('l', [0] )
    genTracks_pt        = array('d', [0.]*15)
    genTracks_eta       = array('d', [0.]*15)
    genTracks_phi       = array('d', [0.]*15)
    genTracks_pdgId     = array('d', [0.]*15)     
    nRecoTracks         = array('l', [ 0 ] )
    recoTracks_pt       = array('d', [0.]*30)
    recoTracks_eta      = array('d', [0.]*30)                
    recoTracks_phi      = array('d', [0.]*30)                
    distance            = array('d', [0])        

    #probeTracksFromSV = array('d', [0])                

    tree.Branch("event",            event,            "event/I")                     
    tree.Branch("nSV",              nSV,                "nSV/I")  
    tree.Branch("pdgID",            pdgID,              "pdgID/I")    
    tree.Branch("pt",               pt,                 "pt/D")
    tree.Branch("eta",              eta,                "eta/D")    
    tree.Branch("phi",              phi,                "phi/D")    
    tree.Branch("vx",               vx,                 "vx/D")
    tree.Branch("vy",               vy,                 "vy/D")
    tree.Branch("vz",               vz,                 "vz/D")
    tree.Branch("displacement",     displacement,       "displacement/D")            
    tree.Branch("nGenTracks",  nGenTracks,    "nGenTracks/I")                
    tree.Branch("genTracks_pt",     genTracks_pt,       "genTracks_pt[nGenTracks]/D")            
    tree.Branch("genTracks_eta",     genTracks_eta,       "genTracks_eta[nGenTracks]/D")            
    tree.Branch("genTracks_phi",     genTracks_phi,       "genTracks_phi[nGenTracks]/D")            
    tree.Branch("genTracks_pdgId",  genTracks_pdgId,    "genTracks_pdgId[nGenTracks]/D")                

    tree.Branch("nRecoTracks",  nRecoTracks,    "nRecoTracks/I")                
    tree.Branch("recoTracks_pt",     recoTracks_pt,       "recoTracks_pt[nRecoTracks]/D")            
    tree.Branch("recoTracks_eta",  recoTracks_eta,    "recoTracks_eta[nRecoTracks]/D")                
    tree.Branch("recoTracks_phi",  recoTracks_phi,    "recoTracks_phi[nRecoTracks]/D")                
    tree.Branch("distance",         distance,           "distance/D")        

    #tree.Branch("probeTracksFromSV",probeTracksFromSV,  "probeTracksFromSV/D")    


    f = uproot.open(fileName)
    tree_ = f['Events']
    branches = tree_.arrays()
    maxEntries = 500 if prova == 1 else tree_.num_entries
    maxEntries = 500 if maxEntries>10000 else maxEntries

    print("%d entries "%maxEntries)
    for ev in range(maxEntries):
        if ev%100==0:
            print(100 * ev/maxEntries)
        input("\n\nNext\n\n")
        if ev == 3:
            sys.exit("Exit")
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
        
        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId_)   # list of pdgId consider as genVertices
        etaMask = abs(GenPart_eta_)<params['max_abs_eta']  # limit the genvertices to GenHadrons within tracker acceptance
        ptMask = GenPart_pt_>params['min_pt_mother']          # require pT of GenHadrons larger than 10 GeV
        # look for genparticles that has the meson as mother and stable
        # loop in mesons taken so far
        # for every genparticles check all the mothers if there is the particular meson
        # count the daughters
        # 
        mesons = np.arange(nGenPart_)[(pdgMask) & (etaMask) & (ptMask) ] #index of the meosns in nGenPart
        nStableDaughters = np.zeros(len(mesons))
        for idx, mesIdx in enumerate(mesons):
            for genPartIdx in np.arange(nGenPart_)[(GenPart_status_==1) & (abs(GenPart_charge)==1) & (GenPart_pt_>1) & (abs(GenPart_eta_)<2.5)]:
                gp = genPartIdx
                while (GenPart_genPartIdxMother_[gp]!=-1):
                    if (GenPart_genPartIdxMother_[gp]==mesIdx):
                        nStableDaughters[idx]=nStableDaughters[idx]+1
                        #break the while and look for the next genPart
                        break

                    gp = GenPart_genPartIdxMother_[gp]
        mesons = mesons[nStableDaughters>1]
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
        print(distances)
        
        distances_filled = distances.copy()
        matchingKey ={}
        while (np.any(distances_filled < 997)):
            print(distances_filled)
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            print(recoIdx, genIdx, " possibly matched")
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)

            #tracks from sv
            recoTracksMask = (ProbeTracks_matchedToSV_==recoIdx) & (ProbeTracks_pt_>1) & (abs(ProbeTracks_eta_)<2.5)

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
        #distances_filled = distances.copy()
        #while (np.any(distances_filled < 997)):
        #    recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
        #    recoIdx = int(recoIdx)
        #    genIdx = int(genIdx)
        #    skip=False
        #    if (genIdx in matchingKey.keys()):
        #        #gen idx already matched so fill the whole column
        #        distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
        #        skip=True
        #    if (recoIdx in matchingKey.values()):
        #        # same for the recoIdx and row
        #        distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
        #        skip=True
        #    if skip:
        #        continue
        #    if ((distances_filled[recoIdx, genIdx]<0.1)):
        #        matchingKey[genIdx]=recoIdx
        #        distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
        #        distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
        #    else:
        #        distances_filled[recoIdx, genIdx]=998



        for genIdx, mes in enumerate(mesons):        
            notNeutrinoMask = (abs(GenPart_pdgId_)!=12) & (abs(GenPart_pdgId_)!=14) & (abs(GenPart_pdgId_)!=16)    
            genDaughtersMask = (GenPart_genPartIdxMother_ == mes) &  (notNeutrinoMask) & (abs(GenPart_eta_)<2.5) & (GenPart_pt_>1)
            genTracks_pt_ = GenPart_pt_[genDaughtersMask]
            genTracks_eta_ = GenPart_eta_[genDaughtersMask]
            genTracks_phi_ = GenPart_phi_[genDaughtersMask]
            genTracks_pdgId_ = GenPart_pdgId_[genDaughtersMask]
            assert len(genDaughtersMask)==len(GenPart_pdgId_)
            #if np.sum(genDaughtersMask)>4:
                #print(GenPart_pdgId_[genDaughtersMask])
            #print(genTracks_pt_)
            if genIdx in matchingKey:
                recoTracksMask = (abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx]) & (ProbeTracks_pt_>1)
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
            else:
                distance[0]=-1
            for gp in range(nGenPart_):
                if (GenPart_genPartIdxMother_[gp] == mes):
                    # found the decay of the meson
                    displacement[0]=distance_3d((GenPart_vx_[gp], GenPart_vy_[gp], GenPart_vz_[gp]), (GenPart_vx_[mes], GenPart_vy_[mes], GenPart_vz_[mes]))
                    break
                        #genTracks_pdgId=genTracks_pdgId_
            # meson variables
            event[0]    = ev
            nSV[0]      = nSV_
            pdgID[0]    = GenPart_pdgId_[mes]
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

            tree.Fill()
    tree.Write()
    file.Close()
    #print("red", red)
    #print("blue", blue)


if __name__ == "__main__":
    fileName = sys.argv[1] if len(sys.argv)>1 else '/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul10/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/240710_125611/0000/TTToHadronic_Run2_mc_2024Jul10_100.root'
    fileNumber = int(sys.argv[2]) if len(sys.argv)>1 else 100
    prova = 0 if len(sys.argv)>1 else 1
    main(fileName, fileNumber, prova)