import ROOT
from array import array
import uproot
import numpy as np
from utilsForScript import distance_3d, getPdgMask
import sys
def main(fileName, fileNumber):
    prova = 0
    if prova==0:
        file = ROOT.TFile("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06Tuple/TTToH_%d.root"%fileNumber, "RECREATE")    
        print("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06Tuple/TTToH_%d.root"%fileNumber)
    else:
        file = ROOT.TFile("/scratch/TTToH_%d.root"%fileNumber, "RECREATE")
        print("Created root file in /scratch/TTToH_%d.root"%fileNumber)
    tree = ROOT.TTree("tree", "Tree")
    params = {
            'max_abs_eta': 2.5,
            'min_pt_mother': 10,
        }

    event = array('l', [0])
    nSV = array('l', [0])
    #nGenPart = array('d', [0])
    #idx               = array('d', [0])    
    pdgID             = array('l', [0])    
    pt                = array('d', [0])
    eta               = array('d', [0])    
    phi               = array('d', [0])    
    vx                = array('d', [0])
    vy                = array('d', [0])
    vz                = array('d', [0])
    displacement      = array('d', [0])            
    nGenTracks       = array('l', [ 0 ] )
    genTracks_pt        = array('d')
    genTracks_pdgId     = array('l')     
    nRecoTracks       = array('l', [ 0 ] )
    recoTracks_pt        = array('d')
    recoTracks_eta     = array('d')                
    recoTracks_phi     = array('d')                
    distance          = array('d', [0])        

    #probeTracksFromSV = array('d', [0])                

    tree.Branch("event",              event,                "event/I")                     
    tree.Branch("nSV",              nSV,                "nSV/I")
    #tree.Branch("idx",              idx,                "idx/D")    
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
    tree.Branch("genTracks_pdgId",  genTracks_pdgId,    "genTracks_pdgId[nGenTracks]/I")                

    tree.Branch("nRecoTracks",  nRecoTracks,    "nRecoTracks/I")                
    tree.Branch("recoTracks_pt",     recoTracks_pt,       "recoTracks_pt[nRecoTracks]/D")            
    tree.Branch("recoTracks_eta",  recoTracks_eta,    "recoTracks_eta[nRecoTracks]/D")                
    tree.Branch("recoTracks_phi",  recoTracks_phi,    "recoTracks_phi[nRecoTracks]/D")                
    tree.Branch("distance",         distance,           "distance/D")        

    #tree.Branch("probeTracksFromSV",probeTracksFromSV,  "probeTracksFromSV/D")    


    f = uproot.open(fileName)
    tree_ = f['Events']
    branches = tree_.arrays()
    maxEntries = 1000 if prova == 0 else tree_.num_entries
    print("%d entries "%maxEntries)
    for ev in range(maxEntries):
        #if ev%100==0:
        #    print(100 * ev/maxEntries)
        #input("next")
        print("event", ev)
        # info extrcted per event
        nSV_                         = branches["nSV"][ev]
        nGenPart_                    = branches["nGenPart"][ev]
        GenPart_genPartIdxMother_    = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId_               = branches["GenPart_pdgId"][ev]
        GenPart_pt_                  = branches["GenPart_pt"][ev]
        GenPart_eta_                 = branches["GenPart_eta"][ev]
        GenPart_phi_                 = branches["GenPart_phi"][ev]
        PV_x_                        = branches["PV_x"][ev]
        PV_y_                        = branches["PV_y"][ev]
        PV_z_                        = branches["PV_z"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]
        GenPart_vx_                  = branches["GenPart_vx"][ev]
        GenPart_vy_                  = branches["GenPart_vy"][ev]
        GenPart_vz_                  = branches["GenPart_vz"][ev]
        ProbeTracks_matchedToSV_     = branches["ProbeTracks_matchedToSV"][ev]
        ProbeTracks_pt_              = branches["ProbeTracks_pt"][ev]
        ProbeTracks_eta_             = branches["ProbeTracks_eta"][ev]
        ProbeTracks_phi_             = branches["ProbeTracks_phi"][ev]

        # SVs anf genVertex
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])

        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId_)   # list of pdgId consider as genVertices
        etaMask = abs(GenPart_eta_)<params['max_abs_eta']  # limit the genvertices to GenHadrons within tracker acceptance
        ptMask = GenPart_pt_>params['min_pt_mother']          # require pT of GenHadrons larger than 10 GeV
        mesons = np.arange(nGenPart_)[pdgMask & etaMask & ptMask ] #index of the meosns in nGenPart
        oneDaughter = []      # one index of a daughter of the meson aligned
        for mes in mesons:
            for gp in range(nGenPart_):
                if (GenPart_genPartIdxMother_[gp] == mes):
                    oneDaughter.append(gp)
                    # fill one daughter per meson
                    break
        #print(oneDaughter)
        genVertices = np.array([(x, y, z) for x, y, z in zip(GenPart_vx_[oneDaughter], GenPart_vy_[oneDaughter], GenPart_vz_[oneDaughter])])
        allDaughters = []      # one index of a daughter of the meson aligned
        for mes in mesons:
            for gp in range(nGenPart_):
                if (GenPart_genPartIdxMother_[gp] == mes):
                    allDaughters.append(gp)


        distances = np.array([[distance_3d(sv, (vx, vy, vz)) for vx, vy, vz in genVertices] for sv in SVs ])
        distances_filled = distances.copy()
        matchingKey ={}
        while (np.any(distances_filled < 997)):
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            #print(recoIdx, genIdx, " possibly matched")
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)
            #print("Distance", distances[recoIdx, genIdx])
            #GenPart_eta_[GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[allDaughters[genIdx]]]
            # if one track reco is similar to a gentrack accept matching

            recoMask = (ProbeTracks_pt_>1) & (abs(ProbeTracks_eta_)<2.5) 
            genMask = (GenPart_pt_>1) & (abs(GenPart_eta_)<2.5)

            eta_recos = ProbeTracks_eta_[abs(ProbeTracks_matchedToSV_==recoIdx) & (recoMask)]
            phi_recos = ProbeTracks_phi_[abs(ProbeTracks_matchedToSV_==recoIdx) & (recoMask)]
            pt_recos = ProbeTracks_pt_[abs(ProbeTracks_matchedToSV_==recoIdx) & (recoMask)]
            
            eta_gens = GenPart_eta_[(GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[allDaughters[genIdx]]) & genMask]
            phi_gens = GenPart_phi_[(GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[allDaughters[genIdx]]) & genMask]
            pt_gens = GenPart_pt_[(GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[allDaughters[genIdx]]) & genMask]
            #print("len(eta_recos) ",len(eta_recos))
            if len(eta_recos)>0:
                deltaRMatrix = np.array([[np.sqrt((eta_reco - eta_gen)**2+(phi_reco - phi_gen)**2) for eta_gen, phi_gen in zip(eta_gens, phi_gens)] for eta_reco, phi_reco in zip(eta_recos, phi_recos)])
                ptMatrix = np.array([[abs(ptReco-ptGen)/(ptReco+ptGen) for ptGen in pt_gens] for ptReco in pt_recos])
                print("2 submatrices defined for the vertices ", recoIdx, genIdx)
                if np.any((deltaRMatrix<0.4) & (ptMatrix<0.1)):
                    # there is atleast one track matching that satisfied both condition
                    #print("deltaR matrix 1 passed")
                    #print("delta R matrix min", np.min(deltaRMatrix))
                    while np.any(deltaRMatrix<0.4):
                        recoDeltaRIdx, genDeltaRIdx = np.unravel_index(np.argmin(deltaRMatrix, axis=None), deltaRMatrix.shape)
                        print("deltaRIdx ", recoDeltaRIdx, genDeltaRIdx)
                        #print("reco, gen, value", recoDeltaRIdx, genDeltaRIdx, deltaRMatrix[recoDeltaRIdx, genDeltaRIdx])
                        if ptMatrix[recoDeltaRIdx, genDeltaRIdx]<0.1:
                            #print("Pt matrix 1 passed", ptMatrix[recoDeltaRIdx, genDeltaRIdx])
                            # look for the second track:
                            deltaRMatrix[recoDeltaRIdx, :] = [999]*deltaRMatrix.shape[1]
                            deltaRMatrix[:, genDeltaRIdx] = [999]*deltaRMatrix.shape[0]
                            ptMatrix[recoDeltaRIdx, :] = [999]*ptMatrix.shape[1]
                            ptMatrix[:, genDeltaRIdx] = [999]*ptMatrix.shape[0]
                            while np.any(deltaRMatrix<0.4):
                                recoDeltaRIdx_2, genDeltaRIdx_2 = np.unravel_index(np.argmin(deltaRMatrix, axis=None), deltaRMatrix.shape)
                                if deltaRMatrix[recoDeltaRIdx_2, genDeltaRIdx_2]<0.4:
                                    #print("deltaR matrix 2 passed", deltaRMatrix[recoDeltaRIdx_2, genDeltaRIdx_2])
                                    if ptMatrix[recoDeltaRIdx_2, genDeltaRIdx_2]<0.1:
                                        #print("pass the pt. vertices matched")
                                        #input('Next')
                                        #print("Ev ", ev)
                                        #print("reco pt ", pt_recos[recoDeltaRIdx_2])
                                        #print("pt gen ", pt_gens[genDeltaRIdx_2])
                                        #print("Pt matrix 2 passed", ptMatrix[recoDeltaRIdx_2, genDeltaRIdx_2])
                                        #
                                        distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
                                        distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                                        matchingKey[genIdx]=recoIdx
                                        # to make all the while stop
                                        deltaRMatrix=deltaRMatrix+100
                                        break
                                    else:
                                        #print("pt 2 track non matched. try another second track")
                                        deltaRMatrix[recoDeltaRIdx_2, genDeltaRIdx_2] = 998
                                        if np.min(deltaRMatrix)>0.4:
                                            #every attempt of second track failed. two vertices will not be matched
                                            distances_filled[recoIdx, genIdx]=998
                                else:
                                    deltaRMatrix[recoDeltaRIdx_2, genDeltaRIdx_2] = 998
                                    print("deltaR2 non passed")
                                    #distances_filled[recoIdx, genIdx]=998
                        else:
                            # if the pt is not satisfied reject the daughters
                            print("pt not satisfied, still a minimum?", np.min(deltaRMatrix))
                            #print(deltaRMatrix, "\n\n")
                            print("going to change this value", deltaRMatrix[recoDeltaRIdx, genDeltaRIdx])
                            deltaRMatrix[recoDeltaRIdx, genDeltaRIdx]=998
                            #print(deltaRMatrix[recoDeltaRIdx, genDeltaRIdx])  
                            if np.all((deltaRMatrix>0.4) | (ptMatrix>0.1)):
                                #the best deltaR does not have pt satisfied but a worse deltaR still <0.4 can satisfy this
                                print("Filling the distances_filled matrix")
                                distances_filled[recoIdx, genIdx]=998
                            print(deltaRMatrix)

                            #continue
                else:
                    # if the minimum deltaR is > 0.4 reject the matching of the SV and gV
                    print("deltaR matrix 1 NON passed")
                    distances_filled[recoIdx, genIdx]=998
                    continue
            else:
                distances_filled[recoIdx, genIdx]=998
                continue
            
            

            #distance[0]=distances[recoIdx, genIdx]




        for genIdx, mes in enumerate(mesons):        
            notNeutrinoMask = (abs(GenPart_pdgId_)!=12) & (abs(GenPart_pdgId_)!=14) & (abs(GenPart_pdgId_)!=16)    
            genTracks_pt_ = GenPart_pt_[(GenPart_genPartIdxMother_ == mes) &  (notNeutrinoMask) & (abs(GenPart_eta_)<2.5)]
            genTracks_pdgId_ = GenPart_pdgId_[(GenPart_genPartIdxMother_ == mes) &  (notNeutrinoMask) & (abs(GenPart_eta_)<2.5)]
            if genIdx in matchingKey:
                nRecoTracks[0]=np.sum((abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx]))
                recoTracks_pt_ = ProbeTracks_pt_[(abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx])]
                recoTracks_eta_ = ProbeTracks_eta_[(abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx])]
                recoTracks_phi_ = ProbeTracks_phi_[(abs(ProbeTracks_eta_)<2.5) & (ProbeTracks_matchedToSV_==matchingKey[genIdx])]
            for i in range(nRecoTracks[0]):
                recoTracks_pt.append(recoTracks_pt_[i])
                recoTracks_eta.append(recoTracks_eta_[i])
                recoTracks_phi.append(recoTracks_phi_[i])


            try:
                distance[0]=distances[matchingKey[genIdx], genIdx]
            except:
                distance[0]=-1
            for gp in range(nGenPart_):
                if (GenPart_genPartIdxMother_[gp] == mes):
                    # found the decay of the meson
                    displacement[0]=distance_3d((GenPart_vx_[gp], GenPart_vy_[gp], GenPart_vz_[gp]), (GenPart_vx_[mes], GenPart_vy_[mes], GenPart_vz_[mes]))
                    break
                        #genTracks_pdgId=genTracks_pdgId_
            # meson variables
            event[0]=ev
            nSV[0]      = nSV_
            pdgID[0]    = GenPart_pdgId_[mes]
            pt[0]       = GenPart_pt_[mes]
            eta[0]      = GenPart_eta_[mes]
            phi[0]      = GenPart_phi_[mes]
            vx[0]       = GenPart_vx_[gp]
            vy[0]       = GenPart_vy_[gp]
            vz[0]       = GenPart_vz_[gp]
            nGenTracks[0]=len(genTracks_pt_)
            for dau in range(nGenTracks[0]):
                genTracks_pt.append(genTracks_pt_[dau])
                genTracks_pdgId.append(genTracks_pdgId_[dau])

            # Fill the nested branch with data
        #    for value in nested_data:
        #        nested_branch.push_back(value)
        #    print(x)
            tree.SetBranchAddress('genTracks_pt', genTracks_pt)
            tree.SetBranchAddress('genTracks_pdgId', genTracks_pdgId)
            tree.SetBranchAddress('recoTracks_pt',  recoTracks_pt)
            tree.SetBranchAddress('recoTracks_eta', recoTracks_eta)
            tree.SetBranchAddress('recoTracks_phi', recoTracks_phi)
            tree.Fill()
    tree.Write()
    file.Close()

if __name__ == "__main__":
    fileName = sys.argv[1] if len(sys.argv)>1 else '/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/crab_TTToHadronic/240506_153739/0000/TTToHadronic_Run2_mc_2024May06_1.root'
    fileNumber = int(sys.argv[2]) if len(sys.argv)>1 else 1
    main(fileName, fileNumber)