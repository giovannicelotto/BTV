import ROOT
from array import array
import uproot
import numpy as np
from helpers.utilsForScript import distance_3d, getPdgMask
from helpers.getGV import getMesons
from helpers.getOneDaughter import getOneDaughter
from helpers.matchingEvent import matchingEvent
from helpers.getParams import getParams
import argparse
import sys

def getSecondMatching():
    secondMatching = False
    return secondMatching
secondMatching = getSecondMatching()


def main(fileName, fileNumber, prova, maxEntries):
    if prova==1:
        fN = "/work/gcelotto/BTV/scripts/tuplizer/TTToH_selPart_%d.root"%fileNumber
        file = ROOT.TFile(fN, "RECREATE")    
        print(fN)
    elif prova==0:
        file = ROOT.TFile("/scratch/TTToH_%d.root"%(fileNumber), "RECREATE")
        print("Created root file in /scratch/TTToH_%d.root"%(fileNumber))
    else:
        print("Invalid value for prova")
    tree = ROOT.TTree("tree", "Tree")
    params = getParams()

    event = array('l', [0])
    nSV = array('l', [0])
    pdgID               = array('i', [0])    
    motherPdgID         = array('i', [0])    
    pt                  = array('d', [0])
    eta                 = array('d', [0])    
    phi                 = array('d', [0])    
    vx                  = array('d', [0])
    vy                  = array('d', [0])
    vz                  = array('d', [0])
    SV_dlenSig          = array('d', [0])
    displacement        = array('d', [0])            
        # status of D mesons and Charmed Baryons: 
    # 0 prompt (from W)
    # 1,(2),[3],... if B [grand(grand)]mother mother of D
    # -2 neither from W nor from B (from hard scattering e.g.)
    # -1 for all the particles different from D mesons and Charmed Baryons
    promptStatus        = array('l', [0])   
        # status of B mesons and Bottom baryons: 
    # -1 not B nor Bottom baryons
    # 0 does not have D below in the chain
    # 1 mother of D or Charmed baryons
    # 2 grandmother of D or Charmed baryons
    # ...
    decayStatus        = array('l', [0])   


    B_eta               = array('d', [0])   # eta of the b meson if the D comes from a B.
    B_phi               = array('d', [0])   # phi of the b meson if the D comes from a B.
    nGenTracks          = array('i', [0] )
    genTracks_pt        = array('d', [-999.]*15)
    genTracks_eta       = array('d', [-999.]*15)
    genTracks_phi       = array('d', [-999.]*15)
    genTracks_pdgId     = array('d', [-999.]*15)     
    genTracks_charge    = array('i', [-999]*15)     
    nRecoTracks         = array('i', [ 0 ] )
    recoTracks_pt       = array('d', [0.]*30)
    recoTracks_eta      = array('d', [0.]*30)                
    recoTracks_phi      = array('d', [0.]*30)                
    recoTracks_pdgID    = array('i', [0]*30)
    distance            = array('d', [0])        
    delta_x            = array('d', [0])        
    delta_y            = array('d', [0])        
    delta_z            = array('d', [0])      
    
             

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
    tree.Branch("recoTracks_pdgID",  recoTracks_pdgID,    "recoTracks_pdgID[nRecoTracks]/I")                
    tree.Branch("distance",         distance,           "distance/D")        
    tree.Branch("delta_x",         delta_x,           "delta_x/D")        
    tree.Branch("delta_y",         delta_y,           "delta_y/D")        
    tree.Branch("delta_z",         delta_z,           "delta_z/D")        
    tree.Branch("SV_dlenSig",         SV_dlenSig,           "SV_dlenSig/D")        



    f = uproot.open(fileName)
    tree_ = f['Events']
    branches = tree_.arrays()
    if maxEntries == -1:
        maxEntries = tree_.num_entries
    print("%d entries "%maxEntries)
    for ev in range(maxEntries):
        if ev % 100 == 0:
            sys.stdout.write(f"\rProgress: {100 * ev / 1000:.2f}%")
            sys.stdout.flush()
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

        nsvDaughters                = branches["nsvDaughters"][ev]
        svDaughters_svIdx           = branches["svDaughters_svIdx"][ev]
        svDaughters_eta             = branches["svDaughters_eta"][ev]
        svDaughters_pt              = branches["svDaughters_pt"][ev]
        svDaughters_phi             = branches["svDaughters_phi"][ev]
        svDaughters_genPartIdx      = branches["svDaughters_genPartIdx"][ev]

        SV_dlenSig_                 = branches["SV_dlenSig"][ev]

        # SVs as reconstructed by IVF
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])

# ***********************************************************************
#                      Gen Vertices definition                         
# ***********************************************************************
        
        mesons, mesonsDaughters = getMesons(params, GenPart_pdgId_, GenPart_eta_, GenPart_pt_, nGenPart_, GenPart_status_, GenPart_charge, GenPart_genPartIdxMother_)     

        oneDaughter = getOneDaughter(mesons, GenPart_genPartIdxMother_, nGenPart_)
        genVertices = np.array([(x, y, z) for x, y, z in zip(GenPart_vx_[oneDaughter], GenPart_vy_[oneDaughter], GenPart_vz_[oneDaughter])])

# ***********************************************************************
#                            Matching                         
# ***********************************************************************
        distances = np.array([[distance_3d(sv, (vx, vy, vz)) for vx, vy, vz in genVertices] for sv in SVs ])
        #if SV_dlenSig_cut is not None:
        #    for recoIdx in range(nSV_):
        #        if SV_dlenSig_[recoIdx]<SV_dlenSig_cut:
        #            distances[recoIdx, :]=[998]*distances.shape[1]
        matchingKey = matchingEvent(distances, svDaughters_svIdx, svDaughters_pt, svDaughters_eta, GenPart_genPartIdxMother_, oneDaughter, svDaughters_genPartIdx, nGenPart_, GenPart_pdgId_, mesonsDaughters)
        #input("Next")
        #print(matchingKey)
         

        for genIdx, mes in enumerate(mesons):        
            #notNeutrinoMask = (abs(GenPart_pdgId_)!=12) & (abs(GenPart_pdgId_)!=14) & (abs(GenPart_pdgId_)!=16)    
            genDaughtersMask = np.isin(np.arange(nGenPart_), mesonsDaughters[genIdx])
            
            
            # No here need to change. Retrieve the information from the matching about the daughters!
            # same daughters used for genVertex definition




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
                recoTracksMask = (abs(svDaughters_eta)<2.5) & (svDaughters_svIdx==matchingKey[genIdx]) & (svDaughters_pt>0.8)
                nRecoTracks[0]=np.sum(recoTracksMask)
                recoTracks_pt_ = svDaughters_pt[recoTracksMask]
                recoTracks_eta_ = svDaughters_eta[recoTracksMask]
                recoTracks_phi_ = svDaughters_phi[recoTracksMask]
                recoTracks_pdgID_ = GenPart_pdgId_[svDaughters_genPartIdx[recoTracksMask]]
                recoTracks_pdgID_ = np.where(svDaughters_genPartIdx[recoTracksMask]>-1, recoTracks_pdgID_, -1)
            else:
                nRecoTracks[0]=0
            for i in range(nRecoTracks[0]):
                recoTracks_pt[i] = recoTracks_pt_[i]
                recoTracks_eta[i] = recoTracks_eta_[i]
                recoTracks_phi[i] = recoTracks_phi_[i]
                recoTracks_pdgID[i] = recoTracks_pdgID_[i]


            if genIdx in matchingKey:
                distance[0]=distances[matchingKey[genIdx], genIdx]
                delta_x[0] = SV_x[matchingKey[genIdx]] - GenPart_vx_[oneDaughter[genIdx]]
                delta_y[0] = SV_y[matchingKey[genIdx]] - GenPart_vy_[oneDaughter[genIdx]]
                delta_z[0] = SV_z[matchingKey[genIdx]] - GenPart_vz_[oneDaughter[genIdx]]
                SV_dlenSig[0] = SV_dlenSig_[matchingKey[genIdx]]
            else:
                distance[0]=-1
                delta_x[0] = -1
                delta_y[0] = -1
                delta_z[0] = -1
                SV_dlenSig[0] = -1
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script.")
    #### Define arguments
    parser.add_argument("-fileName", "--fileName", type=str, default="/work/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_Run2_mc_124X.root")
    parser.add_argument("-fileNumber", "--fileNumber", type=int, default=100)
    parser.add_argument("-p", "--prova", type=int, help="prova True (1) False (0)", default=1)
    parser.add_argument("-mE", "--maxEntries", type=int, help="maxEntries per File", default=500)
    args = parser.parse_args()
    main(args.fileName, args.fileNumber, args.prova, args.maxEntries)