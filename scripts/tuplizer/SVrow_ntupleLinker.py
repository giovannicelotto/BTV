import ROOT
from array import array
import uproot
import numpy as np
from utilsForScript import distance_3d, getPdgMask
import sys
from ntupleLinker import matchingEvent, getOneDaughter, getParams, getMesons, getSecondMatching

def main(fileName, fileNumber, prova):

    if prova==1:
        fN = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/cutSVOnly/flatSV/TTToH_SVrow_%d.root"%fileNumber
        file = ROOT.TFile(fN, "RECREATE")    
        print(fN)
    elif prova==0:
        file = ROOT.TFile("/scratch/TTToH_SVrow_%d.root"%fileNumber, "RECREATE")
        print("Created root file in /scratch/TTToH_SVrow_%d.root"%fileNumber)
    else:
        print("Invalid value for prova")
    tree = ROOT.TTree("tree", "Tree")
    params = getParams()
    print("Params done")

    event               = array('l', [0])
    nSV                 = array('l', [0])
    SV_idx               = array('l', [0])    
    SV_matched         = array('l', [0])    
    #SV_pt               = array('d', [0])
    #SV_eta              = array('d', [0])    
    #SV_phi              = array('d', [0])    
    SV_x                = array('d', [0])
    SV_y                = array('d', [0])
    SV_z                = array('d', [0])
    SV_dlen             = array('d', [0])
    SV_dxy              = array('d', [0])
    SV_dxySig           = array('d', [0])
    SV_dlenSig          = array('d', [0])   
                

    tree.Branch("event",              event,            "event/I")                     
    tree.Branch("nSV",                nSV,              "nSV/I")  
    tree.Branch("SV_idx",             SV_idx,           "SV_idx/I")      
    tree.Branch("SV_matched",         SV_matched,       "SV_matched/I")      
    #tree.Branch("SV_pt",              SV_pt,                 "SV_pt/D")
    #tree.Branch("SV_eta",             SV_eta,                "SV_eta/D")    
    #tree.Branch("SV_phi",             SV_phi,                "SV_phi/D")    
    tree.Branch("SV_x",               SV_x,             "SV_x/D")
    tree.Branch("SV_y",               SV_y,             "SV_y/D")
    tree.Branch("SV_z",               SV_z,             "SV_z/D")
    tree.Branch("SV_dlen",         SV_dlen,           "SV_dlen/D")
    tree.Branch("SV_dxy",           SV_dxy,            "SV_dxy/D")
    tree.Branch("SV_dxySig",     SV_dxySig,         "SV_dxySig/D")
    tree.Branch("SV_dlenSig",         SV_dlenSig,           "SV_dlenSig/D")        


    branches_to_open = ["nSV", "nGenPart", "GenPart_genPartIdxMother", "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_status", "GenPart_charge", "SV_x", "SV_y", "SV_z", "GenPart_vx", "GenPart_vy", "GenPart_vz", "ProbeTracks_matchedToSV", "ProbeTracks_pt", "ProbeTracks_eta", "ProbeTracks_phi", "ProbeTracks_genPartIdx", "SV_dlen", "SV_dlenSig", "SV_dxy", "SV_dxySig"]
    f = uproot.open(fileName)
    tree_ = f['Events']
    branches = tree_.arrays(branches_to_open)
    maxEntries = tree_.num_entries if prova == 1 else tree_.num_entries
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
        GenPart_status_             = branches["GenPart_status"][ev]
        GenPart_charge              = branches["GenPart_charge"][ev]
        SV_x_                       = branches["SV_x"][ev]
        SV_y_                       = branches["SV_y"][ev]
        SV_z_                       = branches["SV_z"][ev]
        GenPart_vx_                 = branches["GenPart_vx"][ev]
        GenPart_vy_                 = branches["GenPart_vy"][ev]
        GenPart_vz_                 = branches["GenPart_vz"][ev]
        ProbeTracks_matchedToSV_    = branches["ProbeTracks_matchedToSV"][ev]
        ProbeTracks_pt_             = branches["ProbeTracks_pt"][ev]
        ProbeTracks_eta_            = branches["ProbeTracks_eta"][ev]
        ProbeTracks_phi_            = branches["ProbeTracks_phi"][ev]
        ProbeTracks_genPartIdx      = branches["ProbeTracks_genPartIdx"][ev]
        SV_dlen_                    = branches["SV_dlen"][ev]
        SV_dlenSig_                 = branches["SV_dlenSig"][ev]
        SV_dxy_                     = branches["SV_dxy"][ev]
        SV_dxySig_                  = branches["SV_dxySig"][ev]

        # SVs as reconstructed by IVF
        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x_, SV_y_, SV_z_)])

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
        #print("Distances ", distances)
        matchingKey = matchingEvent(distances, ProbeTracks_matchedToSV_, ProbeTracks_pt_, ProbeTracks_eta_, GenPart_genPartIdxMother_, oneDaughter, ProbeTracks_genPartIdx, nGenPart_)
        # matchingKey[genIdx]=recoIdx
        
         
        #print("Mathcing key")
        #print(matchingKey)
        for svIdx in range(nSV_):        
            event[0] = ev
            nSV[0] = nSV_
            SV_idx[0] = svIdx
            if svIdx not in matchingKey.values():
                # SV was not matched
                SV_matched[0] = 0
            else:
                SV_matched[0] = 1
            SV_x[0]         = SV_x_[svIdx]
            SV_y[0]         = SV_y_[svIdx]
            SV_z[0]         = SV_z_[svIdx]
            SV_dlen[0]      = SV_dlen_[svIdx]
            SV_dlenSig[0]   = SV_dlenSig_[svIdx]
            SV_dxy[0]       = SV_dxy_[svIdx]
            SV_dxySig[0]    = SV_dxySig_[svIdx]
            
            tree.Fill()
    tree.Write()
    file.Close()


if __name__ == "__main__":
    fileName = sys.argv[1] if len(sys.argv)>1 else '/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/cutSVOnly/nano/TTToHadronic_Run2_cutSVOnly.root'
    fileNumber = int(sys.argv[2]) if len(sys.argv)>1 else 100
    prova = 0 if len(sys.argv)>1 else 1
    main(fileName, fileNumber, prova)