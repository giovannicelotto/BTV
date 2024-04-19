import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math, sys
import awkward as ak
def distance_3d(point1, point2):
    """
    Calculates the distance between two points in three-dimensional space.
    
    Args:
        point1 (tuple): The coordinates of the first point (x, y, z).
        point2 (tuple): The coordinates of the second point (x, y, z).
        
    Returns:
        float: The distance between the two points.
    """
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

def main():
    f = uproot.open("/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/TTTo2L2Nu_Run2_mc_124X.root")
    tree = f['Events']
    branches = tree.arrays()
    GenPart_pdgId               = branches["GenPart_pdgId"]
    Jet_pt                   = branches["Jet_pt"]
    mask1 = ak.any((abs(GenPart_pdgId) == 521) | (abs(GenPart_pdgId) == 511), axis=-1)
    mask2 = ak.any(Jet_pt>100, axis=-1)
    mask = (mask1) & (mask2) 
    
    for ev in np.arange(tree.num_entries)[mask]:
        print("Event ", ev)
        #sys.exit("exit")
        nJet                        = branches["nJet"][ev]
        Jet_pt                   = branches["Jet_pt"][ev]
        nGenJet                     = branches["nGenJet"][ev]
        GenJet_pt                   = branches["GenJet_pt"][ev]
        GenJet_eta                  = branches["GenJet_eta"][ev]
        GenJet_phi                  = branches["GenJet_phi"][ev]
        GenJet_mass                 = branches["GenJet_mass"][ev]

        GenJet_partonIdx            = branches["GenJet_partonIdx"][ev]
        GenJet_partonFlavour        = branches["GenJet_partonFlavour"][ev]
        GenJet_partonMotherPdgId    = branches["GenJet_partonMotherPdgId"][ev]
        GenJet_partonMotherIdx      = branches["GenJet_partonMotherIdx"][ev]
        Jet_genJetIdx               = branches["Jet_genJetIdx"][ev]
        SV_chi2                     = branches["SV_chi2"][ev]
        SV_pt                       = branches["SV_pt"][ev]
        SV_eta                      = branches["SV_eta"][ev]
        SV_phi                      = branches["SV_phi"][ev]
        SV_mass                     = branches["SV_mass"][ev]
        nSV                         = branches["nSV"][ev]
        SV_dlen                     = branches["SV_dlen"][ev]
        SV_dlenSig                  = branches["SV_dlenSig"][ev]
        SV_dxy                      = branches["SV_dxy"][ev]
        SV_dxySig                   = branches["SV_dxySig"][ev]
        SV_pAngle                   = branches["SV_pAngle"][ev]
        SV_charge                   = branches["SV_charge"][ev]
        SV_ntracks                  = branches["SV_ntracks"][ev]
        SV_ndof                     = branches["SV_ndof"][ev]        

        nGenPart                    = branches["nGenPart"][ev]
        GenPart_pt                  = branches["GenPart_pt"][ev]
        GenPart_eta                 = branches["GenPart_eta"][ev]
        GenPart_phi                 = branches["GenPart_phi"][ev]
        GenPart_mass                = branches["GenPart_mass"][ev]
        GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
        GenPart_pdgId               = branches["GenPart_pdgId"][ev]

        PV_x                        = branches["PV_x"][ev]
        PV_y                        = branches["PV_y"][ev]
        PV_z                        = branches["PV_z"][ev]
        SV_x                        = branches["SV_x"][ev]
        SV_y                        = branches["SV_y"][ev]
        SV_z                        = branches["SV_z"][ev]

        GenPart_vx                  = branches["GenPart_vx"][ev]
        GenPart_vy                  = branches["GenPart_vy"][ev]
        GenPart_vz                  = branches["GenPart_vz"][ev]
        
        for jetIdx in np.arange(nJet)[Jet_pt>100]:
            matched=False
            print("Jet_pt", Jet_pt[jetIdx])
            if GenJet_partonFlavour[jetIdx]!=0:
                assert GenPart_pdgId[GenJet_partonIdx[jetIdx]]==GenJet_partonFlavour[jetIdx]
                if abs(GenPart_pdgId[GenJet_partonIdx[jetIdx]])==5:
                    if abs(GenPart_pdgId[GenPart_genPartIdxMother[GenJet_partonIdx[jetIdx]]])==6:
                        matched=True
                    elif abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenJet_partonIdx[jetIdx]]]])==6:
                        matched=True
                    else:
                        matched=False
                if matched:
                    print("PDG mothers", GenPart_pdgId[GenJet_partonIdx[jetIdx]], GenPart_pdgId[GenPart_genPartIdxMother[GenJet_partonIdx[jetIdx]]], GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenJet_partonIdx[jetIdx]]]])

                    daughters=[]    # all the daughters of 5 degrees
                    motherRef = [GenJet_partonIdx[jetIdx]]  # all the mothers we want to find daughters of

                    for LL in range(10):
                        print(LL, motherRef)
                        daughterLevel=[]
                        for gp in range(nGenPart):
                            if (GenPart_genPartIdxMother[gp] in motherRef):
                                daughterLevel.append(gp)
                        
                        
                        motherRef = daughterLevel
                        
                        daughters.append(daughterLevel)
                    #print(daughters)
                    print(SV_x)
                    for d in daughters:
                        print(GenPart_pdgId[d], GenPart_vx[d])
                    
                            
                            

                else:
                    print("No top in sequence")
            else:
                print("No matching")
            


            



        #print(nSV, nJet)
#        svs = []
        
        #print("(%.4f, %.4f, %.4f)  %d SVs"%(PV_x, PV_y, PV_z, nSV))
        for svIdx in range(nSV):
            print("SV ", svIdx, (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx]) )
            for genPartIdx in np.arange(nGenPart):
                if abs(GenPart_pdgId[GenPart_genPartIdxMother[genPartIdx]])==521:
                    print("GenPart %d "%(GenPart_pdgId[genPartIdx]), (GenPart_vx[genPartIdx], GenPart_vy[genPartIdx], GenPart_vz[genPartIdx]))
           #svs.append( (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx]))
        BMesonDaughters = abs(GenPart_pdgId[GenPart_genPartIdxMother])==521
        print(BMesonDaughters)
        
        for bmesDau in np.arange(nGenPart)[BMesonDaughters]:
            print(GenPart_vx[bmesDau])
            dist=999
            for jetIdx in np.arange(nJet):
                if Jet_genJetIdx[jetIdx]!=-1:
                    partonIdx = GenJet_partonIdx[Jet_genJetIdx[jetIdx]]
                    if partonIdx!=-1:
                        partonPdgId = GenPart_pdgId[partonIdx]
                        vertex = (GenPart_vx[partonIdx], GenPart_vy[partonIdx], GenPart_vz[partonIdx])
                        print(GenPart_pdgId[GenPart_genPartIdxMother[partonIdx]] ,partonIdx, partonPdgId, "True : ", vertex)
                        print(GenPart_vx[GenPart_genPartIdxMother[partonIdx]])
                        print(distance_3d(vertex, (SV_x[svIdx], SV_y[svIdx], SV_z[svIdx])))
            

            #input("\n\nNext\n")

                #mathced
                


    return 



if __name__ =="__main__":
    main()