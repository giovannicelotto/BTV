import numpy as np
import uproot
import glob
import sys
from efficiency_differential import map_to_groups
from efficiency_differential import efficiencyVsVariable, plotPtEtaMatchedVsNonMatched
import os
map_to_groups_vec = np.vectorize(map_to_groups)
def main(nFiles):
    rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06Tuple"
    fileNames = glob.glob(rootDir+"/*.root")
    if len(fileNames)<nFiles:
        nFiles = len(fileNames)
    else:
        fileNames=fileNames[:nFiles]
    print("Taking %d files"%len(fileNames))
    for fileIdx, fileName in enumerate(fileNames):
        print("Opening %s"%fileName)
        
        f = uproot.open(fileName)
        try:
            tree=f['tree']
        except:
            print("Skipping file. Broken")
            os.remove(fileName)
            continue
        print("File %d"%fileIdx)
        if fileIdx == 0:
            
            distance = np.array(tree['distance'])
            pdgID   =   np.array(tree['pdgID']).astype(int)
            eta   =   np.array(tree['eta'])
            pt   =   np.array(tree['pt'])
            displacement   =   np.array(tree['displacement'])
            nGenTracks   =   np.array(tree['nGenTracks'])
        else:
            
            distance = np.concatenate((distance, np.array(tree['distance'])))
            pdgID = np.concatenate((pdgID, np.array(tree['pdgID']).astype(int)))
            eta = np.concatenate((eta, np.array(tree['eta'])))
            pt = np.concatenate((pt, np.array(tree['pt'])))
            displacement   =   np.concatenate((displacement, tree['displacement']))
            nGenTracks   =   np.concatenate((nGenTracks, tree['nGenTracks']))
            
        #matched = np.array(tree['matched'])


    #print(distance[distance>0])
    mesons = map_to_groups_vec(pdgID)
    commonMask = displacement>0.


    efficiencyVsVariable(num=eta[(distance>0) & (commonMask)], den=eta[commonMask],
                        bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/etaDiff.png", title = "All particles",)

    efficiencyVsVariable(num=pt[(distance>0) & (commonMask)], den=pt[commonMask],
                        bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/ptDiff.png", title = "All particles",)

    efficiencyVsVariable(num=displacement[(distance>0) & (commonMask)], den=displacement[commonMask],
                        bins=np.linspace(0, 1, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/displacementDiff.png", title = "All particles",)

    efficiencyVsVariable(num=nGenTracks[(distance>0) & (commonMask)], den=nGenTracks[commonMask],
                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/nGenTracksDiff.png", title = "All particles",)

    efficiencyVsVariable(num=mesons[(distance>0) & (commonMask)], den=mesons[commonMask],
                        bins=np.linspace(-1.5, 4.5, 7), xlabel="Hadron PDG Class",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/pdgDiff.png", title = "All particles",
                        tick_positions = [0, 1, 2, 3, 4],
                        tick_labels = ['B mesons', 'D mesons', 'Strange Baryons', 'Charmed Baryons', 'Bottom Baryons'],)


    # B only

    efficiencyVsVariable(num=eta[(distance>0) & (mesons==0) & (commonMask)], den=eta[(mesons==0) & (commonMask)],
                            bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/etaDiff.png", title="B only",)
    efficiencyVsVariable(num=pt[(distance>0) & (mesons==0) & (commonMask)], den=pt[(mesons==0) & (commonMask)],
                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/ptDiff.png", title="B only",)
    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==0) & (commonMask)], den=displacement[(commonMask) & (mesons==0)],
                        bins=np.linspace(0, 1, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/displacementDiff.png", title = "B only",)
    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==0) & (commonMask)], den=nGenTracks[(commonMask) & (mesons==0)],
                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/nGenTracksDiff.png", title = "B only",)


    # D only
    efficiencyVsVariable(num=eta[(distance>0) & (mesons==1) & (commonMask)], den=eta[(mesons==1) & (commonMask)],
                            bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/etaDiff.png", title="D only")
    efficiencyVsVariable(num=pt[(distance>0) & (mesons==1) & (commonMask)], den=pt[(mesons==1) & (commonMask)],
                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/ptDiff.png", title="D only")
    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==1) & (commonMask)], den=displacement[(commonMask) & (mesons==1)],
                        bins=np.linspace(0, 0.5, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/displacementDiff.png", title = "D only",)
    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==1) & (commonMask)], den=nGenTracks[(commonMask) & (mesons==1)],
                        bins=np.arange(2, 10), xlabel="Number of hadron daughters",
                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/nGenTracksDiff.png", title = "D only",)

    print("***********\nYou should add an option for the overflow True False in the efficiencyVsVariable function\n****************")


    print("Overall efficiency: %.2f"%(np.sum(distance>0)/len(distance)*100))





    return 0


if __name__=="__main__":
    nFiles = int(sys.argv[1]) if len(sys.argv)>1 else 1
    main(nFiles=nFiles)

