import numpy as np
import uproot
import glob
import sys
#from efficiency_differential import map_to_groups
#from efficiency_differential import efficiencyVsVariable, plotPtEtaMatchedVsNonMatched
import os
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
eps = 0.000001
class_to_label = {
    0: 'B',
    1: 'D',
    2: 'SB',
    3: 'CB',
    4: 'BB'
}


def map_to_groups(value):
    if abs(value) in [511, 521, 531, 541]:
        return 0
    elif abs(value) in [411, 421, 431]:
        return 1
    elif ((abs(value) > 3000) & (abs(value) < 4000)):
        return 2
    elif ((abs(value) > 4000) & (abs(value) < 5000)):
        return 3
    elif ((abs(value) > 5000) & (abs(value) < 6000)):
        return 4
    # Add more conditions as needed
    else:
        return -1  # or any default value for unmatched cases
#
# 0 -> B
# 1 -> D
# 2 -> SB
# 3 -> CB
# 4 -> BB
map_to_groups_vec = np.vectorize(map_to_groups)
def efficiencyVsVariable(num, den, bins, xlabel, outName, title=None, tick_positions=None, tick_labels=None):

    fig,ax = plt.subplots(1, 1)
    matched = np.histogram(num, bins=bins)[0]
    total =   np.histogram(den, bins=bins)[0]
    ax.errorbar((bins[:-1]+bins[1:])/2, matched/(total+eps), xerr=np.diff(bins)/2,yerr=np.sqrt(matched)/(total+eps), marker='o', linestyle='none', color='black')
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Efficiency [%]")
    ax.set_ylim(0, 1)
    ax.text(x=0.95, y=0.94, s="Matched Entries %d"%(len(num)), transform=ax.transAxes, horizontalalignment='right')
    ax.text(x=0.95, y=0.89, s=title, transform=ax.transAxes, horizontalalignment='right')
    #print(tick_positions)
    if tick_positions is not None:
        #print("tick_positions customized")
        ax.set_xticks(tick_positions, tick_labels,  rotation=30)
    #if title is not None:
    hep.cms.label()
    fig.savefig(outName, bbox_inches='tight')
    print("Saved %s"%outName)
    return
def flightMatchedDistance(distance, displacement, bins, outFile):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    matchedMask = distance>-1
    cM = np.histogram(displacement[matchedMask], bins=bins)[0]
    cM=cM/np.sum(cM)

    cNM = np.histogram(displacement[~matchedMask], bins=bins)[0]
    cNM=cNM/np.sum(cNM)

    ax.hist(bins[:-1], bins=bins, weights=cM, alpha=0.4, label='Matched', linewidth=3)
    ax.hist(bins[:-1], bins=bins, weights=cNM, alpha=0.4, label='Non Matched', linewidth=3)
    ax.hist(bins[:-1], bins=bins, weights=cM, histtype=u'step', linewidth=3, color='C0')
    ax.hist(bins[:-1], bins=bins, weights=cNM, histtype=u'step', linewidth=3, color='C1')
    ax.legend()
    hep.cms.label()
    ax.set_xlabel("Distance of flight [cm]")
    fig.savefig(outFile, bbox_inches='tight')
    print("Saving ", outFile)
    return

def distanceVsPDGClass(distance, pdgClass, bins, outFile):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    uniqueClasses = set(pdgClass)
    fig, ax = plt.subplots(1, 1)

    for myclass in uniqueClasses:
        counts = np.histogram(distance[pdgClass==myclass], bins=bins)[0]
        counts=counts/np.sum(counts)
        label = class_to_label[myclass]
        ax.hist(bins[:-1], bins=bins, weights=counts, histtype=u'step', label=label)

    hep.cms.label()
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlabel("Distance [cm]")
    ax.set_ylabel("Normalized Events")
    ax.legend()
    # Define the region you want to zoom in on
    x1, x2, y1, y2 = 0, 0.2, 0.1, 0.35
    axins = ax.inset_axes([0.3, 0.5, 0.4, 0.4], xlim=(x1, x2), ylim=(y1, y2),  xticklabels=[], yticklabels=[])
    for myclass in uniqueClasses:
        counts = np.histogram(distance[pdgClass==myclass], bins=bins)[0]
        counts=counts/np.sum(counts)
        axins.hist(bins[:-1], bins=bins, weights=counts, histtype=u'step')
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    ax.indicate_inset_zoom(axins, edgecolor='black')



    fig.savefig(outFile)
    print("Saved ", outFile)

def cumulativeDistance(distance, pdgClass, bins, outFile):
    # Second plot:
    # Cumulative distribution as a function of the distance
    uniqueClasses = set(pdgClass)
    fig, ax = plt.subplots(1, 2, figsize=(15, 8))
    #x = np.linspace(0, 2, 200)
    x = np.logspace(np.log10(1e-4), np.log10(5), num=200)
    for myclass in uniqueClasses:
        y = []
        y2 = []
        mask = (pdgClass==myclass) & (distance>0)
        maskDen = (pdgClass==myclass)
        for i in x:
            y.append(np.sum(distance[mask]<i)/len(distance[mask]))
            y2.append(np.sum(distance[mask]<i)/len(distance[maskDen]))
        label = class_to_label[myclass]
        
        ax[0].plot(x, y, label = label)
        ax[1].plot(x, y2, label = label)
    ax[0].set_xlabel("distance [cm]")
    ax[0].set_ylabel("Matched / total Matched")
    ax[0].legend()
    ax[0].set_ylim(0.5, 1)
    ax[1].set_ylim(0., 0.6)
    ax[1].legend()
    ax[1].set_xlabel("distance [cm]")
    ax[1].set_ylabel("Matched/total")
    fig.savefig(outFile, bbox_inches='tight')
    print("Saved ", outFile)

        


def main(nFiles):
    rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple"
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


#    efficiencyVsVariable(num=eta[(distance>0) & (commonMask)], den=eta[commonMask],
#                        bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/etaDiff.png", title = "All particles",)
#
#    efficiencyVsVariable(num=pt[(distance>0) & (commonMask)], den=pt[commonMask],
#                        bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/ptDiff.png", title = "All particles",)
#
#    efficiencyVsVariable(num=displacement[(distance>0) & (commonMask)], den=displacement[commonMask],
#                        bins=np.linspace(0, 1, 20), xlabel="Hadron distance of flight [cm]",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/displacementDiff.png", title = "All particles",)
#
#    efficiencyVsVariable(num=nGenTracks[(distance>0) & (commonMask)], den=nGenTracks[commonMask],
#                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/nGenTracksDiff.png", title = "All particles",)
#
#    efficiencyVsVariable(num=mesons[(distance>0) & (commonMask)], den=mesons[commonMask],
#                        bins=np.linspace(-1.5, 4.5, 7), xlabel="Hadron PDG Class",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/all/pdgDiff.png", title = "All particles",
#                        tick_positions = [0, 1, 2, 3, 4],
#                        tick_labels = ['B mesons', 'D mesons', 'Strange Baryons', 'Charmed Baryons', 'Bottom Baryons'],)
#
#
#    # B only
#
#    efficiencyVsVariable(num=eta[(distance>0) & (mesons==0) & (commonMask)], den=eta[(mesons==0) & (commonMask)],
#                            bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
#                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/etaDiff.png", title="B only",)
#    efficiencyVsVariable(num=pt[(distance>0) & (mesons==0) & (commonMask)], den=pt[(mesons==0) & (commonMask)],
#                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
#                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/ptDiff.png", title="B only",)
#    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==0) & (commonMask)], den=displacement[(commonMask) & (mesons==0)],
#                        bins=np.linspace(0, 1, 20), xlabel="Hadron distance of flight [cm]",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/displacementDiff.png", title = "B only",)
#    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==0) & (commonMask)], den=nGenTracks[(commonMask) & (mesons==0)],
#                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Bonly/nGenTracksDiff.png", title = "B only",)
#
#
#    # D only
#    efficiencyVsVariable(num=eta[(distance>0) & (mesons==1) & (commonMask)], den=eta[(mesons==1) & (commonMask)],
#                            bins=np.linspace(-2.7, 2.7, 20), xlabel="Hadron eta",
#                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/etaDiff.png", title="D only")
#    efficiencyVsVariable(num=pt[(distance>0) & (mesons==1) & (commonMask)], den=pt[(mesons==1) & (commonMask)],
#                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
#                            outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/ptDiff.png", title="D only")
#    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==1) & (commonMask)], den=displacement[(commonMask) & (mesons==1)],
#                        bins=np.linspace(0, 0.5, 20), xlabel="Hadron distance of flight [cm]",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/displacementDiff.png", title = "D only",)
#    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==1) & (commonMask)], den=nGenTracks[(commonMask) & (mesons==1)],
#                        bins=np.arange(2, 10), xlabel="Number of hadron daughters",
#                        outName="/t3home/gcelotto/BTV/plots/newTuple/efficiency/Donly/nGenTracksDiff.png", title = "D only",)
#    
#
    flightMatchedDistance(distance=distance, displacement=displacement, bins=np.linspace(0, .3, 31),outFile="/t3home/gcelotto/BTV/plots/newTuple/displacementM_vs_NM.png")

    distanceVsPDGClass(distance=distance, pdgClass=mesons, bins=np.linspace(0, 2, 100), outFile="/t3home/gcelotto/BTV/plots/newTuple/distancePDG.png")
    cumulativeDistance(distance=distance, pdgClass=mesons, bins=np.linspace(0, 2, 100), outFile="/t3home/gcelotto/BTV/plots/newTuple/cumDistancePDG.png")
    
    print("***********\nYou should add an option for the overflow True False in the efficiencyVsVariable function\n****************")


    print("Overall efficiency: %.2f"%(np.sum(distance>0)/len(distance)*100))





    return 0


if __name__=="__main__":
    nFiles = int(sys.argv[1]) if len(sys.argv)>1 else 1
    main(nFiles=nFiles)

