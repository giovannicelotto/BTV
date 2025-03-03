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
sys.path.append('/work/gcelotto/BTV/scripts/tuplizer')
from ntupleLinker import getSecondMatching
eps = 0.000001

def getClassToLabel():
    class_to_label = {
        0: 'B',
        1: 'D',
        2: 'SB',
        3: 'CB',
        4: 'BB',
        -200: 'D others',
        -2000: 'CB others',
        -1: 'Others',
        100: 'D prompt',
        101: 'D from B-1',
        102: 'D from B-2',
        103: 'D from B-3',
        104: 'D from B-4',
        105: 'D from B-5',
        1000: 'CB prompt',
        1001: 'CB from B-1',
        1002: 'CB from B-2',
        1003: 'CB from B-3',
        1004: 'CB from B-4',
        1005: 'CB from B-5',

    }
    return class_to_label
class_to_label = getClassToLabel()


def getClassToLabelForLegend():
    class_to_label = {
        1: 'D Mesons',
        0: 'B Mesons',
        #2: 'SB',
        3: 'Charmed Baryons',
        4: 'Bottom Baryons',
    }
    return class_to_label
class_to_labelForLegend = getClassToLabelForLegend()

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
    
def map_to_groups_prompt(value, promptStatus):
    if abs(value) in [511, 521, 531, 541]:
        return 0
    elif abs(value) in [411, 421, 431]:
        if promptStatus==-2: # neither W nor B in the chain
            return promptStatus*100 # -200
        elif promptStatus>=0:
            return 100+promptStatus
        else:
            assert False
    elif ((abs(value) > 3000) & (abs(value) < 4000)):
        return 2
    elif ((abs(value) > 4000) & (abs(value) < 5000)):
        if promptStatus==-2: # neither W nor B in the chain
            return promptStatus*1000 # -2000
        elif promptStatus>=0:
            return 1000+promptStatus
        else:
            assert False
    elif ((abs(value) > 5000) & (abs(value) < 6000)):
        return 4
    # Add more conditions as needed
    else:
        return -1
#
# 0 -> B
# 1 -> D
# 2 -> SB
# 3 -> CB
# 4 -> BB
map_to_groups_vec = np.vectorize(map_to_groups)
def efficiencyVsVariable(num, den, bins, xlabel, outName, title=None, underOverFlow = False, tick_positions=None, tick_labels=None):

    fig,ax = plt.subplots(1, 1)
    matched = np.histogram(num, bins=bins)[0]
    total =   np.histogram(den, bins=bins)[0]

    if underOverFlow:
        matched = np.histogram(np.clip(num, bins[0], bins[-1]), bins=bins)[0]
        total =   np.histogram(np.clip(den, bins[0], bins[-1]), bins=bins)[0]


    bins_center = (bins[:-1]+bins[1:])/2
    ax.errorbar(bins_center, matched/(total+eps), xerr=np.diff(bins)/2,yerr=np.sqrt(matched)/(total+eps), marker='o', linestyle='none', color='black')
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Efficiency [%]")
    ax.set_ylim(0, 1.2)
    ax.text(x=0.95, y=0.94, s="Matched Entries %d"%(len(num)), transform=ax.transAxes, horizontalalignment='right')
    ax.text(x=0.95, y=0.89, s=title, transform=ax.transAxes, horizontalalignment='right')
    #print(tick_positions)
    if tick_positions is not None:
        ax.set_xticks(tick_positions, tick_labels,  rotation=30)
        for x,y in zip(bins_center, matched/(total+eps)):
            if y!=0:
                ax.text(x=x, y=y+0.05, s="%.1f%%"%(y*100), ha='center')

    #if title is not None:
    hep.cms.label()
    fig.savefig(outName, bbox_inches='tight')
    print("Saved %s"%outName)
    plt.close()
    return
def flightMatchedDistance(distance, displacement, bins, outFile):
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
    ax.set_ylabel("Normalized Events")
    fig.savefig(outFile, bbox_inches='tight')
    print("Saved ", outFile)
    return

def distanceVsPDGClass(distance, pdgClass, bins, outFile, zoom=True):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    uniqueClasses = set(pdgClass)
    fig, ax = plt.subplots(1, 1)

    for myclass in uniqueClasses:

        counts = np.histogram(np.clip(distance[pdgClass==myclass], bins[0], bins[-1]), bins=bins)[0]
        counts=counts/np.sum(counts)
        label = class_to_labelForLegend[myclass]
        ax.hist(bins[:-1], bins=bins, weights=counts, histtype=u'step', linewidth=2, label=label)

    hep.cms.label()
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlabel("Distance [cm]")
    ax.set_ylabel("Normalized Events")
    ax.legend()
    # Define the region you want to zoom in on
    if zoom:
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

def plotNTracks_Prompt(nRecoTracks, promptStatus, distance, bins, xlabel, outName):
    possiblePrompts = [0, 1, 2, 3]
    fig, ax = plt.subplots(1, 1)
    for pp in possiblePrompts:
        if distance is not None:
            mask = (promptStatus==pp) & (distance>0)
        else:
            mask = (promptStatus==pp)
        c = np.histogram(np.clip(nRecoTracks[mask], bins[0], bins[-1]), bins)[0]
        c = c / np.sum(c)
        label = 'Prompt' if pp==0 else "From B-%d"%pp
        ax.hist(bins[:-1], bins=bins, weights=c, label=label, histtype=u'step')
    hep.cms.label()
    ax.set_xlabel(xlabel)
    ax.text(x=0.9, y=0.65, s='D and CB only', ha='right', transform=ax.transAxes)
    if distance is not None:
        ax.text(x=0.9, y=0.6, s='Matched vertices', ha='right', transform=ax.transAxes)
    ax.legend()
    fig.savefig(outName)
    print("Saved", outName)

def plotNTracks_Decay(nRecoTracks, decayStatus, distance, outName):
    possibleDecay = [0, 1, 2, 3]
    bins = np.arange(8)
    fig, ax = plt.subplots(1, 1)
    for pd in possibleDecay:
        mask = (decayStatus==pd) & (distance>0)
        c = np.histogram(np.clip(nRecoTracks[mask], bins[0], bins[-1]), bins)[0]
        c = c / np.sum(c)
        label = 'To other' if pd==0 else "To D-%d"%pd
        ax.hist(bins[:-1], bins=bins, weights=c, label=label, histtype=u'step')
    hep.cms.label()
    ax.set_xlabel("nRecoTracks Per SV")
    ax.text(x=0.9, y=0.65, s='B and BB only', ha='right', transform=ax.transAxes)
    ax.text(x=0.9, y=0.6, s='Matched vertices', ha='right', transform=ax.transAxes)
    ax.legend()
    fig.savefig(outName)
    print("Saved", outName)



def cumulativeDistance(distance, pdgClass, bins, labels, outFile):
    # Second plot:
    # Cumulative distribution as a function of the distance
    uniqueClasses = set(pdgClass)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    x = bins
    hep.cms.label()
    
    for id, myclass in enumerate(uniqueClasses):
        y2 = []
        mask = (pdgClass==myclass) & (distance>0)
        maskDen = (pdgClass==myclass) & (distance>0)
        label = class_to_labelForLegend[myclass]
        #if label not in labels: 
        #    continue

        for i in x:
            y2.append(np.sum(distance[mask]<i)/len(distance[maskDen]))
        
        linestyle = 'solid'
        color=None


        ax.plot(x, y2, label = label, color=color, linestyle=linestyle)
    ax.set_ylim(0., 1)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    sorted_handles_labels = sorted(zip(labels, handles), key=lambda x: x[0])
    labels, handles = zip(*sorted_handles_labels)
    #ax.legend(handles, labels, bbox_to_anchor=(1,1))
    ax.legend(handles, labels)

    ax.set_xlabel("distance [cm]")
    ax.set_ylabel("Fraction of matched within distance")
    fig.savefig(outFile, bbox_inches='tight')
    print("Saved ", outFile)

        


def main(nFiles):
    secondMatching = getSecondMatching()
    if secondMatching:
        rootDir="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul23Tuple/secondMatching/pt0p8" # change here
    else:
        #rootDir="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/TTToHadronic_withCutAtNano" # change here
        rootDir = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples"
        rootDir = "/work/gcelotto/BTV/scripts/tuplizer"
        #rootDir="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/flattuple/TTToHadronic2024Dec13Tuple_noCutSV_dlenSig" # change here
    fileNames = glob.glob(rootDir+"/*selPart_100.root")
    if nFiles == -1:
        pass
    elif len(fileNames)<nFiles:
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
            
            distance        = np.array(tree['distance'])
            pdgID           = np.array(tree['pdgID']).astype(int)
            eta             = np.array(tree['eta'])
            pt              = np.array(tree['pt'])
            displacement    = np.array(tree['displacement'])
            nGenTracks      = np.array(tree['nGenTracks'])
            decayStatus     = np.array(tree['decayStatus'])
            promptStatus    = np.array(tree['promptStatus'])
            nRecoTracks     = np.array(tree['nRecoTracks'])
            nSV             = np.array(tree['nSV'])
            delta_x             = np.array(tree['delta_x'])
            delta_y             = np.array(tree['delta_y'])
            delta_z             = np.array(tree['delta_z'])
        else:
            
            distance        = np.concatenate((distance, np.array(tree['distance'])))
            pdgID           = np.concatenate((pdgID, np.array(tree['pdgID']).astype(int)))
            eta             = np.concatenate((eta, np.array(tree['eta'])))
            pt              = np.concatenate((pt, np.array(tree['pt'])))
            displacement    = np.concatenate((displacement, tree['displacement']))
            nGenTracks      = np.concatenate((nGenTracks, tree['nGenTracks']))
            promptStatus    = np.concatenate((promptStatus, tree['promptStatus']))
            decayStatus     = np.concatenate((decayStatus, tree['decayStatus']))
            nRecoTracks     = np.concatenate((nRecoTracks, tree['nRecoTracks']))
            nSV             = np.concatenate((nSV, tree['nSV']))
            delta_x             = np.concatenate((delta_x, tree['delta_x']))
            delta_y             = np.concatenate((delta_y, tree['delta_y']))
            delta_z             = np.concatenate((delta_z, tree['delta_z']))
            
        #matched = np.array(tree['matched'])


    #print(distance[distance>0])
    mesons = map_to_groups_vec(pdgID)
    # B mes  : 0
    # D mes  : 1
    # S bar  : 2
    # C bar  : 3
    # B bar  : 4    

    commonMask      = displacement>-1
    mesons          = mesons[commonMask]
    distance        = distance[commonMask]
    pdgID           = pdgID[commonMask]
    eta             = eta[commonMask]
    pt              = pt[commonMask]
    displacement    = displacement[commonMask]
    nGenTracks      = nGenTracks[commonMask]
    promptStatus    = promptStatus[commonMask]
    decayStatus     = decayStatus[commonMask]
    nRecoTracks     = nRecoTracks[commonMask]
    nSV             = nSV[commonMask]
    delta_x         = delta_x[commonMask]
    delta_y         = delta_y[commonMask]
    delta_z         = delta_z[commonMask]

    print(np.sqrt(delta_x**2 + delta_y**2 + delta_z**2), distance)

    import pandas as pd
    df = pd.DataFrame({
        'distance':distance,
        'pdgClass':pdgID,
        'promptStatus':promptStatus
    })
    
    df['mesonsPrompt'] = df.apply(lambda row: map_to_groups_prompt(row['pdgClass'], row['promptStatus']), axis=1)

# PDG Efficiency
    efficiencyVsVariable(num=mesons[(distance>0)], den=mesons,
                        bins=np.linspace(-1.5, 4.5, 7), xlabel="Hadron PDG Class",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/all/eff_pdgDiff.png", title = "All particles",
                        tick_positions = [0, 1, 2, 3, 4],
                        tick_labels = ['B mesons', 'D mesons', 'Strange Baryons', 'Charmed Baryons', 'Bottom Baryons'],)
# Efficiency cumulative Distance
    cumulativeDistance(distance=distance, pdgClass=mesons, bins=np.logspace(np.log10(1e-4), np.log10(1), num=200), labels = ['B', 'D', 'BB', 'CB'], outFile="/work/gcelotto/BTV/plots/newTuple/efficiency/all/cumDistancePDG.png")

# Efficiency in B vs Final states
    efficiencyVsVariable(num=decayStatus[(distance>0)], den=decayStatus,
                        bins=np.linspace(-0.5, 3.5, 5), xlabel="DecayStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_BVsFinalStates.png", title = "B mesons and Bottom Baryons",
                        tick_positions = [ 0, 1, 2, 3],
                        tick_labels = [ 'B To Other', 'B To D', 'B To X To D', 'B To X To Y To D'])
    efficiencyVsVariable(num=decayStatus[(distance>0) & (mesons==0)], den=decayStatus[(mesons==0)],
                        bins=np.linspace(-0.5, 3.5, 5), xlabel="DecayStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_BmesVsFinalStates.png", title = "B mesons only",
                        tick_positions = [ 0, 1, 2, 3],
                        tick_labels = [ 'B To Other', 'B To D', 'B To X To D', 'B To X To Y To D'])
    efficiencyVsVariable(num=decayStatus[(distance>0) & (mesons==4)], den=decayStatus[(mesons==4)],
                        bins=np.linspace(-0.5, 3.5, 5), xlabel="DecayStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_BBaryonsVsFinalStates.png", title = "Bottom Baryons only",
                        tick_positions = [ 0, 1, 2, 3],
                        tick_labels = [ 'B To Other', 'B To D', 'B To X To D', 'B To X To Y To D'])
    # Do the same plot with only prompt non prompt
    efficiencyVsVariable(num=decayStatus[(distance>0) & (mesons==0)], den=decayStatus[(mesons==0)],
                        bins=np.linspace(-0.5, 1.5, 3), xlabel="DecayStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_BmesVsFinalStates_bool.png", title = "B mesons only",
                        tick_positions = [ 0, 1, ],
                        underOverFlow=True, 
                        tick_labels = [ 'B To Other', 'B To D'])
# Dmesons
    efficiencyVsVariable(num=promptStatus[(distance>0) & ((mesons==1) | (mesons==3))], den=promptStatus[((mesons==1) | (mesons==3))],
                        bins=np.linspace(-2.5, 3.5, 7), xlabel="PromptStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_D_PromptOrFromB_2.png", title = "D mesons and Charmed Baryons",
                        tick_positions = [-2,  0, 1, 2, 3],
                        tick_labels = [ 'From Other', 'Prompt', 'D From B' ,'D From X From B', 'D From X From Y From B'])
    efficiencyVsVariable(num=promptStatus[(distance>0) & (mesons==1)], den=promptStatus[(mesons==1)],
                        bins=np.linspace(-2.5, 3.5, 7), xlabel="PromptStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_Dmesons_PromptOrFromB.png", title = "D mesons only",
                        tick_positions = [-2,  0, 1, 2, 3],
                        tick_labels = [ 'From Other', 'Prompt', 'D From B' ,'D From X From B', 'D From X From Y From B'])
    efficiencyVsVariable(num=promptStatus[(distance>0) & (mesons==3)], den=promptStatus[(mesons==3)],
                        bins=np.linspace(-2.5, 3.5, 7), xlabel="PromptStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_DBaryons_PromptOrFromB.png", title = "Charmed Baryons only",
                        tick_positions = [-2,  0, 1, 2, 3],
                        tick_labels = [ 'From Other', 'Prompt', 'D From B' ,'D From X From B', 'D From X From Y From B'])
    
    efficiencyVsVariable(num=promptStatus[(distance>0) & (mesons==1)], den=promptStatus[(mesons==1)],
                        bins=np.linspace(-2.5, 1.5, 5), xlabel="PromptStatus",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_Dmesons_PromptOrFromB_bool.png", title = "D mesons only",
                        tick_positions = [-2,  0, 1],
                        underOverFlow=True,
                        tick_labels = [ 'From Other', 'Prompt', 'D From B'])
    

    efficiencyVsVariable(num=eta[(distance>0)], den=eta,
                        bins=np.linspace(-2.5, 2.5, 20), xlabel="Hadron eta",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/all/eff_etaDiff.png", title = "All particles",)

    efficiencyVsVariable(num=pt[(distance>0)], den=pt,
                        bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/all/eff_ptDiff.png", title = "All particles",)

    efficiencyVsVariable(num=displacement[(distance>0)], den=displacement,
                        bins=np.linspace(0, 0.2, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/all/eff_displacementDiff.png", title = "All particles",)

    efficiencyVsVariable(num=nGenTracks[(distance>0)], den=nGenTracks,
                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/all/eff_nGenTracksDiff.png", title = "All particles",)

    
    

#
#    # B only
#
    efficiencyVsVariable(num=eta[(distance>0) & (mesons==0)], den=eta[(mesons==0)],
                            bins=np.linspace(-2.5, 2.5, 20), xlabel="Hadron eta",
                            outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_etaDiff.png", title="B only",)
    efficiencyVsVariable(num=pt[(distance>0) & (mesons==0)], den=pt[(mesons==0)],
                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                            outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_ptDiff.png", title="B only",)
    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==0)], den=displacement[(mesons==0)],
                        bins=np.linspace(0, 1, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_displacementDiff.png", title = "B only",)
    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==0)], den=nGenTracks[(mesons==0)],
                        bins=np.arange(0, 10), xlabel="Number of hadron daughters",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_nGenTracksDiff.png", title = "B only",)

    
    
#
#
#    # D only
    efficiencyVsVariable(num=eta[(distance>0) & (mesons==1)], den=eta[(mesons==1)],
                            bins=np.linspace(-2.5, 2.5, 20), xlabel="Hadron eta",
                            outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_etaDiff.png", title="D only")
    efficiencyVsVariable(num=pt[(distance>0) & (mesons==1)], den=pt[(mesons==1)],
                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                            outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_ptDiff.png", title="D only")
    efficiencyVsVariable(num=pt[(distance>0) & (np.array(df.mesonsPrompt)==100)], den=pt[(np.array(df.mesonsPrompt)==100)],
                            bins=np.linspace(10, 150, 20), xlabel="Hadron pT [GeV]",
                            outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_nonPrompt_2_ptDiff.png", title="D np2 only")
    efficiencyVsVariable(num=displacement[(distance>0) & (mesons==1)], den=displacement[(mesons==1)],
                        bins=np.linspace(0, 0.5, 20), xlabel="Hadron distance of flight [cm]",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_displacementDiff.png", title = "D only",)
    efficiencyVsVariable(num=nGenTracks[(distance>0) & (mesons==1)], den=nGenTracks[(mesons==1)],
                        bins=np.arange(2, 10), xlabel="Number of hadron daughters",
                        outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/eff_nGenTracksDiff.png", title = "D only",)
    plotNTracks_Prompt(nRecoTracks, promptStatus, distance,  bins=np.arange(8), xlabel='nRecoTracks Per SV', outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/nRecoTracks_vs_PromptStatus.png")
    plotNTracks_Prompt(nSV, promptStatus, distance=None,  bins=np.arange(10), xlabel='nSV Per Event', outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/nSV_vs_PromptStatus.png")
    plotNTracks_Decay(nRecoTracks, decayStatus, distance,  outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Donly/nRecoTracks_vs_DecayStatus.png")
    distanceVsPDGClass(distance, pdgClass=mesons, bins=np.linspace(0, 2, 31), outFile="/work/gcelotto/BTV/plots/newTuple/efficiency/all/distance_vs_pdg.png", zoom=False)

    plt.close('all')
    newMask = (distance >-0)
    fig, ax = plt.subplots(1, 1)
    bins = np.linspace(-0.1, 0.1, 100)
    bins_z = np.linspace(-0.1, 0.1, 100)

    counts = np.histogram(np.clip(delta_x[newMask], bins[0], bins[-1]), bins)[0]
    ax.text(x=0.95, y=0.5, s="Mean : %.4f"%(np.mean(delta_x[newMask])), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.45, s="Err  : %.4f"%(np.std(delta_x[newMask])/np.sqrt(len(delta_x[newMask]))), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.4, s="Std  : %.2f"%(np.std(delta_x[newMask])), transform=ax.transAxes, ha='right')
    counts = counts/np.sum(counts)
    ax.hist(bins[:-1], bins, weights=counts)
    ax.set_xlabel("SVx - GVx [cm]")
    ax.set_ylabel("Normalized Events")
    hep.cms.label()
    fig.savefig("/work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_x.png", bbox_inches='tight')
    print("Saved /work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_x.png")

    ax.clear()
    counts = np.histogram(np.clip(delta_y[newMask], bins[0], bins[-1]), bins)[0]
    ax.text(x=0.95, y=0.5, s="Mean : %.4f"%(np.mean(delta_y[newMask])), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.45, s="Err  : %.4f"%(np.std(delta_y[newMask])/np.sqrt(len(delta_y[newMask]))), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.4, s="Std  : %.2f"%(np.std(delta_y[newMask])), transform=ax.transAxes, ha='right')
    counts = counts/np.sum(counts)
    ax.hist(bins[:-1], bins, weights=counts)
    ax.set_xlabel("SVy - GVy [cm]")
    ax.set_ylabel("Normalized Events")
    hep.cms.label()
    fig.savefig("/work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_y.png", bbox_inches='tight')
    print("Saved /work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_y.png")

    ax.clear()
    
    counts = np.histogram(np.clip(delta_z[newMask], bins_z[0], bins_z[-1]), bins_z)[0]
    ax.text(x=0.95, y=0.5, s="Mean : %.5f"%(np.mean(delta_z[newMask])), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.45, s="Err  : %.5f"%(np.std(delta_z[newMask])/np.sqrt(len(delta_z[newMask]))), transform=ax.transAxes, ha='right')
    ax.text(x=0.95, y=0.4, s="Std  : %.2f"%(np.std(delta_z[newMask])), transform=ax.transAxes, ha='right')
    counts = counts/np.sum(counts)
    ax.hist(bins_z[:-1], bins_z, weights=counts)
    ax.set_xlabel("SVz - GVz [cm]")
    ax.set_ylabel("Normalized Events")
    hep.cms.label()
    fig.savefig("/work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_z.png", bbox_inches='tight')
    print("Saved /work/gcelotto/BTV/plots/newTuple/efficiency/all/resolution_z.png")




    # eff pt with b mesons
    bins=np.linspace(10, 150, 20)
    xlabel="Hadron pT [GeV]"
    outName="/work/gcelotto/BTV/plots/newTuple/efficiency/Bonly/eff_pt_AllvsB.png"
    title="All particles vs B mesons"
    underOverFlow = False
    fig,ax = plt.subplots(1, 1)
    num=pt[(distance>0)]
    den=pt
    matched = np.histogram(num, bins=bins)[0]
    total =   np.histogram(den, bins=bins)[0]
    num_B=pt[(distance>0) & (mesons==0)]
    den_B=pt[mesons==0]
    matched_B = np.histogram(num_B, bins=bins)[0]
    total_B =   np.histogram(den_B, bins=bins)[0]

    num_D=pt[(distance>0) & (mesons==1)]
    den_D=pt[mesons==1]
    matched_D = np.histogram(num_D, bins=bins)[0]
    total_D =   np.histogram(den_D, bins=bins)[0]

    if underOverFlow:
        matched = np.histogram(np.clip(num, bins[0], bins[-1]), bins=bins)[0]
        total =   np.histogram(np.clip(den, bins[0], bins[-1]), bins=bins)[0]


    bins_center = (bins[:-1]+bins[1:])/2
    ax.errorbar(bins_center, matched/(total+eps), xerr=np.diff(bins)/2,yerr=np.sqrt(matched)/(total+eps), marker='o', linestyle='none', color='black', label='All particles')
    ax.errorbar(bins_center, matched_B/(total_B+eps), xerr=np.diff(bins)/2,yerr=np.sqrt(matched_B)/(total_B+eps), marker='o', linestyle='none', color='red', label='B mesons only')
    ax.errorbar(bins_center, matched_D/(total_D+eps), xerr=np.diff(bins)/2,yerr=np.sqrt(matched_D)/(total_D+eps), marker='o', linestyle='none', color='blue', label='D mesons only')
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Efficiency [%]")
    ax.set_ylim(0, 1.2)
    ax.text(x=0.95, y=0.94, s="Matched Entries %d"%(len(num)), transform=ax.transAxes, horizontalalignment='right')
    #ax.text(x=0.95, y=0.89, s=title, transform=ax.transAxes, horizontalalignment='right')
    #print(tick_positions)
    #if tick_positions is not None:
    #    ax.set_xticks(tick_positions, tick_labels,  rotation=30)
    #    for x,y in zip(bins_center, matched/(total+eps)):
    #        if y!=0:
    #            ax.text(x=x, y=y+0.05, s="%.1f%%"%(y*100), ha='center')

    #if title is not None:
    hep.cms.label()
    fig.savefig(outName, bbox_inches='tight')
    print("Saved %s"%outName)
    plt.close()
    
    

    flightMatchedDistance(distance=distance, displacement=displacement, bins=np.linspace(0, .3, 31),outFile="/work/gcelotto/BTV/plots/newTuple/efficiency/all/displacementM_vs_NM.png")

    
    
    print("***********\nYou should add an option for the overflow True False in the efficiencyVsVariable function\n****************")
    print("Overall efficiency: %.2f"%(np.sum(distance>0)/len(distance)*100))
    return 0


if __name__=="__main__":
    nFiles = int(sys.argv[1]) if len(sys.argv)>1 else 1
    main(nFiles=nFiles)

