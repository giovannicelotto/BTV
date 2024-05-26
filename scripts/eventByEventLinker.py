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
from tuplizer.utilsForScript import distance_3d, getPdgMask
from helpers import getTreeAndBranches, criterion0, criterion1, eventDisplay, getTreeAndBranches


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



def main():

    tree, branches = getTreeAndBranches(fileName = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/TTToHadronic_2_Run2_mc_124X.root")
    nEvents = tree.num_entries

    for ev in np.arange(nEvents):
        input("Next ev %d"%ev)
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
        ProbeTracks_matchedToSV     = branches["ProbeTracks_matchedToSV"][ev]
        ProbeTracks_pt              = branches["ProbeTracks_pt"][ev]
        ProbeTracks_genPartIdx      = branches["ProbeTracks_genPartIdx"][ev]

        SVs = np.array([(x, y, z) for x, y, z in zip(SV_x, SV_y, SV_z)])

        pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId)
        etaMask = abs(GenPart_eta)<2.5
        ptMask = GenPart_pt>10
        mesons = np.arange(nGenPart)[pdgMask & etaMask & ptMask]


        oneDaughter = []      # one index of a daughter of the meson aligned
        for mes in mesons:
            foundDaughter = -1
            for gp in range(nGenPart):
                if (GenPart_genPartIdxMother[gp] == mes):
                    oneDaughter.append(gp)
                    foundDaughter=gp
                    # fill one daughter per meson
                    break
            if foundDaughter==-1:
                oneDaughter.append(-1)
        assert len(mesons)==len(oneDaughter)

        genVertices = np.array([(x, y, z) for x, y, z in zip(GenPart_vx[oneDaughter], GenPart_vy[oneDaughter], GenPart_vz[oneDaughter])])
        allDaughters = []      # one index of a daughter of the meson aligned
        for mes in mesons:
            for gp in range(nGenPart):
                if (GenPart_genPartIdxMother[gp] == mes):
                    allDaughters.append(gp)


        distances = np.array([[distance_3d(sv, (vx, vy, vz)) for vx, vy, vz in genVertices] for sv in SVs ])
        distances_filled = distances.copy()
        matchingKey ={}
        while (np.any(distances_filled < 997)):
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            #print(recoIdx, genIdx, " possibly matched")
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)

            #tracks from sv
            recoTracksMask = (ProbeTracks_matchedToSV==recoIdx)

            #tracks from gv
            genTracksMask = (GenPart_genPartIdxMother==GenPart_genPartIdxMother[oneDaughter[genIdx]])
            genTracksMask = genTracksMask | (GenPart_genPartIdxMother==GenPart_genPartIdxMother[GenPart_genPartIdxMother[oneDaughter[genIdx]]])
            genTracksMask = genTracksMask | (GenPart_genPartIdxMother==GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[oneDaughter[genIdx]]]])
            #check if genTracks linked to recoTracks are in genVertex
            commonTracks = np.sum(np.in1d(ProbeTracks_genPartIdx[recoTracksMask], np.arange(nGenPart)[genTracksMask]))
            if commonTracks>=1:
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                matchingKey[genIdx]=recoIdx
            else:
                distances_filled[recoIdx, genIdx]=998
        

        # matching is done

        fig, ax = plt.subplots(1, 1)
        ax.scatter(GenPart_vx[oneDaughter], GenPart_vy[oneDaughter], label='GenVertices', s=80)
        for gp in range(len(oneDaughter)):
            ax.text(x=GenPart_vx[oneDaughter][gp], y=GenPart_vy[oneDaughter][gp], s=map_to_groups_letter(GenPart_pdgId[mesons][gp]), fontsize=18)
        print(len(oneDaughter))
        for dau in oneDaughter:

            if GenPart_genPartIdxMother[dau]!=-1: # the meson
                x = [GenPart_vx[GenPart_genPartIdxMother[dau]], GenPart_vx[dau]]
                y = [GenPart_vy[GenPart_genPartIdxMother[dau]], GenPart_vy[dau]]
                #print(x, y)
                ax.plot(x, y,linestyle='dotted',alpha=0.8, marker='none')
            else:
                print("-1 found")

        #ax.text(x=1.02, y=0.45, s="x", ha='center', fontsize=14, transform=ax.transAxes)
        #ax.text(x=1.08, y=0.45, s="y", ha='center', fontsize=14, transform=ax.transAxes)
        #ax.text(x=1.14, y=0.45, s="z", ha='center', fontsize=14, transform=ax.transAxes)
        #ax.text(x=1.2, y=0.45, s="PDG", ha='center', fontsize=14, transform=ax.transAxes)
        #y=0.4
        #for row in range(len(df_event)):
        #    ax.text(x=df_event.vx.iloc[row], y=df_event.vy.iloc[row], s=df_event.pdgClass.iloc[row], fontsize=12)
        #    ax.text(x=1.02, y=y, s="%.2f"%(df_event.vx.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.08, y=y, s="%.2f"%(df_event.vy.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.14, y=y, s="%.2f"%(df_event.vz.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.2, y=y, s="%s"%(df_event.pdgClass[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        #    y=y-0.05
        #    if GenPart_genPartIdxMother[df_event.idx.iloc[row]]!=-1:
        #        ax.plot([df_event.vx.iloc[row], GenPart_vx[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]], [df_event.vy.iloc[row], GenPart_vy[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]],linestyle='dotted',alpha=0.8, marker='none')
        if SVs.shape[0]>0:
        #    ax.text(x=1.02, y=0.95, s="x", ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.08, y=0.95, s="y", ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.14, y=0.95, s="z", ha='center', fontsize=14, transform=ax.transAxes)
        #    ax.text(x=1.2, y=0.95, s="chi2", ha='center', fontsize=14, transform=ax.transAxes)
        #    y=0.9
        #    for row in range(len(SVs)):
        #        ax.text(x=1.02, y=y, s="%.2f"%(SVs[row][0]), ha='center', fontsize=14, transform=ax.transAxes)
        #        ax.text(x=1.08, y=y, s="%.2f"%(SVs[row][1]), ha='center', fontsize=14, transform=ax.transAxes)
        #        ax.text(x=1.14, y=y, s="%.2f"%(SVs[row][2]), ha='center', fontsize=14, transform=ax.transAxes)
        #        ax.text(x=1.20, y=y, s="%.2f"%(SV_chi2[row]), ha='center', fontsize=14, transform=ax.transAxes)
        #        y=y-0.05
#
            ax.scatter(SVs[:,0], SVs[:,1], label='Reco SV', marker="s", color='C1', s=80)
        else:
            pass
            assert len(SVs)==0
        ax.set_xlabel("X [cm]")
        ax.set_ylabel("Y [cm]")

        ax.scatter(PV_x, PV_y, label='PV', color='red')
        for genIdx in matchingKey:
            x_values = [SVs[matchingKey[genIdx]][0], GenPart_vx[oneDaughter][genIdx]]
            y_values = [SVs[matchingKey[genIdx]][1], GenPart_vy[oneDaughter][genIdx]]
            ax.plot(x_values, y_values, color='black', marker='none')


        #theta = np.linspace(0, 2*np.pi, 100)
        #x = PV_x + 0.1 * np.cos(theta)
        #y = PV_y + 0.1 * np.sin(theta)
        #ax.plot(x, y, color='C0')
        # mathcing
        ax.legend()
        fig.savefig("/t3home/gcelotto/BTV/plots/newEventDisplay.png", bbox_inches='tight')
        plt.close()
        print(matchingKey)


if __name__=="__main__":
    main()


        