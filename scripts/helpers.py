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
def getTreeAndBranches(fileName):
    print("Opening ", fileName)
    f = uproot.open(fileName)
    tree = f['Events']
    branches = tree.arrays()
    return tree, branches
def eventDisplay(df_event, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy):
    fig, ax = plt.subplots(1, 1)
    
    ax.scatter(df_event.vx, df_event.vy, label='GenVertices')

    ax.text(x=1.02, y=0.45, s="x", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.08, y=0.45, s="y", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.14, y=0.45, s="z", ha='center', fontsize=14, transform=ax.transAxes)
    ax.text(x=1.2, y=0.45, s="PDG", ha='center', fontsize=14, transform=ax.transAxes)
    y=0.4
    for row in range(len(df_event)):
        ax.text(x=df_event.vx.iloc[row], y=df_event.vy.iloc[row], s=df_event.pdgClass.iloc[row], fontsize=12)
        ax.text(x=1.02, y=y, s="%.2f"%(df_event.vx.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.08, y=y, s="%.2f"%(df_event.vy.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.14, y=y, s="%.2f"%(df_event.vz.iloc[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.2, y=y, s="%s"%(df_event.pdgClass[row] ), ha='center', fontsize=14, transform=ax.transAxes)
        y=y-0.05
        if GenPart_genPartIdxMother[df_event.idx.iloc[row]]!=-1:
            ax.plot([df_event.vx.iloc[row], GenPart_vx[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]], [df_event.vy.iloc[row], GenPart_vy[GenPart_genPartIdxMother[df_event.idx.iloc[row]]]],linestyle='dotted',alpha=0.8, marker='none')
    if SVs.shape[0]>0:
        ax.text(x=1.02, y=0.95, s="x", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.08, y=0.95, s="y", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.14, y=0.95, s="z", ha='center', fontsize=14, transform=ax.transAxes)
        ax.text(x=1.2, y=0.95, s="chi2", ha='center', fontsize=14, transform=ax.transAxes)
        y=0.9
        for row in range(len(SVs)):
            ax.text(x=1.02, y=y, s="%.2f"%(SVs[row][0]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.08, y=y, s="%.2f"%(SVs[row][1]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.14, y=y, s="%.2f"%(SVs[row][2]), ha='center', fontsize=14, transform=ax.transAxes)
            ax.text(x=1.20, y=y, s="%.2f"%(SV_chi2[row]), ha='center', fontsize=14, transform=ax.transAxes)
            y=y-0.05

        ax.scatter(SVs[:,0], SVs[:,1], label='Reco SV', marker="s", color='C1')
    else:
        pass
        assert len(SVs)==0
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")

    ax.scatter(PV_x, PV_y, label='PV', color='red')


    theta = np.linspace(0, 2*np.pi, 100)
    x = PV_x + 0.1 * np.cos(theta)
    y = PV_y + 0.1 * np.sin(theta)
    ax.plot(x, y, color='C0')
    # mathcing
    ax.legend()
    return fig, ax


def criterion0_new(distances, rows, tracksCounters, ProbeTracks_pt):

    while (np.any(distances < 997)):
        recoIdx, genIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)[:2]
        if distances[recoIdx, genIdx]*10 - np.array(rows[genIdx]['displacement']) > 5:
            # dont match the gen part
            distances[recoIdx, genIdx]=998
            continue
        # put rows distance and matched
        rows[genIdx]['distance']=distances[recoIdx, genIdx]
        rows[genIdx]['matched']=True

        # instead of deleting matched vertices, replace row and columns with high values
        distances[recoIdx, :]=[998]*distances.shape[1]
        distances[:, genIdx]=[998]*distances.shape[0]
    return rows


def criterion0(distances, distances_normalized, df_event, display, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy, tracksCounters, ProbeTracks_pt,
ProbeTracks_matchedToSV, col_mask, row_mask):
    if display:
        fig, ax = eventDisplay(df_event=df_event, SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx, GenPart_vy=GenPart_vy)
    
    while (np.any(distances < 997)):
        minIdx = np.unravel_index(np.argmin(distances, axis=None), distances.shape)
        #try:

        if distances[minIdx[0], minIdx[1]]*10 - np.array(df_event.displacement[col_mask])[minIdx[1]] > 5:
            # dont match the gen part
            distances[minIdx[0], minIdx[1]]=998
            continue
        df_event.loc[minIdx[1], 'matched']=True
        df_event.loc[minIdx[1], 'distance']=distances[minIdx[0], minIdx[1]]
        df_event.loc[minIdx[1], 'normDistance']=distances_normalized[minIdx[0], minIdx[1]]
        df_event.loc[minIdx[1], 'probeTracksFromSV']=tracksCounters[minIdx[0]]
        
        
        #df_event.loc[minIdx[1], 'probeTracks_pt'] = probePt if len(probePt)>0 else [None]
        
        if display:
            x_values = [SVs[row_mask][minIdx[0]][0], np.array(df_event.vx)[col_mask][minIdx[1]]]
            y_values = [SVs[row_mask][minIdx[0]][1], np.array(df_event.vy)[col_mask][minIdx[1]]]
            #plot matching
            ax.plot(x_values, y_values, color='black', marker='none')

        # instead of deleting matched vertices, replace row and columns with high values
        distances[minIdx[0], :]=[998]*distances.shape[1]
        distances[:, minIdx[1]]=[998]*distances.shape[0]

        distances_normalized[minIdx[0], :] = [998]*distances.shape[1]
        distances_normalized[:, minIdx[1]] = [998]*distances.shape[0]
    if display:
        hep.cms.label()
        fig.savefig("/t3home/gcelotto/BTV/plots/EventGenReco.png", bbox_inches='tight')
        plt.close('all')
    return df_event

def criterion1(distances, distances_normalized, df_event, display, SVs, SV_chi2, PV_x, PV_y, GenPart_genPartIdxMother, GenPart_vx, GenPart_vy, col_mask, row_mask):
    if display:
        fig, ax = eventDisplay(df_event=df_event, SVs=SVs, SV_chi2=SV_chi2, PV_x=PV_x, PV_y=PV_y, GenPart_genPartIdxMother=GenPart_genPartIdxMother, GenPart_vx=GenPart_vx, GenPart_vy=GenPart_vy)
    row_ind, col_ind = linear_sum_assignment(distances)
    for i in range(len(row_ind)):
        df_event.loc[col_ind[i], 'matched']=True
        df_event.loc[col_ind[i], 'distance']=distances[row_ind[i], col_ind[i]]
        df_event.loc[col_ind[i], 'normDistance']=distances_normalized[row_ind[i], col_ind[i]]
        x_values = [SVs[row_mask][row_ind[i]][0], np.array(df_event.vx)[col_mask][col_ind[i]]]
        y_values = [SVs[row_mask][row_ind[i]][1], np.array(df_event.vy)[col_mask][col_ind[i]]]
        ax.plot(x_values, y_values, color='black', marker='none')
    return df_event