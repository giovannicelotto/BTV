import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
import sys
from BTV.scripts.plotScripts.old.criterion1_differential import getSuffix
eps = 0.000001
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
def plot2d(x, y, x_bins, y_bins, xlabel, ylabel, outName, title):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig, ax_main = plt.subplots(figsize=(8, 8))
    divider = make_axes_locatable(ax_main)
    ax_top = divider.append_axes("top", 1.2, pad=0.2, sharex=ax_main)
    ax_right = divider.append_axes("right", 1.2, pad=0.2, sharey=ax_main)

# Plot the 2D histogram in the main axes
    hist, x_edges, y_edges = np.histogram2d(x=np.clip(x, x_bins[0], x_bins[-1]), y=np.clip(y, y_bins[0], y_bins[-1]), bins=[x_bins, y_bins])
    import matplotlib.colors as mcolors
    ax_main.imshow(hist.T, origin='lower', extent=(x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()), aspect='auto', cmap='viridis', norm=mcolors.LogNorm())
    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)

    # Plot the marginalized histogram on top
    ax_top.hist(np.clip(x, x_bins[0], x_bins[-1]), bins=x_bins, color='lightblue', edgecolor='black')
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_yticks([])
    ax_top.xaxis.tick_top()

    # Plot the marginalized histogram on the right
    ax_right.hist(np.clip(y, y_bins[0], y_bins[-1]), bins=y_bins, color='lightblue', edgecolor='black', orientation='horizontal')#lightcoral
    ax_right.set_ylim(ax_main.get_ylim())
    ax_right.set_xticks([])
    ax_right.yaxis.tick_right()

    ax_top.set_title(title)
    ax_main.plot(x_bins, (2*x_bins+1)/5, linestyle='dotted', color='red')

    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)
    return





def plotPtEtaMatchedVsNonMatched(df, folder, suffix,title):
    #Bmask = (df.pdgID == 511) | (df.pdgID == 521) | (df.pdgID == 531) | (df.pdgID == 541)
    matchedMask = df.matched == True
    
    fig, ax = plt.subplots(2, 1, constrained_layout=True)
    bins = np.linspace(10, 150, 15)
    cmatched = np.histogram(df[matchedMask].pt, bins=bins)[0]
    cnonmatched = np.histogram(df[~matchedMask].pt, bins=bins)[0]
    cmatched, cnonmatched = cmatched/np.sum(cmatched), cnonmatched/np.sum(cnonmatched)
    ax[0].hist(bins[:-1], bins=bins, weights=cmatched, histtype=u'step', label='matched')
    ax[0].hist(bins[:-1], bins=bins, weights=cnonmatched, histtype=u'step', label='non matched')
    ax[0].legend()
    ax[0].set_xlabel("GenPart_pt [GeV]")
    if title is not None:
        ax[0].set_title(title)
    bins=np.linspace(-3, 3, 20)
    cmatched            = np.histogram(df[matchedMask].eta, bins=bins)[0]
    cnonmatched         = np.histogram(df[~matchedMask].eta, bins=bins)[0]
    cmatched, cnonmatched = cmatched/np.sum(cmatched), cnonmatched/np.sum(cnonmatched)
    ax[1].hist(bins[:-1], bins=bins, weights=cmatched, histtype=u'step', label='matched')
    ax[1].hist(bins[:-1], bins=bins, weights=cnonmatched, histtype=u'step', label='non matched')
    ax[1].legend()
    ax[1].set_xlabel("GenPart_eta ")
    fig.savefig(folder + "/ptEta_%s.png"%suffix, bbox_inches='tight')
    print("Saved ", folder + "/ptEta_%s.png"%suffix)
    
    plt.close('all')


    fig, ax = plt.subplots(1, 1)
    ax.scatter(df[matchedMask].vx, df[matchedMask].vy, label='matched', color='blue')
    ax.scatter(df[~matchedMask].vx, df[~matchedMask].vy, label='non matched', color='red')
    ax.set_xlabel("GenPart_Vx [cm]")
    ax.set_ylabel("GenPart_Vy [cm]")

    ax.legend(frameon=True, edgecolor='black')
    fig.savefig(folder + "/VxVy_%s.png"%suffix)
    print(folder + "/VxVy_%s.png"%suffix)

    ax.clear()
    bins=np.linspace(0, 0.5, 20)
    cmatched            = np.histogram(df[matchedMask].displacement, bins=bins)[0]
    cnonmatched         = np.histogram(df[~matchedMask].displacement, bins=bins)[0]
    cmatched, cnonmatched = cmatched/np.sum(cmatched), cnonmatched/np.sum(cnonmatched)
    ax.hist(bins[:-1], bins=bins, weights=cmatched, histtype=u'step', label='matched')
    ax.hist(bins[:-1], bins=bins, weights=cnonmatched, histtype=u'step', label='non matched')
    ax.legend()
    ax.set_xlabel("GenPart_displacement ")
    fig.savefig(folder + "/displacement_%s.png"%suffix, bbox_inches='tight')
    print(folder + "/displacement_%s.png"%suffix)


    #ax.clear()
    #bins=np.linspace(0, 900, 20)
    #cmatched            = np.histogram(np.clip(df[matchedMask & Bmask].displacement/df[matchedMask & Bmask].distance, bins[0], bins[-1]), bins=bins)[0]
    #ax.hist(bins[:-1], bins=bins, weights=cmatched, histtype=u'step', label='matched')
    #ax.set_yscale('log')
    #ax.legend()
    #ax.set_xlabel("GenMatch_distance / GenPart_displacement ")
    #fig.savefig("/t3home/gcelotto/BTV/plots/displacementDist_%s.png"%suffix, bbox_inches='tight')
    #print("/t3home/gcelotto/BTV/plots/displacementDist_%s.png"%suffix)
    #ax.clear()
    
    
    
def main(criterion, threshold):
    suffix=getSuffix(criterion, threshold)
    df = pd.read_parquet("/t3home/gcelotto/BTV/output/df_"+suffix+".parquet")
    df['mesons'] = df['pdgID'].apply(map_to_groups)
    df['daughters_pt'] = df['daughters_pt'].apply(lambda x: list(map(float, x.split(','))) if len(x)>0 else None)
    df['daughters_pdgId'] = df['daughters_pdgId'].apply(lambda x: list(map(int, x.split(','))) if len(x)>0 else None)
    df['genDaughters'] = df['daughters_pt'].apply(lambda x: len(x) if x is not None else 0 )

    # example of how to get mask of the daughters properties
    df['genDaughters_pt30Mask'] = df['daughters_pt'].apply(lambda x: np.array(x)>30 if x is not None else None )
    df['genDaughters_pdgMask'] = df['daughters_pdgId'].apply(lambda x: np.isin(x, [111, 211, -211]) if x is not None else None )
    
    
    import os
    folders = ["/t3home/gcelotto/BTV/plots/%s"%suffix,
                "/t3home/gcelotto/BTV/plots/%s/matchedVsNonMatched"%suffix,
                "/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Bonly",
                "/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Donly",
                "/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/all",
                "/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Bonly",
                "/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Donly",
                "/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/all"]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
        
    if True:
        fig,ax = plt.subplots(nrows=1, ncols=1, sharex=True, layout='constrained')
        bins = np.linspace(-5.5, 5.5, 12)
        c = np.histogram(df[(df.matched==True) & (df.mesons==0)].genDaughters - df[(df.matched==True) & (df.mesons==0)].probeTracksFromSV, bins=bins)[0]
        ax.hist(bins[:-1], bins=bins, weights=c, histtype=u'step', linewidth=4, color='black', label='Inclusive')[0]
        close = np.histogram(df[(df.matched==True) & (df.mesons==0) & (df.distance<0.01)].genDaughters - df[(df.matched==True) & (df.mesons==0) & (df.distance<0.01)].probeTracksFromSV, bins=bins)[0]
        close = ax.hist(bins[:-1], bins=bins, weights=close, label='0.00 < d < 0.01')[0]
        medium = np.histogram(df[(df.matched==True) & (df.mesons==0) & (df.distance>0.01) & (df.distance<0.1) ].genDaughters - df[(df.matched==True) & (df.mesons==0) & (df.distance>0.01) & (df.distance<0.1)].probeTracksFromSV, bins=bins)[0]
        medium = ax.hist(bins[:-1], bins=bins, weights=medium, bottom=close, label='0.01 < d < 0.10')[0]
        high = np.histogram(df[(df.matched==True) & (df.mesons==0) & (df.distance>0.1) ].genDaughters - df[(df.matched==True) & (df.mesons==0) & (df.distance>0.1)].probeTracksFromSV, bins=bins)[0]
        high = ax.hist(bins[:-1], bins=bins, weights=high, bottom=(close+medium), label='0.10 < d < Inf')[0]
        ax.set_xlabel("#GenDaughters - #ProbeMatched")
        ax.set_ylabel("")
        ax.legend()
        outName = "/t3home/gcelotto/BTV/plots/cr0_t-1/temp.png"
        hep.cms.label()
        fig.savefig(outName, bbox_inches='tight')
        print("Saved %s"%outName)

        ax.clear()
        bins = np.linspace(0, 0.2, 30)
        for delta in [0, 1, 2]:
            m = (df[(df.matched==True) & (df.mesons==0)].genDaughters - df[(df.matched==True) & (df.mesons==0)].probeTracksFromSV)==delta
            c = np.histogram(np.clip(df[(df.matched==True) & (df.mesons==0)][m].distance, bins[0], bins[-1]), bins=bins, )[0]
            cerr = np.sqrt(c)
            c=c/np.sum(c)
            cerr= cerr/np.sum(c)
            ax.hist(bins[:-1], bins=bins, weights=c, histtype=u'step',linewidth=2, label='Delta NTracks %d'%delta)[0]
        #ax.set_yscale('log')
        ax.set_xlabel("Distance Gen-Reco matching")
        ax.set_ylabel("")
        ax.legend()
        outName = "/t3home/gcelotto/BTV/plots/cr0_t-1/temp_2.png"
        hep.cms.label()
        fig.savefig(outName, bbox_inches='tight')
        print("Saved %s"%outName)

        ax.clear()
        bins = np.linspace(0, 5, 30)
        for delta in [0, 1, 2]:
            m = (df[(df.matched==True) & (df.mesons==0)].genDaughters - df[(df.matched==True) & (df.mesons==0)].probeTracksFromSV)==delta
            c = np.histogram(np.clip(df[(df.matched==True) & (df.mesons==0)][m].displacement, bins[0], bins[-1]), bins=bins, )[0]
            cerr = np.sqrt(c)
            c=c/np.sum(c)
            cerr= cerr/np.sum(c)
            ax.hist(bins[:-1], bins=bins, weights=c, histtype=u'step',linewidth=2, label='Delta NTracks %d'%delta)[0]
        #ax.set_yscale('log')
        ax.set_xlabel("GenPart Travel")
        ax.set_ylabel("")
        ax.legend()
        outName = "/t3home/gcelotto/BTV/plots/cr0_t-1/temp_3.png"
        hep.cms.label()
        fig.savefig(outName, bbox_inches='tight')
        print("Saved %s"%outName)

    plot2d(x=df[df.matched == True].displacement, y=df[df.matched == True].distance,
           x_bins=np.linspace(0, 10, 20), y_bins=np.linspace(0, 5, 20),
           xlabel="displacement [cm]", ylabel="GenMatch_distance [cm]",
           outName="/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/all/displacementVsDistance_%s.png"%suffix,
           title = "All Particles")
    
    plot2d(x=df[(df.matched == True) & (df.mesons==0)].displacement, y=df[(df.matched == True) & (df.mesons==0)].distance,
           x_bins=np.linspace(0, 10, 20), y_bins=np.linspace(0, 5, 20),
           xlabel="displacement [cm]", ylabel="GenMatch_distance [cm]",
           outName="/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Bonly/displacementVsDistance_%s.png"%suffix,
           title = "B only")
    plot2d(x=df[(df.matched == True) & (df.mesons==1)].displacement, y=df[(df.matched == True) & (df.mesons==1)].distance,
           x_bins=np.linspace(0, 10, 20), y_bins=np.linspace(0, 5, 20),
           xlabel="displacement [cm]", ylabel="GenMatch_distance [cm]",
           outName="/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Donly/displacementVsDistance_%s.png"%suffix,
           title = "D only")

    efficiencyVsVariable(num=df.eta[df.matched==True], den=df.eta,
                         bins=np.linspace(-2.7, 2.7, 20), xlabel="GenPart_eta",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/all/etaDiff_"+suffix+".png", title = "All particles",)
    efficiencyVsVariable(num=df.pt[df.matched==True], den=df.pt,
                         bins=np.linspace(10, 150, 20), xlabel="GenPart_pt",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/all/ptDiff_"+suffix+".png", title = "All particles",)
    efficiencyVsVariable(num=df.mesons[df.matched==True], den=df.mesons,
                         bins=np.linspace(-1.5, 4.5, 7), xlabel="GenPart_pdgID",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/all/pdgDiff_"+suffix+".png",
                        tick_positions = [0, 1, 2, 3, 4],
                        tick_labels = ['B mesons', 'D mesons', 'Strange Baryons', 'Charmed Baryons', 'Bottom Baryons'],
                        )
    # B only
    efficiencyVsVariable(num=df.eta[(df.matched==True) & (df.mesons==0)], den=df.eta[df.mesons==0],
                         bins=np.linspace(-2.7, 2.7, 20), xlabel="GenPart_eta",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Bonly/etaDiff_"+suffix+".png", title="B only",)
    efficiencyVsVariable(num=df.pt[(df.matched==True) & (df.mesons==0)], den=df.pt[df.mesons==0],
                         bins=np.linspace(10, 150, 20), xlabel="GenPart_pt",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Bonly/ptDiff_"+suffix+".png", title="B only",)
    plotPtEtaMatchedVsNonMatched(df[df.mesons==0], folder="/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Bonly", suffix=suffix, title="B only")

    # D only
    efficiencyVsVariable(num=df.eta[(df.matched==True) & (df.mesons==1)], den=df.eta[df.mesons==1],
                         bins=np.linspace(-2.7, 2.7, 20), xlabel="GenPart_eta",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Donly/etaDiff_"+suffix+".png", title="D only")
    efficiencyVsVariable(num=df.pt[(df.matched==True) & (df.mesons==1)], den=df.pt[df.mesons==1],
                         bins=np.linspace(10, 150, 20), xlabel="GenPart_pt",
                         outName="/t3home/gcelotto/BTV/plots/"+suffix+"/efficiency/Donly/ptDiff_"+suffix+".png", title="D only")
    plotPtEtaMatchedVsNonMatched(df[df.mesons==1], folder="/t3home/gcelotto/BTV/plots/"+suffix+"/matchedVsNonMatched/Donly", suffix=suffix, title="D only")




    
    
    



if __name__ =="__main__":
    criterion=int(sys.argv[1])
    t=int(sys.argv[2])
    print("Threshold = %d"%t)
    main(criterion=criterion, threshold=t)
    