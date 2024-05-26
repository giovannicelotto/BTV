import uproot
import glob
import matplotlib.pyplot as plt
import numpy as np
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mplhep as hep
from matplotlib.colors import LogNorm
hep.style.use("CMS")
path ="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06"
fileNames = glob.glob(path+"/**/*.root", recursive=True)
X,Y,Z = [],[],[]
for fileName in fileNames[:1]:
    print("opening ", fileName)
    f = uproot.open(fileName)
    tree = f['Events']
    branches = tree.arrays()
    maxEntries = tree.num_entries
    #for ev in range(maxEntries):
    ProbeTracks_genPartIdx  = branches["ProbeTracks_genPartIdx"]
    ProbeTracks_pt          = branches["ProbeTracks_pt"]
    ProbeTracks_eta          = branches["ProbeTracks_eta"]
    ProbeTracks_phi          = branches["ProbeTracks_phi"]
    GenPart_pt              = branches["GenPart_pt"]
    GenPart_eta              = branches["GenPart_eta"]
    GenPart_phi              = branches["GenPart_phi"]

    matched = ProbeTracks_genPartIdx>-1
    y = -(GenPart_pt[ProbeTracks_genPartIdx[matched]] - ProbeTracks_pt[matched])/GenPart_pt[ProbeTracks_genPartIdx[matched]]
    z = np.sqrt((GenPart_eta[ProbeTracks_genPartIdx[matched]] - ProbeTracks_eta[matched])**2 + (GenPart_phi[ProbeTracks_genPartIdx[matched]] - ProbeTracks_phi[matched])**2)
    x = GenPart_pt[ProbeTracks_genPartIdx[matched]]
    y = list(itertools.chain.from_iterable(y))
    z = list(itertools.chain.from_iterable(z))
    x = list(itertools.chain.from_iterable(x))
    Y.append(y)
    X.append(x)
    Z.append(z)

X = list(itertools.chain.from_iterable(X))
Y = list(itertools.chain.from_iterable(Y))
Z= list(itertools.chain.from_iterable(Z))

fig, ax = plt.subplots(1, 1)
ax.hist(Y, bins=np.linspace(-.25, .25, 100))
ax.set_xlabel(r"$(pT_{Reco} - pT_{Gen})/pT_{Gen}$")
fig.savefig("/t3home/gcelotto/BTV/plots/pTResolution.png", bbox_inches='tight')

fig, ax = plt.subplots(1, 1)
ax.hist(Z, bins=np.linspace(0, .2, 100))
ax.set_xlabel(r"$(\Delta R (Reco,Gen)$")
fig.savefig("/t3home/gcelotto/BTV/plots/dRResolution.png", bbox_inches='tight')

x_bins, y_bins = np.linspace(0, 25, 15), np.linspace(-.5, .5, 20)
fig, ax_main = plt.subplots(figsize=(8, 8))
divider = make_axes_locatable(ax_main)
ax_top = divider.append_axes("top", 1.2, pad=0.2, sharex=ax_main)
ax_right = divider.append_axes("right", 1.2, pad=0.2, sharey=ax_main)

# Plot the 2D histogram in the main axes
hist, x_edges, y_edges = np.histogram2d(x=X, y=Y, bins=[x_bins, y_bins])
img = ax_main.imshow(hist.T, origin='lower', extent=(x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()), aspect='auto', cmap='viridis', norm=LogNorm(vmin=100))
cax = divider.append_axes("right", size="5%", pad=.8)
cbar = plt.colorbar(img, cax=cax)
cbar.set_label('Counts')
#ax_main.scatter(X, Y, s=2)
ax_main.set_xlim(x_bins[0],x_bins[-1])
ax_main.set_xlabel("Track pT [GeV]")
ax_main.set_ylabel("Resolution")

# Plot the marginalized histogram on top
ax_top.hist(X, bins=x_bins, color='lightblue', edgecolor='black')
ax_top.set_xlim(ax_main.get_xlim())
ax_top.set_yticks([])
ax_top.xaxis.tick_top()

# Plot the marginalized histogram on the right
ax_right.hist(Y, bins=y_bins, color='lightblue', edgecolor='black', orientation='horizontal')#lightcoral
ax_right.set_ylim(ax_main.get_ylim())
ax_right.set_xticks([])
ax_right.yaxis.tick_right()

outName = "/t3home/gcelotto/BTV/plots/pTResolution2D.png"
fig.savefig(outName, bbox_inches='tight')
print("Saving in ", outName)
