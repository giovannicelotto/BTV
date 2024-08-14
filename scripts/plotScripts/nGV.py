import matplotlib.pyplot as plt
import urpoot
import glob
path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024Jul10Tuple"
fileNames = glob.glob(path + "/*.root")[:10]
for fileName in fileNames:
    f =urpoot.open(fileName)
    tree = f['tree']
    branches = tree.arrays()
    maxentries = tree.num_entries

    event = branches["event"]
    pdgID = branches["pdgID"]