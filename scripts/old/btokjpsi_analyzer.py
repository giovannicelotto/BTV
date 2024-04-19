import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot


path = "/t3home/gcelotto/BTV/CMSSW_12_4_8/src/PhysicsTools/BParkingNano/test/BuToJpsiK_Run2_mc_124X.root"
f = uproot.open(path)
tree = f['Events']
branches = tree.arrays()
maxEntries = 100 #tree.num_entries
for ev in range(maxEntries):
    print("Ev %d"%ev)
    GenPart_pdgId = branches['GenPart_pdgId'][ev]
    if 13 not in abs(GenPart_pdgId):
        assert False

    print(GenPart_pdgId)



