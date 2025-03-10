import numpy as np
from helpers.utilsForScript import getPdgMask
secondMatching = False
def matchingEvent(distances, svDaughters_svIdx, svDaughters_pt, svDaughters_eta, GenPart_genPartIdxMother_, oneDaughter, svDaughters_genPartIdx, nGenPart_, GenPart_pdgId_, mesonsDaughters):
    distances_filled = distances.copy()
    matchingKey ={}
    while (np.any(distances_filled < 997)):
        recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
        recoIdx = int(recoIdx)
        genIdx = int(genIdx)

        #tracks from sv
        recoTracksMask = (svDaughters_svIdx==recoIdx) & (svDaughters_pt>0.8) & (abs(svDaughters_eta)<2.5)


        genTracksMask =  np.in1d(np.arange(nGenPart_), mesonsDaughters[genIdx])
        #genTracksMask = (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[oneDaughter[genIdx]])
        ##genTracksMask = (genTracksMask) |  (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]) 
        ##genTracksMask = (genTracksMask) |  (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]]) 
        #genTracksMask = (genTracksMask) | ( (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]) &  (~getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]])))
        #genTracksMask = (genTracksMask) | ( (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]]) & (~getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]])))

        commonTracks = np.sum(np.in1d(svDaughters_genPartIdx[recoTracksMask], np.arange(nGenPart_)[genTracksMask]))
        if commonTracks>=1:
            distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
            matchingKey[genIdx]=recoIdx
        else:
            distances_filled[recoIdx, genIdx]=998

        
    # second matching
    if secondMatching:
        distances_filled = distances.copy()
        while (np.any(distances_filled < 997)):
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)
            skip=False
            if (genIdx in matchingKey.keys()):
                #gen idx already matched so fill the whole column
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                skip=True
            if (recoIdx in matchingKey.values()):
                # same for the recoIdx and row
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
                skip=True
            if skip:
                continue
            if ((distances_filled[recoIdx, genIdx]<0.1)):
                matchingKey[genIdx]=recoIdx
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            else:
                distances_filled[recoIdx, genIdx]=998
    return matchingKey



def matchingEvent_old(distances, svDaughters_svIdx, svDaughters_pt, svDaughters_eta, GenPart_genPartIdxMother_, oneDaughter, svDaughters_genPartIdx, nGenPart_):
    distances_filled = distances.copy()
    matchingKey ={}
    while (np.any(distances_filled < 997)):
        recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
        recoIdx = int(recoIdx)
        genIdx = int(genIdx)

        #tracks from sv
        recoTracksMask = (svDaughters_svIdx==recoIdx) & (svDaughters_pt>0.8) & (abs(svDaughters_eta)<2.5)



        genTracksMask = (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[oneDaughter[genIdx]])
        genTracksMask = (genTracksMask) |  (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]) 
        genTracksMask = (genTracksMask) |  (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]]) 
        #genTracksMask = (genTracksMask) | ( (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]) &  (~getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]])))
        #genTracksMask = (genTracksMask) | ( (GenPart_genPartIdxMother_==GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]]) & (~getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[GenPart_genPartIdxMother_[oneDaughter[genIdx]]]])))

        commonTracks = np.sum(np.in1d(svDaughters_genPartIdx[recoTracksMask], np.arange(nGenPart_)[genTracksMask]))
        if commonTracks>=1:
            distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
            matchingKey[genIdx]=recoIdx
        else:
            distances_filled[recoIdx, genIdx]=998

        
    # second matching
    if secondMatching:
        distances_filled = distances.copy()
        while (np.any(distances_filled < 997)):
            recoIdx, genIdx = np.unravel_index(np.argmin(distances_filled, axis=None), distances.shape)
            recoIdx = int(recoIdx)
            genIdx = int(genIdx)
            skip=False
            if (genIdx in matchingKey.keys()):
                #gen idx already matched so fill the whole column
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                skip=True
            if (recoIdx in matchingKey.values()):
                # same for the recoIdx and row
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
                skip=True
            if skip:
                continue
            if ((distances_filled[recoIdx, genIdx]<0.1)):
                matchingKey[genIdx]=recoIdx
                distances_filled[:, genIdx]=[998]*distances_filled.shape[0]
                distances_filled[recoIdx, :]=[998]*distances_filled.shape[1]
            else:
                distances_filled[recoIdx, genIdx]=998
    return matchingKey