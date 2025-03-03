import numpy as np
import sys
sys.path.append("/work/gcelotto/BTV/scripts/tuplizer/helpers")
from utilsForScript import getPdgMask
import awkward as ak
def getMesons(params, GenPart_pdgId_, GenPart_eta_, GenPart_pt_, nGenPart_, GenPart_status_, GenPart_charge, GenPart_genPartIdxMother_):
    """
    Identifies mesons from generated particle data based on PDG ID, kinematic constraints, and stability criteria.

    Parameters:
        params (dict): Dictionary containing selection criteria such as max_abs_eta and min_pt_mother.
        GenPart_pdgId_ (array-like): Array of PDG IDs of GenParticles.
        GenPart_eta_ (array-like): Array of pseudorapidities (η) of GenParticles.
        GenPart_pt_ (array-like): Array of transverse momenta (pT) of GenParticles.
        nGenPart_ (int): Total number of GenParticles.
        GenPart_status_ (array-like): Status flags indicating if a particle is stable.
        GenPart_charge (array-like): Charges of GenParticles.
        GenPart_genPartIdxMother_ (array-like): Indices of mother particles for each GenParticles.

    Returns:
        numpy.ndarray: Array of indices corresponding to baryons that pass the selection criteria.
        Baryons with
            * pT > 10
            * |eta| <2.5
            * pdgID in list
            * with >=2 charged stable daughters with pT> 0.8, |eta| < 2.5
    """

    pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId_)   # list of pdgId consider as genVertices
    etaMask = abs(GenPart_eta_)<params['max_abs_eta']  # limit the genvertices to GenHadrons within tracker acceptance
    ptMask = GenPart_pt_>params['min_pt_mother']          # require pT of GenHadrons larger than 10 GeV
    mesons = np.arange(nGenPart_)[(pdgMask) & (etaMask) & (ptMask) ] #index of the mesons in nGenPart
    mesonsDaughters =[]
    nStableDaughters = np.zeros(len(mesons))
    for idx, mesIdx in enumerate(mesons):
        currentMesonDaughter = []
        for genPartIdx in np.arange(nGenPart_)[(GenPart_status_==1) & (abs(GenPart_charge)==1) & (GenPart_pt_>0.8) & (abs(GenPart_eta_)<2.5)]:
            #input("next dau %d"%(GenPart_pdgId_[genPartIdx]))
            gp = genPartIdx
            #if (GenPart_genPartIdxMother_[gp]==mesIdx):
            #    nStableDaughters[idx]=nStableDaughters[idx]+1
    # Definition of 01/03
            while (GenPart_genPartIdxMother_[gp]!=-1):
                if (GenPart_genPartIdxMother_[gp]==mesIdx):
                    nStableDaughters[idx]=nStableDaughters[idx]+1
                    currentMesonDaughter.append(genPartIdx)
                    #break the while and look for the next genPart
                    break
            # if the mother has displacement dont go back! It will be another vertex
            # Check also if the pdgID is different
                if (getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[gp]])) & (GenPart_pdgId_[GenPart_genPartIdxMother_[gp]] != GenPart_pdgId_[gp]):
                    break
                else:
                    gp = GenPart_genPartIdxMother_[gp]
        mesonsDaughters.append(currentMesonDaughter)
    # Limit the found baryons to those with at least 2 charged stable daughters pt>0.8 adn |eta|<2.5
    mesons = mesons[nStableDaughters>1]
    mesonsDaughters = ak.from_iter(mesonsDaughters)

    if len(mesonsDaughters)>0:
        mesonsDaughters = mesonsDaughters[nStableDaughters>1, :]
        
    return mesons, mesonsDaughters


def getMesons_old(params, GenPart_pdgId_, GenPart_eta_, GenPart_pt_, nGenPart_, GenPart_status_, GenPart_charge, GenPart_genPartIdxMother_):
    """
    Identifies mesons from generated particle data based on PDG ID, kinematic constraints, and stability criteria.

    Parameters:
        params (dict): Dictionary containing selection criteria such as max_abs_eta and min_pt_mother.
        GenPart_pdgId_ (array-like): Array of PDG IDs of GenParticles.
        GenPart_eta_ (array-like): Array of pseudorapidities (η) of GenParticles.
        GenPart_pt_ (array-like): Array of transverse momenta (pT) of GenParticles.
        nGenPart_ (int): Total number of GenParticles.
        GenPart_status_ (array-like): Status flags indicating if a particle is stable.
        GenPart_charge (array-like): Charges of GenParticles.
        GenPart_genPartIdxMother_ (array-like): Indices of mother particles for each GenParticles.

    Returns:
        numpy.ndarray: Array of indices corresponding to baryons that pass the selection criteria.
        Baryons with
            * pT > 10
            * |eta| <2.5
            * pdgID in list
            * with >=2 charged stable daughters with pT> 0.8, |eta| < 2.5
    """

    pdgMask = getPdgMask(GenPart_pdgId=GenPart_pdgId_)   # list of pdgId consider as genVertices
    etaMask = abs(GenPart_eta_)<params['max_abs_eta']  # limit the genvertices to GenHadrons within tracker acceptance
    ptMask = GenPart_pt_>params['min_pt_mother']          # require pT of GenHadrons larger than 10 GeV
    mesons = np.arange(nGenPart_)[(pdgMask) & (etaMask) & (ptMask) ] #index of the mesons in nGenPart
    #mesonsDaughters =[]
    nStableDaughters = np.zeros(len(mesons))
    for idx, mesIdx in enumerate(mesons):
        #currentMesonDaughter = []
        for genPartIdx in np.arange(nGenPart_)[(GenPart_status_==1) & (abs(GenPart_charge)==1) & (GenPart_pt_>0.8) & (abs(GenPart_eta_)<2.5)]:
            #input("next dau %d"%(GenPart_pdgId_[genPartIdx]))
            gp = genPartIdx
            #if (GenPart_genPartIdxMother_[gp]==mesIdx):
            #    nStableDaughters[idx]=nStableDaughters[idx]+1
    # Definition of 01/03
            while (GenPart_genPartIdxMother_[gp]!=-1):
                if (GenPart_genPartIdxMother_[gp]==mesIdx):
                    nStableDaughters[idx]=nStableDaughters[idx]+1
                    #currentMesonDaughter.append(genPartIdx)
                    #break the while and look for the next genPart
                    break
            # if the mother has displacement dont go back! It will be another vertex
            # Check also if the pdgID is different
                #if (getPdgMask(GenPart_pdgId_[GenPart_genPartIdxMother_[gp]])) & (GenPart_pdgId_[GenPart_genPartIdxMother_[gp]] != GenPart_pdgId_[gp]):
                #    break
                #else:
                gp = GenPart_genPartIdxMother_[gp]
    # Limit the found baryons to those with at least 2 charged stable daughters pt>0.8 adn |eta|<2.5
    mesons = mesons[nStableDaughters>1]     
    return mesons