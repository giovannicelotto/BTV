import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

tracksBPark = cms.EDProducer('TrackMerger',
                             beamSpot   = cms.InputTag("offlineBeamSpot"),
                             trgLepton    = cms.InputTag("slimmedMuons"),
                             tracks     = cms.InputTag("packedPFCandidates"),
                             #tracks     = cms.InputTag("generalTracks"),
                             lostTracks = cms.InputTag("lostTracks"),
                             trkPtCut = cms.double(0.5),
                             muons      = cms.InputTag("slimmedMuons"),
                             pfElectrons= cms.InputTag("slimmedElectrons"),
                             vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             lowPtElectrons=cms.InputTag("slimmedLowPtElectrons"),
                             trkEtaCut = cms.double(2.5),
                             filterTrack = cms.bool(False),
                             dzTrg_cleaning = cms.double(-1), #1.
                             drTrg_Cleaning = cms.double(-1), #0.03
                             dcaSig = cms.double(-1), #-100000
                             trkNormChiMin = cms.int32(-1),
                             trkNormChiMax = cms.int32(-1),
                             svSrc = cms.InputTag("slimmedSecondaryVertices"),      #GC
                             svCut = cms.string(""),                                #GC
                             dlenMin = cms.double(0),
                             dlenSigMin = cms.double(3),
                            )


trackBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksBPark:SelectedTracks"),
    cut = cms.string(""),
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection probe side for BPark after basic selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False), 
    variables = cms.PSet(
         CandVars,
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        isLostTrk = Var("userInt('isLostTrk')",int,doc="track from lostTrack collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        dzTrg = Var("userFloat('dzTrg')", float,doc="dz from the corresponding trigger lepton, in cm", precision=10),
        isMatchedToMuon = Var("userInt('isMatchedToMuon')",bool,doc="track was used to build a muon", precision=10),
        isMatchedToLooseMuon = Var("userInt('isMatchedToLooseMuon')",bool,doc="track was used to build a muon passing LooseID", precision=10),
        isMatchedToSoftMuon = Var("userInt('isMatchedToSoftMuon')",bool,doc="track was used to build a muon passing softID", precision=10),
        isMatchedToMediumMuon = Var("userInt('isMatchedToMediumMuon')",bool,doc="track was used to build a muon passing mediumID", precision=10),
        isMatchedToEle = Var("userInt('isMatchedToEle')",bool,doc="track was used to build a PF ele", precision=10),
        isMatchedToLowPtEle = Var("userInt('isMatchedToLowPtEle')",bool,doc="track was used to build a low-pT ele", precision=10),
        nValidHits = Var("userInt('nValidHits')", int,doc="Number of valid hits on track", precision=10),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        skipTrack = Var("userInt('skipTrack')",bool,doc="Is track skipped (due to small dR or large dZ w.r.t. trigger)?"),
        matchedToSV = Var("userInt('matchedToSV')", int, doc="Index of SV the track is matched to", precision=10)
        ),
)


tracksBParkMCMatchForTable = cms.EDProducer("MCMatcher",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = trackBParkTable.src,                     # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),  # final mc-truth particle collection
    mcPdgId     = cms.vint32(321,211),                     # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),             # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
)

tracksBParkMCMatchEmbedded = cms.EDProducer(
    'CompositeCandidateMatchEmbedder',
    src = trackBParkTable.src,
    matching = cms.InputTag("tracksBParkMCMatchForTable")
)

tracksBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = tracksBParkMCMatchForTable.src,
    mcMap   = cms.InputTag("tracksBParkMCMatchForTable"),
    objName = trackBParkTable.name,
    objType = trackBParkTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 kaons or pions"),
)


tracksBParkSequence = cms.Sequence(tracksBPark)
tracksBParkTables = cms.Sequence(trackBParkTable)
tracksBParkMC = cms.Sequence(tracksBParkSequence + tracksBParkMCMatchForTable + tracksBParkMCMatchEmbedded + tracksBParkMCTable)

###########
# Modifiers
###########

from PhysicsTools.BParkingNano.modifiers_cff import *

_modifiers = BToKMuMu_OpenConfig | BToKEE_OpenConfig
_modifiers.toModify(tracksBPark,
                    trkPtCut=0.5,
                    trkEtaCut=2.5,
                    trkNormChiMin=-1,
                    trkNormChiMax=-1,
                    dcaSig=-100000,
                    #dzTrg_cleaning=-1.,
                    #drTrg_Cleaning=-1.,
                    filterTrack=False)

BToKEE_DiEle.toModify(tracksBPark,
                      trgLepton = "electronTrgSelector:trgElectrons",
                      lowPtElectrons = "") # don't use "slimmedLowPtElectrons"
