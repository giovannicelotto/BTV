import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *



##################### User floats producers, selectors ##########################


##################### Tables for final output and docs ##########################
vertexTable = cms.EDProducer("VertexTableProducer",
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), 
    svSrc = cms.InputTag("slimmedSecondaryVertices"),  #slimmedSecondaryVertices GC
    svCut = cms.string(""),
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(3),
    pvName = cms.string("PV"),
    svName = cms.string("SV"),
    svDaughtersName = cms.string("svDaughters"),
    svDoc  = cms.string("secondary vertices from IVF algorithm"),
)

svCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("vertexTable:selCandSv"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("SV"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("position().x()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("position().y()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("position().z()", float, doc = "secondary vertex Z position, in cm",precision=14),
        ndof    = Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
        chi2    = Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
        ntracks = Var("numberOfDaughters()", "uint8", doc = "number of tracks"),
    ),
)

svDaughtersCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("vertexTable:selCandSvDaughters"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("svDaughters"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(
        P4Vars,
    ),
)


svDaughtersMCMatchForTable = cms.EDProducer("MCMatcher",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = svDaughtersCandidateTable.src,                     # final reco collection
    matched     = cms.InputTag("finalGenParticles"),  # final mc-truth particle collection
    mcPdgId     = cms.vint32(3334, 3322, 3312, 3222, 3122, 3112, 2212, 2112, 321, 211, 13, 11),                     # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),             # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
)

svDaughtersMCMatchEmbedded = cms.EDProducer(
    'CompositeCandidateMatchEmbedder',
    src = svDaughtersCandidateTable.src,
    matching = cms.InputTag("svDaughtersMCMatchForTable")
)

svDaughtersMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = svDaughtersMCMatchForTable.src,
    mcMap   = cms.InputTag("svDaughtersMCMatchForTable"),
    objName = svDaughtersCandidateTable.name,
    objType = svDaughtersCandidateTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 kaons or pions"),
)


svCandidateTable.variables.pt.precision=10
svCandidateTable.variables.phi.precision=12


#before cross linking
vertexTask = cms.Task()
#after cross linkining
vertexTablesTask = cms.Task( vertexTable, svCandidateTable,svDaughtersCandidateTable,svDaughtersMCMatchForTable, svDaughtersMCMatchEmbedded, svDaughtersMCTable  )