import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          l1results  = cms.InputTag("gtStage2Digis::RECO"),
                          #add interesting paths
                          paths = cms.vstring("HLT_Mu7_IP4",
                                              "HLT_Mu8_IP6",
                                              "HLT_Mu8_IP5",
                                              "HLT_Mu8_IP3",
                                              "HLT_Mu8p5_IP3p5",
                                              "HLT_Mu9_IP6",
                                              "HLT_Mu9_IP5",
                                              "HLT_Mu9_IP4",
                                              "HLT_Mu10p5_IP3p5",
                                              "HLT_Mu12_IP6",
                                          ),
                           #add interesting seeds
                           seeds = cms.vstring("L1_SingleMu7er1p5",
                                               "L1_SingleMu8er1p5",
                                               "L1_SingleMu9er1p5",
                                               "L1_SingleMu10er1p5",
                                               "L1_SingleMu12er1p5",
                                               "L1_SingleMu22",
                                           ),
                            
)

trgTables = cms.Sequence(trgTable)

###########
# Modifiers
###########

from PhysicsTools.BParkingNano.modifiers_cff import *

BToKEE_DiEle.toModify(trgTable,
                      paths = ['HLT_DoubleEle10_eta1p22_mMax6',
                               'HLT_DoubleEle9p5_eta1p22_mMax6',
                               'HLT_DoubleEle9_eta1p22_mMax6',
                               'HLT_DoubleEle8p5_eta1p22_mMax6',
                               'HLT_DoubleEle8_eta1p22_mMax6',
                               'HLT_DoubleEle7p5_eta1p22_mMax6',
                               'HLT_DoubleEle7_eta1p22_mMax6',
                               'HLT_DoubleEle6p5_eta1p22_mMax6',
                               'HLT_DoubleEle6_eta1p22_mMax6',
                               'HLT_DoubleEle5p5_eta1p22_mMax6',
                               'HLT_DoubleEle5_eta1p22_mMax6',
                               'HLT_DoubleEle4p5_eta1p22_mMax6',
                               'HLT_DoubleEle4_eta1p22_mMax6',],
                      seeds = ['L1_DoubleEG11_er1p2_dR_Max0p6',
                               'L1_DoubleEG10p5_er1p2_dR_Max0p6',
                               'L1_DoubleEG10_er1p2_dR_Max0p6',
                               'L1_DoubleEG9p5_er1p2_dR_Max0p6',
                               'L1_DoubleEG9_er1p2_dR_Max0p7',
                               'L1_DoubleEG8p5_er1p2_dR_Max0p7',
                               'L1_DoubleEG8_er1p2_dR_Max0p7',
                               'L1_DoubleEG7p5_er1p2_dR_Max0p7',
                               'L1_DoubleEG7_er1p2_dR_Max0p8',
                               'L1_DoubleEG6p5_er1p2_dR_Max0p8',
                               'L1_DoubleEG6_er1p2_dR_Max0p8',
                               'L1_DoubleEG5p5_er1p2_dR_Max0p8',
                               'L1_DoubleEG5_er1p2_dR_Max0p9',
                               'L1_DoubleEG4p5_er1p2_dR_Max0p9',
                               'L1_DoubleEG4_er1p2_dR_Max0p9',],
)

BToKMuMu_DiMuon.toModify(trgTable,
                         paths = ["HLT_DoubleMu4_JpsiTrk_Displaced"],
                         seeds = ["L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
                                  "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
                                  "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
                                  "L1_DoubleMu4_SQ_OS_dR_Max1p2",],
)
