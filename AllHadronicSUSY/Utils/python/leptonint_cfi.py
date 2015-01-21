import FWCore.ParameterSet.Config as cms

leptonint = cms.EDProducer('LeptonInt',
srcEle = cms.InputTag("selectedIDIsoElectrons"),
srcMuon = cms.InputTag("slimmedMuons"),
srcPV=cms.InputTag("offlineSlimmedPrimaryVertices"),
LeptonTag = cms.VInputTag(cms.InputTag('selectedIDIsoMuons'),cms.InputTag('selectedIDIsoElectrons')),
)
