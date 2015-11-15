import FWCore.ParameterSet.Config as cms

toyPFCandidates = cms.EDProducer("ToyPFCandidateProducer",
    src = cms.InputTag('genParticles')
)
