import FWCore.ParameterSet.Config as cms

import copy
import os
import re

process = cms.Process("makeToyTaus")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff') 
process.load('Configuration.Geometry.GeometryExtended2017_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
##from Configuration.AlCa.autoCond import autoCond
##print "autoCond = %s" % autoCond
##process.GlobalTag.globaltag = cms.string( autoCond[ 'upgrade2017' ] )
process.GlobalTag.globaltag = cms.string('DES17_61_V5::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(                                
        '/store/user/veelken/Phys14/AODs/tth_HiggsToTauTau_AOD_gen_2b_4tauh.root'                                
    ),
)

process.dumpPATTauSequence = cms.Sequence()

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.dumpPATTauSequence += process.tauGenJets
process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
    'oneProng2Pi0', 
    'oneProngOther',
    'threeProng0Pi0', 
    'threeProng1Pi0', 
    'threeProngOther', 
    'rare'
)
process.tauGenJetsSelectorAllHadrons.filter = cms.bool(True)
process.dumpPATTauSequence += process.tauGenJetsSelectorAllHadrons

process.selectedTauGenJets = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    cut = cms.string("pt > 20. & abs(eta) < 2.3"),
    filter = cms.bool(True)
)
process.dumpPATTauSequence += process.selectedTauGenJets

process.dumpGenTaus = cms.EDAnalyzer("DumpGenTauSummary",
    src = cms.InputTag('genParticles')
)
process.dumpPATTauSequence += process.dumpGenTaus

#--------------------------------------------------------------------------------
# produce collection of "toy" PFCandidates 

process.toyPFCandidates = cms.EDProducer("ToyPFCandidateProducer",
    src = cms.InputTag('genParticles')
)
process.dumpPATTauSequence += process.toyPFCandidates
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run jet reconstruction on "toy" PFCandidates

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets

process.ak4ToyPFJets = ak4PFJets.clone(
   src = cms.InputTag('toyPFCandidates', 'pfCandidates')
)
process.dumpPATTauSequence += process.ak4ToyPFJets
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run tau reconstruction on "toy" PFCandidates

process.genTauMatchedPFJets = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak4ToyPFJets'),
    srcNotToBeFiltered = cms.VInputTag('selectedTauGenJets'),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.dumpPATTauSequence += process.genTauMatchedPFJets

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag 
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('particleFlow'), cms.InputTag('toyPFCandidates', 'pfCandidates'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak4PFJets'), cms.InputTag('genTauMatchedPFJets'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('generalTracks'), cms.InputTag('toyPFCandidates', 'tracks'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('offlinePrimaryVertices'), cms.InputTag('toyPFCandidates', 'vertices'))

process.dumpPATTauSequence += process.PFTau
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.cleanPatTaus.preselection = cms.string('')
process.cleanPatTaus.checkOverlaps = cms.PSet()
process.cleanPatTaus.finalCut = cms.string("")

process.dumpPATTauSequence += process.makePatTaus
process.dumpPATTauSequence += process.selectedPatTaus
process.dumpPATTauSequence += process.cleanPatTaus
#--------------------------------------------------------------------------------

process.dumpPFCandidates = cms.EDAnalyzer("DumpPFCandidatesByROI",
    srcROIs = cms.InputTag('genTauMatchedPFJets'),
    srcPFCandidates = cms.InputTag('toyPFCandidates', 'pfCandidates'),
    srcVertex = cms.InputTag('toyPFCandidates', 'vertices'),
    dRcone = cms.double(0.3),
    minPt = cms.double(-1.)
)
##process.dumpPATTauSequence += process.dumpPFCandidates

process.dumpPATTaus = cms.EDAnalyzer("DumpPATTauSummary",
    src = cms.InputTag('patTaus')
)
process.dumpPATTauSequence += process.dumpPATTaus

process.p = cms.Path(process.dumpPATTauSequence)

processDumpFile = open('makeToyTaus.dump', 'w')
print >> processDumpFile, process.dumpPython()



