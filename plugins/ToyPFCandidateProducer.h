#ifndef RecoTauTag_GenTauAnalysisTools_ToyPFCandidateProducer_h
#define RecoTauTag_GenTauAnalysisTools_ToyPFCandidateProducer_h

/** \class ToyPFCandidateProducer
 *
 * Produce "toy MC" reco::PFCandidate objects 
 * corresponding to Monte Carlo "truth" information 
 * 
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

class ToyPFCandidateProducer : public edm::EDProducer
{
 public:
  // constructor 
  explicit ToyPFCandidateProducer(const edm::ParameterSet&);
  ~ToyPFCandidateProducer();
  
  void produce(edm::Event&, const edm::EventSetup&);
  
 private:
  edm::EDGetTokenT<reco::GenParticleCollection> tokenGenParticles_;
};

#endif
