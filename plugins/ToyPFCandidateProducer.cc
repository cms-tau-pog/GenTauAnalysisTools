#include "RecoTauTag/GenTauAnalysisTools/plugins/ToyPFCandidateProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include <TMath.h>

#include <vector>

ToyPFCandidateProducer::ToyPFCandidateProducer(const edm::ParameterSet& cfg)
{ 
  edm::InputTag srcGenParticles = cfg.getParameter<edm::InputTag>("src");
  tokenGenParticles_ = consumes<reco::GenParticleCollection>(srcGenParticles);
  produces<reco::PFCandidateCollection>("pfCandidates");
  produces<reco::VertexCollection>("vertices");
  produces<reco::TrackCollection>("tracks");
}

ToyPFCandidateProducer::~ToyPFCandidateProducer()
{
// nothing to be done yet...
}

std::vector<reco::PFCandidatePtr> getPFCandidatePtrs(reco::PFCandidateRefProd& pfCandidateRefs, const std::vector<unsigned>& indices)
{
  std::vector<reco::PFCandidatePtr> retVal;
  for ( std::vector<unsigned>::const_iterator idx = indices.begin();
	idx != indices.end(); ++idx ){
    retVal.push_back(edm::refToPtr(reco::PFCandidateRef(pfCandidateRefs, *idx)));
  }
  return retVal;
}

void ToyPFCandidateProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ToyPFCandidateProducer::produce>:" << std::endl;

  std::auto_ptr<reco::VertexCollection> vertices(new reco::VertexCollection());
  reco::VertexRefProd vertexRefs = evt.getRefBeforePut<reco::VertexCollection>("vertices");
  //std::cout << " vertexRefs: productId = " << vertexRefs.id().processIndex() << ":" << vertexRefs.id().productIndex() << std::endl;
  
  std::auto_ptr<reco::TrackCollection> tracks(new reco::TrackCollection());
  reco::TrackRefProd trackRefs = evt.getRefBeforePut<reco::TrackCollection>("tracks");
  //std::cout << " trackRefs: productId = " << trackRefs.id().processIndex() << ":" << trackRefs.id().productIndex() << std::endl;

  std::auto_ptr<reco::PFCandidateCollection> pfCandidates(new reco::PFCandidateCollection());

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(tokenGenParticles_, genParticles);

  size_t idxTrack = 0;

  size_t numGenParticles = genParticles->size();
  for ( size_t idxGenParticle = 0; idxGenParticle < numGenParticles; ++idxGenParticle ) {
    reco::GenParticleRef genParticle(genParticles, idxGenParticle);
    int genParticle_pdgId = TMath::Abs(genParticle->pdgId());
    if ( genParticle_pdgId == 12 || genParticle_pdgId == 14 || genParticle_pdgId == 16 ) continue;
    int genParticle_status = genParticle->status();
    if ( genParticle_status != 1 ) continue;
    if ( vertices->size() == 0 ) {
      reco::Vertex::Point vertex_point(genParticle->vertex().x(), genParticle->vertex().y(), genParticle->vertex().z());
      reco::Vertex::Error vertex_cov;
      double vertex_chi2 = 10.;
      double vertex_ndof = 20.;
      reco::Vertex vertex(vertex_point, vertex_cov, vertex_chi2, vertex_ndof, 0);
      vertices->push_back(vertex);
    }
    double genParticle_charge = genParticle->charge();
    reco::PFCandidate::ParticleType pfCandidate_type = reco::PFCandidate::X;
    int pfCandidate_charge = 0;
    bool isCharged = false;
    if ( TMath::Abs(genParticle_charge) > 0.5 && TMath::Abs(genParticle->eta()) < 2.5 ) {
      if      ( genParticle_pdgId ==  11 ) pfCandidate_type = reco::PFCandidate::e;
      else if ( genParticle_pdgId ==  13 ) pfCandidate_type = reco::PFCandidate::mu;
      else                                 pfCandidate_type = reco::PFCandidate::h;
      pfCandidate_charge = TMath::Nint(genParticle->charge());
      isCharged = true;
    } else if ( TMath::Abs(genParticle->eta()) < 3.0 ) {
      if      ( genParticle_pdgId ==  22 ) pfCandidate_type = reco::PFCandidate::gamma;
      else if ( genParticle_pdgId == 111 ) pfCandidate_type = reco::PFCandidate::gamma;
      else if ( genParticle_pdgId ==  11 ) pfCandidate_type = reco::PFCandidate::gamma;
      else                                 pfCandidate_type = reco::PFCandidate::h0;
    } else {
      if      ( genParticle_pdgId ==  22 ) pfCandidate_type = reco::PFCandidate::egamma_HF;
      else if ( genParticle_pdgId == 111 ) pfCandidate_type = reco::PFCandidate::egamma_HF;
      else if ( genParticle_pdgId ==  11 ) pfCandidate_type = reco::PFCandidate::egamma_HF;
      else                                 pfCandidate_type = reco::PFCandidate::h_HF;
    }
    reco::PFCandidate pfCandidate(pfCandidate_charge, genParticle->p4(), pfCandidate_type);
    double ecalEnergy = 0.;
    double hcalEnergy = 0.;
    if ( pfCandidate_type == reco::PFCandidate::mu ) {
      ecalEnergy = 0.5;
      if ( genParticle->energy() < ecalEnergy ) ecalEnergy = genParticle->energy();
      hcalEnergy = 2.5;
      if ( (genParticle->energy() - ecalEnergy) < hcalEnergy ) hcalEnergy = genParticle->energy() - ecalEnergy;
    } else if ( pfCandidate_type == reco::PFCandidate::e         ||
		pfCandidate_type == reco::PFCandidate::gamma     ||
		pfCandidate_type == reco::PFCandidate::egamma_HF ) {
      ecalEnergy = genParticle->energy();
      hcalEnergy = 0.;
    } else { 
      ecalEnergy = 0.2*genParticle->energy();
      hcalEnergy = 0.8*genParticle->energy();
    }
    pfCandidate.setEcalEnergy(ecalEnergy, ecalEnergy);
    pfCandidate.setHcalEnergy(hcalEnergy, hcalEnergy);
    if ( isCharged ) {
      assert(vertices->size() >= 1);
      reco::Vertex& vertex = vertices->front();
      double track_chi2 = 5.;
      double track_ndof = 10.;
      reco::TrackBase::Point track_vtx(vertex.position().x(), vertex.position().y(), vertex.position().z());
      reco::TrackBase::Vector track_momentum(genParticle->px(), genParticle->py(), genParticle->pz());
      int track_charge = TMath::Nint(genParticle->charge());
      reco::TrackBase::CovarianceMatrix track_cov;
      track_cov(0,0) = 0.1;
      track_cov(1,1) = 0.1;
      track_cov(2,2) = 0.1;
      track_cov(3,3) = 0.1;
      track_cov(4,4) = 0.1;
      reco::TrackBase::TrackAlgorithm track_algorithm = reco::TrackBase::ctf;
      reco::TrackBase::TrackQuality track_quality = reco::TrackBase::highPurity;
      reco::Track track(track_chi2, track_ndof, track_vtx, track_momentum, track_charge, track_cov, track_algorithm, track_quality);
      track.appendHitPattern(PXBDetId(1, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(PXBDetId(2, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(PXBDetId(3, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TIBDetId(1, 1, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TIBDetId(2, 1, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TIBDetId(3, 1, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TIBDetId(4, 1, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(1, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(2, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(3, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(4, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(5, 1, 0, 0, 0), TrackingRecHit::valid);
      track.appendHitPattern(TOBDetId(6, 1, 0, 0, 0), TrackingRecHit::valid);
      tracks->push_back(track);
      reco::TrackRef trackRef(trackRefs, idxTrack);
      pfCandidate.setTrackRef(trackRef);
      pfCandidate.setVertex(vertex.position());
      //reco::TrackBaseRef trackBaseRef(trackRefs, idxTrack);
      //vertex.add(trackBaseRef);
      ++idxTrack;
    }
    pfCandidates->push_back(pfCandidate);
  }

  evt.put(pfCandidates, "pfCandidates");
  evt.put(vertices, "vertices");
  evt.put(tracks, "tracks");
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ToyPFCandidateProducer);
