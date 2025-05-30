// Merges the PFPackedCandidates and Lost tracks
// beam spot readout in case dcasig to be calculated wrt beam spot
// currently computed wrt triggeringMuon vertex

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h" //GC to include Svcut

#include "helper.h"

class TrackMerger : public edm::global::EDProducer<> {


public:

  //would it be useful to give this a bit more standard structure?
  explicit TrackMerger(const edm::ParameterSet &cfg):
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    beamSpotSrc_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgLeptonToken_(consumes<edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("trgLepton"))),
    muonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    eleToken_(consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("pfElectrons"))),
    vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 
    svs_(consumes<edm::View<reco::VertexCompositePtrCandidate>>(cfg.getParameter<edm::InputTag>("svSrc"))), //GC to add linker tracks->SV
    svCut_(cfg.getParameter<std::string>("svCut"), true),                                                   // GC
    lowptele_(),
    lowpteleTag_(cfg.getParameter<edm::InputTag>("lowPtElectrons")),
    trkPtCut_(cfg.getParameter<double>("trkPtCut")), //0.5
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dzTrg_cleaning_(cfg.getParameter<double>("dzTrg_cleaning")),
    drTrg_Cleaning_(cfg.getParameter<double>("drTrg_Cleaning")),
    dcaSig_(cfg.getParameter<double>("dcaSig")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")),
    filterTrack_(cfg.getParameter<bool>("filterTrack")),
    dlenMin_(cfg.getParameter<double>("dlenMin")),      //GC
    dlenSigMin_(cfg.getParameter<double>("dlenSigMin")) //GC
{
  if ( !lowpteleTag_.label().empty() ) {
    lowptele_ = consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("lowPtElectrons"));
  }
    produces<pat::CompositeCandidateCollection>("SelectedTracks");  
    produces<TransientTrackCollection>("SelectedTransientTracks");  
}

  ~TrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<edm::View<reco::Candidate> > trgLeptonToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> eleToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svs_;  // GC needed to add linker tracks->SV
  const StringCutObjectSelector<reco::Candidate> svCut_;                // GC
  edm::EDGetTokenT<pat::ElectronCollection> lowptele_;
  const edm::InputTag lowpteleTag_;

  //selections                                                                 
  const double trkPtCut_;
  const double trkEtaCut_;
  const double dzTrg_cleaning_;
  const double drTrg_Cleaning_;
  const double dcaSig_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
  const bool filterTrack_;
  const double dlenMin_;    //GC
  const double dlenSigMin_; //GC
};






void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  //input
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from Event" ;
  }  
  const reco::BeamSpot& beamSpot = *beamSpotHandle;
  const auto& bField = stp.getData(bFieldToken_);
  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);
  edm::Handle<edm::View<reco::Candidate> > trgLeptons;
  evt.getByToken(trgLeptonToken_, trgLeptons);
  edm::Handle<edm::View<reco::VertexCompositePtrCandidate>> svsIn;  //GC
  evt.getByToken(svs_, svsIn);                                   //GC

  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);
  edm::Handle<pat::ElectronCollection> pfele;
  evt.getByToken(eleToken_, pfele);
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  edm::Handle<pat::ElectronCollection> lowptele;
  if ( !lowpteleTag_.label().empty() ) { evt.getByToken(lowptele_, lowptele); }

  //for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();

  //ok this was CompositeCandidateCollection 
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out      (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection>          trans_tracks_out(new TransientTrackCollection);

   std::vector< std::pair<pat::CompositeCandidate,reco::TransientTrack> > vectrk_ttrk; 
  //try topreserve same logic avoiding the copy of the full collection
  /*
  //correct logic but a bit convoluted -> changing to smthn simpler
   std::vector<pat::PackedCandidate> totalTracks(*tracks);
   totalTracks.insert(totalTracks.end(),lostTracks->begin(),lostTracks->end());
  */
 
  // for loop is better to be range based - especially for large ensembles  
  for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];

    //arranging cuts for speed
    if(!trk.hasTrackDetails()) continue;
    if(abs(trk.pdgId()) != 211) continue; //do we want also to keep muons?
    if(trk.pt() < trkPtCut_ ) continue;
    if(fabs(trk.eta()) > trkEtaCut_) continue;

    if( (trk.bestTrack()->normalizedChi2() < trkNormChiMin_ &&
         trkNormChiMin_>=0 ) ||
        (trk.bestTrack()->normalizedChi2() > trkNormChiMax_ &&
         trkNormChiMax_>0)  )    continue; 

    bool skipTrack=true;
    float dzTrg = 0.0;
    for (const auto & trg: *trgLeptons){
      //remove tracks inside trg leptons jet
      if(reco::deltaR(trk, trg) < drTrg_Cleaning_ && drTrg_Cleaning_ >0) 
        continue;
      //if dz is negative it is deactivated
      if((fabs(trk.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0))
        continue;
      skipTrack=false;
      dzTrg = trk.vz() - trg.vz();
      break; // at least for one trg lepton to pass this cuts
    }
    // if (trgLeptons->empty()) { skipTrack=false; } //@@ needed???

    // if track is closer to at least a triggering lepton keep it
    if (filterTrack_ && skipTrack ) continue;

    // high purity requirment applied only in packedCands
    if( iTrk < nTracks && !trk.trackHighPurity()) continue;
    const reco::TransientTrack trackTT( (*trk.bestTrack()) , &bField);
    //distance closest approach in x,y wrt beam spot
    std::pair<double,double> DCA = computeDCA(trackTT, beamSpot);
    float DCABS = DCA.first;
    float DCABSErr = DCA.second;
    float DCASig = (DCABSErr != 0 && float(DCABSErr) == DCABSErr) ? fabs(DCABS/DCABSErr) : -1;
    if (DCASig >  dcaSig_  && dcaSig_ >0) continue;

    // clean tracks wrt to all muons
    int matchedToMuon       = 0;
    int matchedToLooseMuon  = 0;
    int matchedToSoftMuon   = 0;
    int matchedToMediumMuon = 0;
    for (const pat::Muon &imutmp : *muons) {
        for (unsigned int i = 0; i < imutmp.numberOfSourceCandidatePtrs(); ++i) {
            if (! ((imutmp.sourceCandidatePtr(i)).isNonnull() && 
                   (imutmp.sourceCandidatePtr(i)).isAvailable())
               )   continue;
            
            const edm::Ptr<reco::Candidate> & source = imutmp.sourceCandidatePtr(i);
            if (source.id() == tracks.id() && source.key() == iTrk){
                matchedToMuon =1;
                if (imutmp.isLooseMuon())    matchedToLooseMuon  = 1;
                if (imutmp.isSoftMuon(PV))   matchedToSoftMuon   = 1;
                if (imutmp.isMediumMuon())   matchedToMediumMuon = 1;
                break;
            }
        }
    }
    int matchedToSV = -1;
    unsigned int svIdx = 0;
    VertexDistance3D vdist;
    for (const auto& sv : *svsIn) {
      if (svCut_(sv)) {
        Measurement1D dl =  vdist.distance(PV, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
        if (sv.numberOfDaughters()==0){
          std::cout<<"No daughters here";
        }
        if (dl.value() > dlenMin_ and dl.significance() > dlenSigMin_) {

          // cut on SV quality
          //std::cout<<"Number of daughters : "<<sv.numberOfDaughters()<<std::endl;
          for (unsigned int id = 0; id < sv.numberOfDaughters(); ++id) {
            // reco::Candidate* daughter is a PF candidate
            // trk is also a PF candidate:
            //     const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];
            // try to change const reco::Track (*)& trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];
            // lock at the continue statements
            const reco::Candidate* daughter = sv.daughter(id);
            //std::cout<<svIdx<<"  "<<daughter->pt()<<" "<<daughter->eta()<<std::endl;
            if (reco::deltaR(trk, *daughter)<0.00001){
              matchedToSV = svIdx;
              //std::cout<<"Track phi "<<trk.phi()<<std::endl;
              //std::cout<<"SV Daughter phi "<<daughter->phi()<<std::endl;
              //std::cout<<"Track Eta "<<trk.eta()<<std::endl;
              //std::cout<<"SV Daughter Eta "<<daughter->eta()<<std::endl;
              //std::cout<<"Track PT "<<trk.pt()<<std::endl;
              //std::cout<<"SV Daughter PT "<<daughter->pt()<<std::endl<<std::endl;
              break;
            }
          }
        ++svIdx;
        }
      }
    }

    // clean tracks wrt to all pf electrons
    int matchedToEle        = 0;
    for (const pat::Electron &ietmp : *pfele) {
        for (unsigned int i = 0; i < ietmp.numberOfSourceCandidatePtrs(); ++i) {
            
            if (! ((ietmp.sourceCandidatePtr(i)).isNonnull() && 
                   (ietmp.sourceCandidatePtr(i)).isAvailable())
               )   continue;
            const edm::Ptr<reco::Candidate> & source = ietmp.sourceCandidatePtr(i);
            if (source.id() == tracks.id() && source.key() == iTrk){
                matchedToEle =1;
                break;
            }        
        }

    }

    // clean tracks wrt to all low pT electrons
    edm::Ptr<pat::PackedCandidate> track;
    if ( iTrk < nTracks ) { track = edm::Ptr<pat::PackedCandidate>(tracks,iTrk); }
    else { track = edm::Ptr<pat::PackedCandidate>(lostTracks,iTrk-nTracks); }
    int matchedToLowPtEle = 0;
    if ( !lowpteleTag_.label().empty() ) {
      for ( auto const& ele : *lowptele ) {
	reco::GsfTrackRef gsf = ele.gsfTrack();
	if ( iTrk < nTracks ) {
	  const edm::Ptr<pat::PackedCandidate>* packed = ele.hasUserData("ele2packed") ?
	    ele.userData<edm::Ptr<pat::PackedCandidate> >("ele2packed") : NULL;
	  if ( packed != NULL && track == *packed ) { matchedToLowPtEle = 1; break; }
	} else {
	  const edm::Ptr<pat::PackedCandidate>* lost = ele.hasUserData("ele2lost") ?
	    ele.userData<edm::Ptr<pat::PackedCandidate> >("ele2lost") : NULL;
	  if ( lost != NULL && track == *lost ) { matchedToLowPtEle = 1; break; }
	}
      }
    }

    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserInt("isPacked", (iTrk < nTracks));
    pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);      
    pcand.addUserFloat("dxy", trk.dxy());
    pcand.addUserFloat("dxyS", trk.dxy()/trk.dxyError());
    pcand.addUserFloat("dz", trk.dz()); 
    pcand.addUserFloat("dzS", trk.dz()/trk.dzError());
    pcand.addUserFloat("DCASig", DCASig);
    pcand.addUserFloat("dzTrg", dzTrg);
    pcand.addUserInt("isMatchedToMuon", matchedToMuon);
    pcand.addUserInt("isMatchedToLooseMuon", matchedToLooseMuon);
    pcand.addUserInt("isMatchedToSoftMuon", matchedToSoftMuon);
    pcand.addUserInt("isMatchedToMediumMuon", matchedToMediumMuon);
    pcand.addUserInt("isMatchedToEle", matchedToEle);
    pcand.addUserInt("isMatchedToLowPtEle", matchedToLowPtEle);
    pcand.addUserInt("nValidHits", trk.bestTrack()->found());
    pcand.addUserInt("keyPacked", iTrk);
    pcand.addUserInt("skipTrack",skipTrack);
    pcand.addUserInt("matchedToSV", matchedToSV);
    //adding the candidate in the composite stuff for fit (need to test)
    if ( iTrk < nTracks )
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( tracks, iTrk ));
    else 
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( lostTracks, iTrk-nTracks ));   
 
  //in order to avoid revoking the sxpensive ttrack builder many times and still have everything sorted, we add them to vector of pairs
   vectrk_ttrk.emplace_back( std::make_pair(pcand,trackTT ) );   
  }

  // sort to be uniform with leptons
  std::sort( vectrk_ttrk.begin(), vectrk_ttrk.end(), 
             [] ( auto & trk1, auto & trk2) -> 
                  bool {return (trk1.first).pt() > (trk2.first).pt();} 
           );

  // finnaly save ttrks and trks to the correct _out vectors
  for ( auto & trk: vectrk_ttrk){
    tracks_out -> emplace_back( trk.first);
    trans_tracks_out -> emplace_back(trk.second);
  }

  evt.put(std::move(tracks_out),       "SelectedTracks");
  evt.put(std::move(trans_tracks_out), "SelectedTransientTracks");
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackMerger);
