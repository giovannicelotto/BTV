// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      VertexTableProducer
//
/**\class VertexTableProducer VertexTableProducer.cc PhysicsTools/VertexTableProducer/plugins/VertexTableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Rizzi
//         Created:  Mon, 28 Aug 2017 09:26:39 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //GC
//#include "FWCore/Framework/interface/global/EDProducer.h" //GC
//
// class declaration
//
//GC new for pat format
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

class VertexTableProducer : public edm::stream::EDProducer<> {
public:
  explicit VertexTableProducer(const edm::ParameterSet&);
  ~VertexTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
  const edm::EDGetTokenT<edm::ValueMap<float>> pvsScore_;
  const edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svs_;
  //const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >    Cands_; //GC
  const StringCutObjectSelector<reco::Candidate> svCut_;
  const StringCutObjectSelector<reco::Vertex> goodPvCut_;
  const std::string goodPvCutString_;
  const std::string pvName_;
  const std::string svName_;
  const std::string svDaughtersName_;
  const std::string svDoc_;
  const double dlenMin_, dlenSigMin_;
};

//
// constructors and destructor
//
VertexTableProducer::VertexTableProducer(const edm::ParameterSet& params):
      pvs_(consumes<std::vector<reco::Vertex>>(params.getParameter<edm::InputTag>("pvSrc"))),
      pvsScore_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("pvSrc"))),
      svs_(consumes<edm::View<reco::VertexCompositePtrCandidate>>(params.getParameter<edm::InputTag>("svSrc"))),
      //Cands_(consumes< std::vector<pat::PackedCandidate> >(params.getParameter<edm::InputTag>("packedCandidates"))), //GC
      svCut_(params.getParameter<std::string>("svCut"), true),
      goodPvCut_(params.getParameter<std::string>("goodPvCut"), true),
      goodPvCutString_(params.getParameter<std::string>("goodPvCut")),
      pvName_(params.getParameter<std::string>("pvName")),
      svName_(params.getParameter<std::string>("svName")),
      svDaughtersName_(params.getParameter<std::string>("svDaughtersName")),
      svDoc_(params.getParameter<std::string>("svDoc")),
      dlenMin_(params.getParameter<double>("dlenMin")),
      dlenSigMin_(params.getParameter<double>("dlenSigMin"))

{
  produces<nanoaod::FlatTable>("pv");
  produces<nanoaod::FlatTable>("otherPVs");
  produces<nanoaod::FlatTable>("svs");
  produces<nanoaod::FlatTable>("svDaughters");
  produces<edm::PtrVector<reco::Candidate>>("selCandSv");
  //produces<pat::CompositeCandidateCollection>("selCandSvDaughters");
  produces<std::vector<edm::Ptr<reco::Candidate>>>("selCandSvDaughters");

}

VertexTableProducer::~VertexTableProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void VertexTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle<edm::ValueMap<float>> pvsScoreIn;
  edm::Handle<std::vector<reco::Vertex>> pvsIn;
  iEvent.getByToken(pvs_, pvsIn);
  iEvent.getByToken(pvsScore_, pvsScoreIn);
  //edm::Handle<std::vector<pat::PackedCandidate> > cands;    //GC
	//iEvent.getByToken(Cands_, cands);                         //GC
  auto pvTable = std::make_unique<nanoaod::FlatTable>(1, pvName_, true);
  pvTable->addColumnValue<float>("ndof", (*pvsIn)[0].ndof(), "main primary vertex number of degree of freedom", 8);
  pvTable->addColumnValue<float>("x", (*pvsIn)[0].position().x(), "main primary vertex position x coordinate", 10);
  pvTable->addColumnValue<float>("y", (*pvsIn)[0].position().y(), "main primary vertex position y coordinate", 10);
  pvTable->addColumnValue<float>("z", (*pvsIn)[0].position().z(), "main primary vertex position z coordinate", 16);
  pvTable->addColumnValue<float>("chi2", (*pvsIn)[0].normalizedChi2(), "main primary vertex reduced chi2", 8);
  int goodPVs = 0;
  for (const auto& pv : *pvsIn)
    if (goodPvCut_(pv))
      goodPVs++;
  pvTable->addColumnValue<int>("npvs", pvsIn->size(), "total number of reconstructed primary vertices");
  pvTable->addColumnValue<int>(
      "npvsGood", goodPVs, "number of good reconstructed primary vertices. selection:" + goodPvCutString_);
  pvTable->addColumnValue<float>(
      "score", pvsScoreIn->get(pvsIn.id(), 0), "main primary vertex score, i.e. sum pt2 of clustered objects", 8);

  auto otherPVsTable =
      std::make_unique<nanoaod::FlatTable>((*pvsIn).size() > 4 ? 3 : (*pvsIn).size() - 1, "Other" + pvName_, false);
  std::vector<float> pvsz;
  std::vector<float> pvscores;
  for (size_t i = 1; i < (*pvsIn).size() && i < 4; i++) {
    pvsz.push_back((*pvsIn)[i].position().z());
    pvscores.push_back((*pvsScoreIn).get(pvsIn.id(), i));
  }
  otherPVsTable->addColumn<float>("z", pvsz, "Z position of other primary vertices, excluding the main PV", 8);
  otherPVsTable->addColumn<float>("score", pvscores, "scores of other primary vertices, excluding the main PV", 8);

  edm::Handle<edm::View<reco::VertexCompositePtrCandidate>> svsIn;
  iEvent.getByToken(svs_, svsIn);
  auto selCandSv = std::make_unique<PtrVector<reco::Candidate>>();
  auto selCandSvDaughters = std::make_unique<std::vector<edm::Ptr<reco::Candidate>>>();

  std::vector<float> dlen, dlenSig, pAngle, dxy, dxySig;
  std::vector<int> charge;
  std::vector<float> daughters_idx;
  //std::vector<float> daughters_pt;
  //std::vector<float> daughters_eta;
  //std::vector<float> daughters_phi;
  VertexDistance3D vdist;
  VertexDistanceXY vdistXY;

  size_t i = 0;
  size_t i_pass = 0; // new variable to count the index of the SVs that pass the cuts GC
  const auto& PV0 = pvsIn->front();
  for (const auto& sv : *svsIn) {
    if (svCut_(sv)) {
      Measurement1D dl =
          vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
      if (dl.value() > dlenMin_ and dl.significance() > dlenSigMin_) {
        dlen.push_back(dl.value());
        dlenSig.push_back(dl.significance());
        edm::Ptr<reco::Candidate> c = svsIn->ptrAt(i);
        selCandSv->push_back(c);
        double dx = (PV0.x() - sv.vx()), dy = (PV0.y() - sv.vy()), dz = (PV0.z() - sv.vz());
        double pdotv = (dx * sv.px() + dy * sv.py() + dz * sv.pz()) / sv.p() / sqrt(dx * dx + dy * dy + dz * dz);
        pAngle.push_back(std::acos(pdotv));
        Measurement1D d2d = vdistXY.distance(
            PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
        dxy.push_back(d2d.value());
        dxySig.push_back(d2d.significance());

        int sum_charge = 0;
        for (unsigned int id = 0; id < sv.numberOfDaughters(); ++id) {
          edm::Ptr<reco::Candidate> daughterPtr = sv.daughterPtr(id);
          selCandSvDaughters->push_back(daughterPtr);
          sum_charge += daughterPtr->charge();
        }
        charge.push_back(sum_charge);

        for (unsigned int id = 0; id < sv.numberOfDaughters(); ++id) {
          daughters_idx.push_back(i_pass); // use the index of the SV (passing the cuts) to know which daughter comes from which SV
          //daughters_eta.push_back(sv.daughter(id)->eta());
          //daughters_phi.push_back(sv.daughter(id)->phi());
          //daughters_pt.push_back(sv.daughter(id)->pt());
        }
        // increment the SV that pass the cut
        i_pass++;
      }
    }
    i++;
  }

  //edm::Handle<std::vector<pat::PackedCandidate> > cands;
	//iEvent.getByToken(Cands_, cands);
  //std::auto_ptr< std::vector<reco::Vertex> > outSv( new std::vector<reco::Vertex> );
  //for(size_t i=0;i< svsIn->size(); i++) {
	//	const reco::VertexCompositePtrCandidate &sv = (*svsIn)[i];	
	//	outSv->push_back(reco::Vertex(sv.vertex(),sv.vertexCovariance(),sv.vertexChi2(),sv.vertexNdof(),0));
  //  std::cout<<"\n\n"<<i<<std::endl;
	//	for(size_t j=0;j<sv.numberOfDaughters();j++){
	//                reco::TrackRef r;
  //                std::cout<<sv.daughterPtr(j).key()<<"  "<<sv.daughterPtr(j)->pt()<<"  "<<sv.daughterPtr(j)->eta()<<"  "<<sv.daughterPtr(j)->phi()<<std::endl;
                  //std::cout<<sv.daughterPtr(j)->pt()<<std::endl;
			//if(sv.daughterPtr(j).id() == cands.id()) {cd
	    //            	 r= TrackRef(oh,trackKeys[sv.daughterPtr(j).key()]); // use trackKeys because cand->track has gaps from neutral
			//} else {
//				std::cout << "vertex " << i << " using lost Track " << sv.daughterPtr(j).key()  << "  " << offsetAdd+sv.daughterPtr(j).key() << std::endl;  
      //                          r=TrackRef(oh,offsetAdd+sv.daughterPtr(j).key());  // use directly the key because addTracks is only charged
			//}
//        	        TrackBaseRef rr(r);
//			outSv->back().add(rr);
//
//		}	
//	}   
//
  auto svsTable = std::make_unique<nanoaod::FlatTable>(selCandSv->size(), svName_, false);
  // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  svsTable->addColumn<float>("dlen", dlen, "decay length in cm", 10);
  svsTable->addColumn<float>("dlenSig", dlenSig, "decay length significance", 10);
  svsTable->addColumn<float>("dxy", dxy, "2D decay length in cm", 10);
  svsTable->addColumn<float>("dxySig", dxySig, "2D decay length significance", 10);
  svsTable->addColumn<float>("pAngle", pAngle, "pointing angle, i.e. acos(p_SV * (SV - PV)) ", 10);
  svsTable->addColumn<int>("charge", charge, "sum of the charge of the SV tracks", 10);
  
  
  // Create a new Table
  auto svDaughtersTable = std::make_unique<nanoaod::FlatTable>(selCandSvDaughters->size(), svDaughtersName_, false);

  // Fill the table with per-daughter information
  //svDaughtersTable->addColumn<float>("eta", daughters_eta, "eta of the daughters of the SV", 10);
  //svDaughtersTable->addColumn<float>("phi", daughters_phi, "phi of the daughters of the SV", 10);
  //svDaughtersTable->addColumn<float>("pt", daughters_pt, "pt of the daughters of the SV", 10);
  svDaughtersTable->addColumn<int>("svIdx", daughters_idx, "index of the SV the daughter comes from", 10);  // This maps daughters to SVs

  // Put the new table into the event
  iEvent.put(std::move(svDaughtersTable), "svDaughters");

  iEvent.put(std::move(pvTable), "pv");
  iEvent.put(std::move(otherPVsTable), "otherPVs");
  iEvent.put(std::move(svsTable), "svs");
  iEvent.put(std::move(selCandSv), "selCandSv");
  iEvent.put(std::move(selCandSvDaughters), "selCandSvDaughters");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void VertexTableProducer::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void VertexTableProducer::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void VertexTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexTableProducer);
