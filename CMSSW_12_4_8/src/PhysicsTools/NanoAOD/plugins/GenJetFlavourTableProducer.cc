// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"               // GC
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"       // GC              
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

class GenJetFlavourTableProducer : public edm::stream::EDProducer<> {
    public:
        explicit GenJetFlavourTableProducer(const edm::ParameterSet &iConfig) :
            name_(iConfig.getParameter<std::string>("name")),
            src_(consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("src"))),
            cut_(iConfig.getParameter<std::string>("cut"), true),
            deltaR_(iConfig.getParameter<double>("deltaR")),
            jetFlavourInfosToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetFlavourInfos")))
        {
            produces<nanoaod::FlatTable>();
        }

        ~GenJetFlavourTableProducer() override {};

        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
            edm::ParameterSetDescription desc;
            desc.add<edm::InputTag>("src")->setComment("input genJet collection");
            desc.add<edm::InputTag>("jetFlavourInfos")->setComment("input flavour info collection");
            desc.add<std::string>("name")->setComment("name of the genJet FlatTable we are extending with flavour information");
            desc.add<std::string>("cut")->setComment("cut on input genJet collection");
            desc.add<double>("deltaR")->setComment("deltaR to match genjets");
            descriptions.add("genJetFlavourTable", desc);
        }

    private:
        void produce(edm::Event&, edm::EventSetup const&) override ;

        std::string name_;
        edm::EDGetTokenT<std::vector<reco::GenJet> > src_;
        const StringCutObjectSelector<reco::GenJet> cut_;
        const double deltaR_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

};

// ------------ method called to produce the data  ------------
void
GenJetFlavourTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenJetCollection> jets;
    iEvent.getByToken(src_, jets);
    
    edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);

    unsigned int ncand = 0;
    std::vector<int> partonFlavour;
    std::vector<int> partonIdx;
    std::vector<int> partonMotherIdx;
    std::vector<int> partonMotherPdgId;
    std::vector<uint8_t> hadronFlavour;
    int motherIdx;
    int motherPdgId;
    
    for (const reco::GenJet & jet : *jets) {
      if (!cut_(jet)) continue;
      ++ncand;
      bool matched = false;
      for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
        if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < deltaR_) {
          partonFlavour.push_back(jetFlavourInfoMatching.second.getPartonFlavour());

    if (jetFlavourInfoMatching.second.getPartons().size()!= 0){
    partonIdx.push_back(jetFlavourInfoMatching.second.getPartons().at(0).key());
	  //std::cout << "Partons size " << jetFlavourInfoMatching.second.getPartons().size()<< std::endl;
    //std::cout << "Parton pdgId " << jetFlavourInfoMatching.second.getPartons().at(0)->pdgId()<< std::endl;
    //std::cout << "Parton key " << jetFlavourInfoMatching.second.getPartons().at(0).key()<< std::endl;
    }

    else{
      partonIdx.push_back(-1);
      //std::cout << "Size = 0"<<std::endl;
    }

        // if (jetFlavourInfoMatching.second.getPartons().size()!=0){
      // partonIdx.push_back(jetFlavourInfoMatching.second.getPartons().at(0));
    // }
    // else{
      // partonIdx.push_back(-1);
    // }
// 

	  //If parton has a mother, fetch idx and pdg Id information
	  if (jetFlavourInfoMatching.second.getPartons().size()!= 0 and jetFlavourInfoMatching.second.getPartons().at(0)->numberOfMothers()>0){ 
		motherIdx = int(jetFlavourInfoMatching.second.getPartons().at(0)->motherRef(0).key());
		motherPdgId = int(jetFlavourInfoMatching.second.getPartons().at(0)->motherRef(0)->pdgId());
	  }
	  else{
	   motherIdx = -1;
	   motherPdgId = -1;
	  }
//	  std::cout << motherIdx << std::endl;
          partonMotherIdx.push_back(motherIdx);
          partonMotherPdgId.push_back(motherPdgId);
//	  std::cout << motherIdx << std::endl;

          hadronFlavour.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
          matched = true;
          break;
        }
      }
      if (!matched) {
        partonFlavour.push_back(0);
        partonIdx.push_back(0);
        partonMotherIdx.push_back(0);
        partonMotherPdgId.push_back(0);
        hadronFlavour.push_back(0);
      }
    }

    auto tab  = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
    tab->addColumn<int>("partonFlavour", partonFlavour, "flavour from parton matching");
    tab->addColumn<int>("partonIdx", partonIdx, "GenPart idx from parton matching");
    tab->addColumn<uint8_t>("hadronFlavour", hadronFlavour, "flavour from hadron ghost clustering");
    tab->addColumn<int>("partonMotherIdx", partonMotherIdx, "parton mother Idx");
    tab->addColumn<int>("partonMotherPdgId", partonMotherPdgId, "parton mother pdgId");

    iEvent.put(std::move(tab));
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(GenJetFlavourTableProducer);