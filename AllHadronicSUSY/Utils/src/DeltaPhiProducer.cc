// -*- C++ -*-
//
// Package:    DeltaPhiProducer
// Class:      DeltaPhiProducer
// 
/**\class DeltaPhiProducer DeltaPhiProducer.cc RA2Classic/DeltaPhiProducer/src/DeltaPhiProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Thu Apr 17 10:49:52 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

class DeltaPhiProducer : public edm::EDProducer {
   public:
      explicit DeltaPhiProducer(const edm::ParameterSet&);
      ~DeltaPhiProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      edm::InputTag mhtTag_, mhtJetsTag_;
   double deltaPhi1, deltaPhi2, deltaPhi3;
   double minDeltaPhi;
   double minJetPt_, maxJetEta_;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
DeltaPhiProducer::DeltaPhiProducer(const edm::ParameterSet& iConfig)
{
  mhtJetsTag_ = iConfig.getParameter<edm::InputTag>("MHTJets");
  mhtTag_ = iConfig.getParameter<edm::InputTag>("MHT");
  minJetPt_    = iConfig.getParameter<double>("MinJetPt");
  maxJetEta_   = iConfig.getParameter<double>("MaxJetEta");
  const std::string string1("DeltaPhi1");
  const std::string string2("DeltaPhi2");
  const std::string string3("DeltaPhi3");
  const std::string string4("minDeltaPhi");
  produces<double> (string1).setBranchAlias(string1);
  produces<double> (string2).setBranchAlias(string2);
  produces<double> (string3).setBranchAlias(string3);
  produces<double> (string3).setBranchAlias(string4);
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


DeltaPhiProducer::~DeltaPhiProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DeltaPhiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  deltaPhi1=-9999;
  deltaPhi2=-9999;
  deltaPhi3=-9999;
  minDeltaPhi=9999;
  edm::Handle< edm::View<reco::Candidate> > mht;
  iEvent.getByLabel(mhtTag_,mht);
  edm::Handle< edm::View<reco::Candidate> > mhtJets;
  iEvent.getByLabel(mhtJetsTag_,mhtJets);
  int iJet=0;
  if( mhtJets.isValid() && mht.isValid() ) {
    for(unsigned int i=0; i < mhtJets->size(); i++)
    {
      double deltaPhi=-999;
      if(mhtJets->at(i).pt()>minJetPt_ && mhtJets->at(i).eta()<maxJetEta_ ) deltaPhi=std::abs(reco::deltaPhi(mhtJets->at(i).phi(),mht->at(0).phi()));
      if (iJet==0) { deltaPhi1=deltaPhi; iJet++; continue;}
      if (iJet==1) { deltaPhi2=deltaPhi; iJet++; continue;}
      if (iJet==2) { deltaPhi3=deltaPhi; iJet++; continue;}
      if(minDeltaPhi>deltaPhi)minDeltaPhi=deltaPhi;
      if (iJet==3) break;
    }
  }
  else std::cout<<"DeltaPhiProducer::Error input mhtJets or mhtMetTag not valid!!! mhtJetTag:"<<mhtJetsTag_.label()<<", mhtMetTag:"<<mhtTag_.label()<<std::endl;
 // std::cout<<"DletaPhi1"<<deltaPhi1<<", deltaPhi2"<<deltaPhi2<<", deltaPhi3"<<deltaPhi3<<std::endl;
  const std::string string1("DeltaPhi1");
  const std::string string2("DeltaPhi2");
  const std::string string3("DeltaPhi3");
  const std::string string4("minDeltaPhi");
  std::auto_ptr<double> htp1(new double(deltaPhi1));
  std::auto_ptr<double> htp2(new double(deltaPhi2));
  std::auto_ptr<double> htp3(new double(deltaPhi3));
  std::auto_ptr<double> htp4(new double(minDeltaPhi));
  iEvent.put(htp1,string1 );
  iEvent.put(htp2,string2 );
  iEvent.put(htp3,string3 );
  iEvent.put(htp4, string4);
}

// ------------ method called once each job just before starting event loop  ------------
void 
DeltaPhiProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DeltaPhiProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DeltaPhiProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DeltaPhiProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DeltaPhiProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DeltaPhiProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeltaPhiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeltaPhiProducer);
