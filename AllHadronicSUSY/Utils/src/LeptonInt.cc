// -*- C++ -*-
//
// Package:    LeptonInt
// Class:      LeptonInt
// 
/**\class LeptonInt LeptonInt.cc RA2Classic/LeptonInt/src/LeptonInt.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Fri Apr 11 16:35:33 CEST 2014
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

class LeptonInt : public edm::EDProducer {
public:
	explicit LeptonInt(const edm::ParameterSet&);
	~LeptonInt();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
        virtual bool EleVeto(const pat::Electron aEle,const reco::Vertex vtx);	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	std::vector<edm::InputTag> leptonTag_;
	edm::InputTag  eleTag_,muonTag_,PrimVtxTag_;
	
	
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
LeptonInt::LeptonInt(const edm::ParameterSet& iConfig)
{
	//register your produc

	leptonTag_ =iConfig.getParameter< std::vector<edm::InputTag> >("LeptonTag");
	PrimVtxTag_=iConfig.getParameter<edm::InputTag>("srcPV");
	eleTag_=iConfig.getParameter<edm::InputTag>("srcEle"); 
	muonTag_=iConfig.getParameter<edm::InputTag>("srcMuon");
	produces<int>("Leptons");
	produces<int>("Electrons");
        produces<int>("Muons");
	/* Examples
	 *   produces<ExampleData2>();
	 * 
	 *   //if do put with a label
	 *   produces<ExampleData2>("label");
	 * 
	 *   //if you want to put into the Run
	 *   produces<ExampleData2,InRun>();
	 */
	//now do what ever other initialization is needed
	
}


LeptonInt::~LeptonInt()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
LeptonInt::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	int Electrons=0;
        int Muons=0;
	edm::Handle<reco::VertexCollection> vtx_h;
        iEvent.getByLabel(PrimVtxTag_, vtx_h);

        edm::Handle<edm::View<pat::Electron> > eleHandle;
        iEvent.getByLabel(eleTag_, eleHandle);
	 if(eleHandle.isValid()){

           for(unsigned int e=0; e<eleHandle->size(); ++e){
         
           if(fabs(eleHandle->at(e).eta())>2.5 ||eleHandle->at(e).pt()<10)continue;
           	if(EleVeto(eleHandle->at(e), vtx_h->at(0)))++Electrons;
           }
         }

	edm::Handle<edm::View<pat::Muon> > muonHandle;
        iEvent.getByLabel(muonTag_, muonHandle);
        if(muonHandle.isValid()){
          for(unsigned int m=0; m<muonHandle->size(); ++m){
          if(muonHandle->at(m).pt()<10 || fabs(muonHandle->at(m).eta())>2.5)continue;
          float ChgIso=muonHandle->at(m).pfIsolationR04().sumChargedHadronPt;
          float ChgPU=muonHandle->at(m).pfIsolationR04().sumPUPt;
          float NeuIso=muonHandle->at(m).pfIsolationR04().sumNeutralHadronEt+
          muonHandle->at(m).pfIsolationR04().sumPhotonEt;
          float dBIsoMu= (ChgIso+std::max(0., NeuIso-0.5*ChgPU))/muonHandle->at(m).pt();
         if(muonHandle->at(m).isTightMuon( vtx_h->at(0)) && dBIsoMu<0.2)++Muons;

          }
        }

	int Leptons=0;
	for(unsigned int i=0; i< leptonTag_.size();i++)
	{
		edm::Handle< edm::View<reco::Candidate> > cands;
		iEvent.getByLabel(leptonTag_.at(i),cands);
		if( cands.isValid() ) 
		{
			Leptons+=cands->size();
		}
		else std::cout<<"LeptonIntProducer::Error tag invalid: "<<leptonTag_[i]<<std::endl;
	}

	std::auto_ptr<int> htp(new int(Leptons));
	iEvent.put(htp, "Leptons");
	
	std::auto_ptr<int> htp1(new int(Electrons));
        iEvent.put(htp1,"Electrons" );

        std::auto_ptr<int> htp2(new int(Muons));
        iEvent.put(htp2,"Muons" );
}

bool LeptonInt::EleVeto(const pat::Electron aEle, const reco::Vertex vtx){
bool vetoed=false;
//float eta=aEle.eta();
  // id variables
     float sieie         = aEle.full5x5_sigmaIetaIeta();
     bool convVeto       = aEle.passConversionVeto();
     int mhits 		 = aEle.gsfTrack()->numberOfLostHits();
     float dEtaIn        = aEle.deltaEtaSuperClusterTrackAtVtx();
     float dPhiIn        = aEle.deltaPhiSuperClusterTrackAtVtx();
     float hoe           = aEle.hadronicOverEm();
     float ooemoop       = (1.0/aEle.ecalEnergy() - aEle.eSuperClusterOverP()/aEle.ecalEnergy());
     float d0vtx         = 0.0;
     float dzvtx         = 0.0;
     reco::GsfElectron::PflowIsolationVariables pfIso = aEle.pfIsolationVariables();
    // // Compute isolation with delta beta correction for PU
    float absiso = pfIso.sumChargedHadronPt
     + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );

          d0vtx = aEle.gsfTrack()->dxy(vtx.position());
	  dzvtx = aEle.gsfTrack()->dz(vtx.position());
   absiso=absiso/aEle.pt(); 

if(aEle.isEB()){
	if(sieie< 0.011100 && fabs(dEtaIn)< 0.016315  && fabs(dPhiIn)< 0.252044 && hoe<0.345843 && absiso<0.164369 && ooemoop< 0.248070 && fabs(d0vtx)< 0.060279 && fabs(dzvtx)< 0.800538 && mhits<=2 && convVeto)vetoed=true;
}
else{
        if(sieie<  0.033987  && fabs(dEtaIn)<  0.010671   && fabs(dPhiIn)<  0.245263  && hoe< 0.134691  && absiso< 0.212604  && ooemoop<  0.157160  && fabs(d0vtx)<  0.273097  && fabs(dzvtx)<  0.885860  && mhits<=3 && convVeto)vetoed=true;
}



return vetoed;
}
// ------------ method called once each job just before starting event loop  ------------
void 
LeptonInt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonInt::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
LeptonInt::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
LeptonInt::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LeptonInt::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LeptonInt::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonInt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonInt);
