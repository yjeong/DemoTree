// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Yongho Jeong
//         Created:  Thu, 04 Feb 2016 16:32:47 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <map>

// user include files
// -------------------------FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//-----------------------CommonTools
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//----------------------DataFormats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
//#include "DataFormats/TrackingRecHit/interface/TrackingRecoHitFwd.h"

#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimMuon/MCTruth/plugins/MuonAssociatorEDProducer.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <string>
#include <sstream>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

// class declaration

const int nMax=24000;//over than cmsRun events
int ngp, nm;
float mllpm;
float mllpmGen;
float Tgen_pt[nMax], Tgen_eta[nMax];
float mu_pt[nMax],mu_eta[nMax];
float mudR;

class DemoAnalyzer : public edm::EDAnalyzer {
	public:
		//Constructor
		explicit DemoAnalyzer(const edm::ParameterSet&);
		//Destructor
		~DemoAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		edm::EDGetTokenT<reco::TrackCollection> generalTracks;
		//edm::EDGetTokenT<reco::MuonCollection> muons;
		edm::EDGetTokenT<reco::TrackCollection> globalMuons;
		edm::EDGetTokenT<reco::GenParticleCollection> genParticles;
		//edm::EDGetTokenT<reco::VertexCollection> offlinePrimaryVertices;

		// ----------member data ---------------------------
		TFile *demo_tree;
		TH1D *demohisto;

		TTree *tree;
		unsigned int minTracks_;
		bool SaveHisto;
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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
	minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
	//	minTracks_(iConfig.getUntrackedParameter<edm::InputTag>("minTracks"))
{
	generalTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks"));
	globalMuons = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Mu"));
	//muons = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("Mu"));
	genParticles = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("Gen"));
	//offlinePrimaryVertices = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("Ver"));

	//now do what ever initialization is needed
	SaveHisto=iConfig.getParameter<bool>("SaveHisto");

	edm::Service<TFileService> fs;
	if(SaveHisto)demo_tree = new TFile("out.root","recreate");
	//if(SaveHisto) tree = fs->make<TTree>("PAT.root","PAT_tree");
	demohisto = new TH1D("tracks","Tracks",100,0,5000);

	tree = new TTree("tree","example_tree");
	tree->Branch("ngp",&ngp,"ngp/I");
	tree->Branch("Tgen_pt",Tgen_pt,"Tgen_pt[ngp]/F");
	tree->Branch("Tgen_eta",Tgen_eta,"Tgen_eta[ngp]/F");
	tree->Branch("mllpmGen",&mllpmGen,"mllpmGen/F");

	tree->Branch("nm",&nm,"nm/I");
	tree->Branch("mu_pt",mu_pt,"mu_pt[nm]/F");
	tree->Branch("mu_eta",mu_eta,"mu_eta[nm]/F");
	tree->Branch("mudR",&mudR,"mudR/F");
	tree->Branch("mllpm",&mllpm,"mllpm/F");
}

DemoAnalyzer::~DemoAnalyzer()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
	void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	vector<float> vpt;

	/*edm::Handle<reco::VertexCollection> Ver;
	iEvent.getByToken(offlinePrimaryVertices,Ver);
	if(Ver->empty()){cout<<"no PV"<<endl;return;}
	auto pv0 = Ver->front();*/

	edm::Handle<reco::GenParticleCollection> Gen;
	iEvent.getByToken(genParticles,Gen);

	reco::CandidateCollection GenCandidate;

	//double minDRP=100;//set variable of minimum "DeltaR" between two particle.
	//double minDRM=100;//set variable of minimum "DeltaR" between two particle.
	unsigned int ngen=0;
	int nGenP=-1, nGenM=-1;
	double geptP=0,geptM=0;
	for(auto gen = Gen->begin();gen != Gen->end(); gen++, ngen++)
	{
		int ndau = gen->numberOfDaughters();
		int id = gen->pdgId();
		int gstat = gen->status();
		if(fabs(gen->pt())>10 && fabs(gen->eta())<2.4 && abs(id)==13 && abs(ndau)!=0 && abs(gstat)!=1)//or status>1
		{
			vpt.push_back(gen->pt());
			if(abs(gen->mother(0)->pdgId())==23)
			{
				if(gen->charge()==1) nGenP=ngen;
				if(gen->charge()==-1) nGenM=ngen;
			}
			if(gen->charge()==1 && geptP<gen->pt()) geptP=gen->pt();
			if(gen->charge()==-1 && geptM<gen->pt()) geptM=gen->pt();
		}
		Tgen_pt[ngen] = gen->pt();
		Tgen_eta[ngen] = gen->eta();
		if(abs(id)==13 && abs(ndau)!=0 && abs(gstat)!=1 && (abs(gen->mother(0)->pdgId())==23))	GenCandidate.push_back(*gen);
	}
	if(GenCandidate.size()==2) mllpmGen=(GenCandidate[0].p4()+GenCandidate[1].p4()).M();
	ngp=ngen;

	TLorentzVector GenP, GenM;
	if(nGenP!=-1 && nGenM!=-1)
	{
		//cout<<nGenP<<endl;
		GenP.SetPxPyPzE((*Gen)[nGenP].px(),(*Gen)[nGenP].py(),(*Gen)[nGenP].pz(),(*Gen)[nGenP].energy());
		GenM.SetPxPyPzE((*Gen)[nGenM].px(),(*Gen)[nGenM].py(),(*Gen)[nGenM].pz(),(*Gen)[nGenM].energy());
	}

	edm::Handle<reco::TrackCollection> Tracks;
	iEvent.getByToken(generalTracks,Tracks);

	if (minTracks_ <=Tracks->size()) {
		LogInfo("Demo")<<"number of Tracks"<<Tracks->size();
	}
	demohisto->Fill(Tracks->size());

	//edm::Handle<reco::MuonCollection> Mu;
	//iEvent.getByToken(muons,Mu);

	edm::Handle<reco::TrackCollection> Mu;
	iEvent.getByToken(globalMuons,Mu);

	unsigned int nmu=0;
	int nmuP = -1, nmuM =-1;
	double muPt_P = 0, muPt_M = 0;
	for(auto mu = Mu->begin(); mu != Mu->end(); mu++, nmu++){
		//V_muons[nmu] = mu->momentum().x(), mu->momentum().y(), mu->momentum().z(), mu->energy();
		if(mu->charge()==1 && muPt_P<mu->pt()){
			nmuP = nmu;
			muPt_P = mu->pt();
		}
		if(mu->charge()==-1 && muPt_M<mu->pt()){
			nmuM = nmu;
			muPt_M = mu->pt();
		}
		mu_pt[nmu] = mu->pt();
		mu_eta[nmu] = mu->eta();
	}
	nm=nmu;

	TLorentzVector muP, muM;
	if(nmuP != -1 && nmuM != 1){
		muP.SetPxPyPzE((*Mu)[nmuP].px(),(*Mu)[nmuP].py(),(*Mu)[nmuP].pz(),(*Mu)[nmuP].p());
		muM.SetPxPyPzE((*Mu)[nmuM].px(),(*Mu)[nmuM].py(),(*Mu)[nmuM].pz(),(*Mu)[nmuM].p());
		double dR = muP.DeltaR(muM);
		mudR = dR;
		mllpm = (muP+muM).M();
	}


	tree->Fill();
}
//-----------------------------------------------------

// ------------ method called once each job just before starting event loop  ------------
	void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
DemoAnalyzer::endJob() 
{


	if(SaveHisto)
	{
		demo_tree->cd();
		demohisto->Write();
		tree->Write();

		demo_tree->Close();
	}
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
