// system include files
#include <memory>
#include <iomanip>
#include <iostream>
#include <string>
#include <set>
#include <utility>
#include <algorithm>  // std::sort, std::swap
#include <math.h>
#include <unordered_set>
#include <utility>

// ROOT
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>



// CMS Framework (FWCore)
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


// DataFormats
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


// Vertexing 
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// Other
#include "PhysicsTools/PatAlgos/plugins/PATSingleVertexSelector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" //for pileup 
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

using namespace pat;

using pat::PATSingleVertexSelector;

using namespace std;
using namespace edm;
using namespace reco;

class testJpsi_AOD : public edm::EDAnalyzer {
 public:
	explicit testJpsi_AOD(const edm::ParameterSet&);
      ~testJpsi_AOD();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
   virtual void beginJob();
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();

 
  std::map<std::string,TProfile*> tprofile_;


 


 //Token for AOD

 virtual void beginRun(edm::Run const&, edm::EventSetup const&);
 virtual void endRun(edm::Run const&, edm::EventSetup const&);
 virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


 int i;
 int run_number = 0, lumi = 0, event_number=0;

 edm::EDGetToken tracks_;
 edm::EDGetToken hVtx_;
 edm::EDGetToken beamSpotHandle_;
 edm::EDGetToken trigger;

 edm::EDGetTokenT<reco::MuonCollection> muonCollToken;
 edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;
 edm::EDGetTokenT<BXVector<l1t::Muon>> l1MuonCollToken;
 edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
 //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
   
  struct pair_hash {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ h2;
        }
    };


  TTree * Tree;
  TTree * Tree2;
  TTree * Tree3;
  TTree * Tree4;
  TTree * treemc;


  //from J/psi
  TH1F* h_RecDiMuonM;
  TH1F* h_RecDiMuonM2;
  TH1F* h_GenDiMuonM;
  TH1F* h_MassRes;
  TH1F* h_MupTRes;

   
    TLorentzVector jpsiP4;


    TLorentzVector jpsi2P4;

    
    reco::Track tr1;
    reco::Track tr2;
   
    int vertex_size, num_tracks;
    bool passTrigger;
    double mu1_pt, mu1_eta, mu1_phi, mu1_mass;
    double amu1_pt, amu1_eta, amu1_phi, amu1_mass; 
    double mu2_pt, mu2_eta, mu2_phi, mu2_mass;
    double amu2_pt, amu2_eta, amu2_phi, amu2_mass; 
    double Psireco1_pt, Psireco1_eta, Psireco1_phi, Psireco1_mass;
    double Psireco2_pt, Psireco2_eta, Psireco2_phi, Psireco2_mass;
    double ctau, ctau2; 
    bool is_Jpsiprompt = false;
    bool is_Jpsiprompt2 = false;

    double PV_x = 0;
    double PV_y = 0;

    double PV2_x = 0;
    double PV2_y = 0;
    
    double r_xy = 0;
    double L_xy = 0;
    
    double x_diff = 0;
    double y_diff = 0;

    double r_xy2 = 0;
    double L_xy2 = 0;
  
    double x_diff2 = 0;
    double y_diff2 = 0;

    double x_diff3 = 0;
    double y_diff3 = 0;

    double r_xy3 = 0;

    double CL = 0;
    double CL2 = 0;

    double Ctau_ok = 0;
    double Ctau2_ok = 0;

    double Jpsi1_x = 0;
    double Jpsi1_y = 0;

    double Jpsi2_x = 0;
    double Jpsi2_y = 0;


    double fourmuons_x = 0;
    double fourmuons_y = 0;


    double Psireco1_px = 0;
    double Psireco1_py  = 0;
    double Psireco2_px = 0;
    double Psireco2_py  = 0;
    double scalar_prod1 = 0;
    double scalar_prod2 = 0;
    double CL_vertex = 0;
 

    int sv=0;
    std::vector<double> lepton1_fromPV, lepton2_fromPV, lepton3_fromPV, lepton4_fromPV;
    std::vector<double> lepton1_pvassq, lepton2_pvassq, lepton3_pvassq, lepton4_pvassq;

    double lepton1frompv = -10; double lepton2frompv = -10; double lepton3frompv = -10; double lepton4frompv = -10;
    double lepton1pvassq = -10; double lepton2pvassq = -10; double lepton3pvassq = -10; double lepton4pvassq = -10;
    int lepton1PVIndex = -1; int lepton2PVIndex = -1; int lepton3PVIndex = -1;
    std::vector<std::string> triggerlist;

};


 testJpsi_AOD:: testJpsi_AOD(const edm::ParameterSet& iConfig)
{

edm::Service<TFileService> fs;

//from J/psi





   Tree = fs->make<TTree>("tree", "tree");
   Tree2 = fs->make<TTree>("tree2", "tree2");
   Tree3 = fs->make<TTree>("tree3", "tree3");
   Tree4 = fs->make<TTree>("tree4", "tree4");

   Tree->Branch("event_number", &event_number);
   Tree->Branch("run_number", &run_number);
   Tree->Branch("lumi", &lumi, "lumi/I");

   
   Tree->Branch("vertex_size",&vertex_size,"vertex_size/I");
   Tree->Branch("passTrigger", &passTrigger, "passTrigger/B"); 

   Tree->Branch("mu1_pt",&mu1_pt,"mu1_pt/D"); //muon
   Tree->Branch("mu1_eta",&mu1_eta,"mu1_eta/D");
   Tree->Branch("mu1_phi",&mu1_phi,"mu1_phi/D");
   Tree->Branch("mu1_mass",&mu1_mass,"mu1_mass/D");

   Tree->Branch("amu1_pt",&amu1_pt,"amu1_pt/D"); //anti-muon
   Tree->Branch("amu1_eta",&amu1_eta,"amu1_eta/D");
   Tree->Branch("amu1_phi",&amu1_phi,"amu1_phi/D");
   Tree->Branch("amu1_mass",&amu1_mass,"amu1_mass/D");

   Tree2->Branch("mu2_pt",&mu2_pt,"mu2_pt/D"); //muon
   Tree2->Branch("mu2_eta",&mu2_eta,"mu2_eta/D");
   Tree2->Branch("mu2_phi",&mu2_phi,"mu2_phi/D");
   Tree2->Branch("mu2_mass",&mu2_mass,"mu2_mass/D");

   Tree2->Branch("amu2_pt",&amu2_pt,"amu2_pt/D"); //anti-muon
   Tree2->Branch("amu2_eta",&amu2_eta,"amu2_eta/D");
   Tree2->Branch("amu2_phi",&amu2_phi,"amu2_phi/D");
   Tree2->Branch("amu2_mass",&amu2_mass,"amu2_mass/D");

   Tree->Branch("Psireco1_pt", &Psireco1_pt, "Psireco1_pt/D");
   Tree->Branch("Psireco1_eta", &Psireco1_eta, "Psireco1_eta/D");
   Tree->Branch("Psireco1_phi", &Psireco1_phi, "Psireco1_phi/D");
   Tree->Branch("Psireco1_mass", &Psireco1_mass, "Psireco1_mass/D");

   Tree2->Branch("Psireco2_pt", &Psireco2_pt, "Psireco2_pt/D");
   Tree2->Branch("Psireco2_eta", &Psireco2_eta, "Psireco2_eta/D");
   Tree2->Branch("Psireco2_phi", &Psireco2_phi, "Psireco2_phi/D");
   Tree2->Branch("Psireco2_mass", &Psireco2_mass, "Psireco2_mass/D");

   Tree->Branch("lepton1_fromPV", &lepton1_fromPV);
   Tree->Branch("lepton1_pvassq", &lepton1_pvassq);
   Tree->Branch("lepton2_fromPV", &lepton2_fromPV);
   Tree->Branch("lepton2_pvassq", &lepton2_pvassq);
   Tree->Branch("lepton3_fromPV", &lepton3_fromPV);
   Tree->Branch("lepton3_pvassq", &lepton3_pvassq);
   Tree->Branch("lepton4_fromPV", &lepton4_fromPV);
   Tree->Branch("lepton4_pvassq", &lepton4_pvassq);


   Tree->Branch("is_Jpsiprompt", &is_Jpsiprompt, "is_Jpsiprompt/B");
   Tree->Branch("is_Jpsiprompt2", &is_Jpsiprompt2, "is_Jpsiprompt2/B");
   Tree->Branch("Jpsi1_x", &Jpsi1_x, "Jpsi1_x/D");
   Tree->Branch("Jpsi1_y", &Jpsi1_y, "Jpsi1_y/D");
   Tree2->Branch("Jpsi2_x", &Jpsi2_x, "Jpsi2_x/D"); 
   Tree2->Branch("Jpsi2_y", &Jpsi2_y, "Jpsi2_y/D");

   Tree->Branch("PV_x", &PV_x, "PV_x/D");
   Tree->Branch("PV_y", &PV_y, "PV_y/D");
   Tree->Branch("Ctau_ok", &Ctau_ok, "Ctau_ok/D");
   Tree2->Branch("PV2_x", &PV2_x, "PV2_x/D");
   Tree2->Branch("PV2_y", &PV2_y, "PV2_y/D");
   Tree2->Branch("Ctau2_ok", &Ctau2_ok, "Ctau2_ok/D");
   


   /* Tree3->Branch("lepton1_fromPV", &lepton1_fromPV, "lepton1_fromPV/D");
    Tree3->Branch("lepton1_pvassq", &lepton1_pvassq, "lepton1_pvassq/D");
    Tree3->Branch("lepton2_fromPV", &lepton2_fromPV, "lepton2_fromPV/D");
    Tree3->Branch("lepton2_pvassq", &lepton2_pvassq, "lepton2_pvassq/D");
    Tree3->Branch("lepton3_fromPV", &lepton3_fromPV, "lepton3_fromPV/D");
    Tree3->Branch("lepton3_pvassq", &lepton3_pvassq, "lepton3_pvassq/D");
    Tree3->Branch("lepton4_fromPV", &lepton4_fromPV, "lepton4_fromPV/D");
    Tree3->Branch("lepton4_pvassq", &lepton4_pvassq, "lepton4_pvassq/D");*/
    Tree3->Branch("sv", &sv, "sv/I");
    Tree3->Branch("ctau", &ctau, "ctau/D");
    Tree3->Branch("ctau2", &ctau2, "ctau2/D");
    Tree3->Branch("fourmuons_x", &fourmuons_x, "fourmuons_x/D");
    Tree3->Branch("fourmuons_y", &fourmuons_y, "fourmuons_y/D");

    Tree3->Branch("r_xy3", &r_xy3, "r_xy3/D");
    Tree3->Branch("x_diff3", &x_diff3, "x_diff3/D");
    Tree3->Branch("y_diff3", &y_diff3, "y_diff3/D");

    Tree4->Branch("CL_vertex", &CL_vertex, "CL_vertex/D");
    
  
tracks_= consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
hVtx_= consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
beamSpotHandle_= consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
trigger = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
muonCollToken = consumes<reco::MuonCollection>(edm::InputTag("muons"));
genCollToken = consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"));
l1MuonCollToken  = consumes<BXVector<l1t::Muon>>(edm::InputTag ("gmtStage2Digis:Muon"));




}

 testJpsi_AOD:: ~testJpsi_AOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}



void testJpsi_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  using MuonPair = std::pair<size_t, size_t>;
  std::unordered_set<MuonPair, pair_hash> usedMuonPairs;
  



 edm::Handle<reco::TrackCollection> tracks;
 iEvent.getByToken(tracks_ , tracks );

 reco::BeamSpot beamSpot;
 edm::Handle<reco::BeamSpot> beamSpotHandle;
 iEvent.getByToken(beamSpotHandle_ , beamSpotHandle);


 edm::ESHandle<TransientTrackBuilder> theB;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


   //from J/psi

  edm::Handle<vector<reco::Muon>> muons;
  iEvent.getByToken(muonCollToken, muons);

  bool muonsize_minimum = false;

  if (muons->size()>3) {
    muonsize_minimum = true;
  }

  edm::Handle <pat::PackedGenParticleCollection> genColl;
  iEvent.getByToken(genCollToken, genColl);

  //edm::Handle<reco::VertexCollection> vertices;
  //iEvent.getByToken(token_vertices, vertices);


// *******
// TRIGGER
// *******

edm::Handle<edm::TriggerResults> trigResults;
 iEvent.getByToken(trigger, trigResults);

const edm::TriggerNames& names = iEvent.triggerNames(*trigResults);
bool passTrig=false;
bool mc = false;


//from J/psi
//edm::Handle<edm::TriggerResults> triggerResults;
//iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), triggerResults);

//const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
//bool passTrig=false;

for (unsigned int i = 0; i < names.size(); i++) {
   //std::cout << "Trigger name: " << names.triggerName(i) << std::endl;
   if (names.triggerName(i).find("HLT_Dimuon0_Jpsi3p5_Muon2")) {
  //if (names.triggerName(i).find("HLT_Dimuon0_Jpsi_Muon90_trigger")) {
        if (trigResults->accept(names.triggerIndex(names.triggerName(i)))) {
            passTrig = true;
           // cout<<"PASSATO TRIGGER"<<endl;
        }
    }
}

// ******
// VERTEX
// ******

//from J/psi

 edm::Handle<std::vector<reco::Vertex>> hVtx;
 iEvent.getByToken(hVtx_ , hVtx); 
 reco::Vertex primVertex;
 if(hVtx->size() > 0){
 TransientVertex v;
 reco::Vertex primaryVertex;
 reco::Vertex primVertex_tmp;
 bool goodvertex = false;

 double bsx = 0, bsy = 0, bsz = 0;
  if ( beamSpotHandle.isValid() ){ beamSpot = *beamSpotHandle;
   bsx = beamSpot.x0();
   bsy = beamSpot.y0();
   bsz = beamSpot.z0();
  }
  
  GlobalPoint BeamSpotGP(bsx, bsy, bsz);
  primaryVertex = hVtx->at(0);
  //const reco::Vertex* firstGoodVertex;
  for(unsigned int t = 0; t<hVtx->size(); t++){
   primVertex_tmp = hVtx->at(t);
   if(!primVertex_tmp.isFake() && primVertex_tmp.isValid() && primVertex_tmp.ndof() > 4 && fabs(primVertex_tmp.z()-bsz) < 10 ){
    goodvertex = true;
    primaryVertex = primVertex_tmp;
    //firstGoodVertex = &primVertex_tmp;
    break;
   } 
  }

if (!goodvertex || !passTrig || !muonsize_minimum) return;  //overall condition


//J/psi reconstruction

/////////////////////////////////////////////////
// Dimuon pairs /////////////////////////////////
/////////////////////////////////////////////////


    std::vector<reco::Muon> muonsPositiveCharge;
    std::vector<reco::Muon> muonsNegativeCharge;

    
    TLorentzVector Jpsi_candidate;
    std::vector<TLorentzVector> Jpsivector;


    
    reco::Muon highestPtMuPlus1;
    reco::Muon highestPtMuMinus1;
    reco::Muon highestPtMuPlus2;
    reco::Muon highestPtMuMinus2;

   
    std::vector<reco::TransientTrack> tks2; 
      
    std::vector<reco::TransientTrack> tks2_save1; 
    std::vector<reco::TransientTrack> tks2_save2; 

    std::vector<reco::TransientTrack> tks2_save3;
    
    reco::TransientTrack  transientMuonTrack;
    reco::TransientTrack  transientMuonTrack2;

    TransientVertex v2;
    TransientVertex v2_check;
    std::vector<double> vJpsi1;
    std::vector<double> vJpsi2;
    int sameVertex = 0;
    bool goodtracks=false;

    bool first_found=false;




std::unordered_set<MuonPair, pair_hash> usedMuonPairs;
std::unordered_set<int> usedMuons;


for (auto mup = muons->cbegin(); mup != muons->cend(); ++mup) {
    if (!mup->isGlobalMuon()) continue;
    if (not (mup->charge() > 0)) continue;
    if ((mup->pt() < 2.5 && (fabs(mup->eta()) < 1.2 && fabs(mup->eta()) > 2.4)) || (mup->pt() < 3.5 && fabs(mup->eta()) >= 1.2)) continue;

    for (auto mum = muons->cbegin(); mum != muons->cend(); ++mum) {
        if (!mum->isGlobalMuon()) continue;
        if (not (mum->charge() < 0)) continue;
        if ((mum->pt() < 2.5 && (fabs(mup->eta()) < 1.2 && fabs(mup->eta()) > 2.4))  || (mum->pt() < 3.5 && fabs(mum->eta()) >= 1.2)) continue;

        int mupIndex = std::distance(muons->cbegin(), mup);
        int mumIndex = std::distance(muons->cbegin(), mum);

        // Controlliamo se uno dei muoni è già stato utilizzato
        if (usedMuons.count(mupIndex) > 0 || usedMuons.count(mumIndex) > 0) continue;

        // Controllo se la coppia attuale è già stata usata
        auto currentPair = std::make_pair(mupIndex, mumIndex);
        if (usedMuonPairs.count(currentPair) > 0) continue;

        // Codice per la coppia di muoni
        TLorentzVector muPlus, muMinus;
        muPlus.SetPtEtaPhiM(mup->pt(), mup->eta(), mup->phi(), mup->mass());
        muMinus.SetPtEtaPhiM(mum->pt(), mum->eta(), mum->phi(), mum->mass());

        double diMuonRecMass = (mup->p4() + mum->p4()).M();
        double diMuonRecEta = (mup->p4() + mum->p4()).Eta();
        double diMuonRecPt = (mup->p4() + mum->p4()).Pt();
        if (diMuonRecMass >= 2.7 && diMuonRecMass <= 3.4 && fabs(diMuonRecEta) < 2.4 && diMuonRecPt > 2.9) {
        //if (diMuonRecMass >= 3.0 && diMuonRecMass <= 3.2 && fabs(diMuonRecEta) < 2.4 && diMuonRecPt > 2.9) {
            goodtracks = true;
            sameVertex = 0;      
            reco::TrackRef muonTrack = mup->track();
            reco::TrackRef muonTrack2 = mum->track(); 
            transientMuonTrack = ((*theB).build(muonTrack));
            transientMuonTrack2 = ((*theB).build(muonTrack2)); 

            tks2.clear();
            tks2.push_back(transientMuonTrack);
            tks2.push_back(transientMuonTrack2);

            if (tks2.size() > 1) {
                KalmanVertexFitter kalman(true);
                if (goodtracks) {
                    v2 = kalman.vertex(tks2);
                    if (v2.isValid()) {
                        sameVertex++;
                        CL = TMath::Prob(v2.totalChiSquared(), (int)v2.degreesOfFreedom());
                        Tree4->Fill();
                        if (CL > 0.01) {
                            

                            if (first_found == false) {  
                                tks2_save1.clear();
                                tks2_save1.push_back(transientMuonTrack);
                                tks2_save1.push_back(transientMuonTrack2);

                                vJpsi1.push_back(v2.position().x());
                                vJpsi1.push_back(v2.position().y());
                                vJpsi1.push_back(v2.position().z());

                                first_found = true;
                                // Aggiungiamo gli indici dei muoni e la coppia all'insieme degli usati
                                usedMuons.insert(mupIndex);
                                usedMuons.insert(mumIndex);
                                usedMuonPairs.insert(currentPair);

                                Jpsi_candidate = muPlus + muMinus;
                                Jpsivector.push_back(Jpsi_candidate);
                                mu1_pt = muMinus.Pt();
                                mu1_eta = muMinus.Eta();
                                mu1_phi = muMinus.Phi();
                                mu1_mass = muMinus.M();
                                amu1_pt = muPlus.Pt();
                                amu1_eta = muPlus.Eta();
                                amu1_phi = muPlus.Phi();
                                amu1_mass = muPlus.M();
                                Psireco1_pt = Jpsi_candidate.Pt();
                                Psireco1_eta = Jpsi_candidate.Eta();
                                Psireco1_phi = Jpsi_candidate.Phi();
                                Psireco1_mass = Jpsi_candidate.M();

                                x_diff = v2.position().x() - primaryVertex.position().x();
                                y_diff = v2.position().y() - primaryVertex.position().y();
                                //r_xy = sqrt(x_diff * x_diff + y_diff * y_diff);
                                Psireco1_px = Psireco1_pt * cos(Psireco1_phi);
                                Psireco1_py = Psireco1_pt * sin(Psireco1_phi);
                                scalar_prod1 = Psireco1_px * x_diff + Psireco1_py * y_diff;
                                L_xy = scalar_prod1 * Psireco1_pt / std::abs(Psireco1_pt);
                                ctau = L_xy * Psireco1_mass / Psireco1_pt;
                                is_Jpsiprompt = (ctau < 0.2 && ctau > -0.05);

                                PV_x = primaryVertex.position().x();
                                PV_y = primaryVertex.position().y();

                                Jpsi1_x = v2.position().x();
                                Jpsi1_y = v2.position().y();

                                Ctau_ok = ctau;

                                //Tree->Fill(); 
                                //std::cout << "index used for the first: " << mupIndex << " " << mumIndex << std::endl;

                            } else if (first_found == true) { 
                                tks2_save2.clear();
                                tks2_save2.push_back(transientMuonTrack);
                                tks2_save2.push_back(transientMuonTrack2);

                                vJpsi2.push_back(v2.position().x());
                                vJpsi2.push_back(v2.position().y());
                                vJpsi2.push_back(v2.position().z());

                                Jpsi_candidate = muPlus + muMinus; 
                                Jpsivector.push_back(Jpsi_candidate);
                                mu2_pt = muMinus.Pt();
                                mu2_eta = muMinus.Eta();
                                mu2_phi = muMinus.Phi();
                                mu2_mass = muMinus.M();
                                amu2_pt = muPlus.Pt();
                                amu2_eta = muPlus.Eta();
                                amu2_phi = muPlus.Phi();
                                mu2_mass = muPlus.M();
                                Psireco2_pt = Jpsi_candidate.Pt();
                                Psireco2_eta = Jpsi_candidate.Eta();
                                Psireco2_phi = Jpsi_candidate.Phi();
                                Psireco2_mass = Jpsi_candidate.M();
                                /*
                                std::cout << "Ps2 mass: " << Psireco2_mass << std::endl;
                                std::cout << "Ps2 pt: " << Psireco2_pt << std::endl;
                                std::cout << "index used for the second: " << mupIndex << " " << mumIndex << std::endl;
                                std::cout << "mu2 pt: " << mu2_pt << std::endl;
                                std::cout << "mu2 eta: " << mu2_eta << std::endl;
                                std::cout << "mu2 phi: " << mu2_phi << std::endl;
                               
                                std::cout << "amu2 pt: " << amu2_pt << std::endl;
                                std::cout << "amu2 eta: " << amu2_eta << std::endl;
                                std::cout << "amu2 phi: " << amu2_phi << std::endl;

                                std::cout <<"mu1 pt: " << mu1_pt << std::endl;
                                std::cout <<"mu1 eta: " << mu1_eta << std::endl;
                                std::cout <<"mu1 phi: " << mu1_phi << std::endl;
                                std::cout <<"amu1 pt: " << amu1_pt << std::endl;
                                std::cout <<"amu1 eta: " << amu1_eta << std::endl;
                                std::cout <<"amu1 phi: " << amu1_phi << std::endl;*/

                              


                                x_diff2 = v2.position().x() - primaryVertex.position().x();
                                y_diff2 = v2.position().y() - primaryVertex.position().y();
                                r_xy2 = sqrt(x_diff2 * x_diff2 + y_diff2 * y_diff2);
                                Psireco2_px = Psireco2_pt * cos(Psireco2_phi);
                                Psireco2_py = Psireco2_pt * sin(Psireco2_phi);
                                scalar_prod2 = Psireco2_px * x_diff + Psireco2_py * y_diff;
                                L_xy2 = scalar_prod2 * Psireco2_pt / std::abs(Psireco2_pt);
                                
                                ctau2 = L_xy2 * Psireco2_mass / Psireco2_pt;
                                is_Jpsiprompt2 = (ctau2 < 0.2 && ctau2 > -0.05);

                                PV2_x = primaryVertex.position().x();
                                PV2_y = primaryVertex.position().y();

                                Jpsi2_x = v2.position().x();
                                Jpsi2_y = v2.position().y();

                                Ctau2_ok = ctau2;


                                bool twofound = true;
                                //Tree2->Fill();  
                                tks2_save3.clear();
                                tks2_save3.push_back(tks2_save1[0]);
                                tks2_save3.push_back(tks2_save1[1]);
                                tks2_save3.push_back(tks2_save2[0]);
                                tks2_save3.push_back(tks2_save2[1]);

                                if (tks2_save3.size() > 3) {
                                    KalmanVertexFitter kalman2(true);
                                    if (goodtracks && twofound && is_Jpsiprompt && is_Jpsiprompt2) {
                                        v2_check = kalman2.vertex(tks2_save3);
                                        if (v2_check.isValid()) {
                                            sameVertex++;
                                            CL2 = TMath::Prob(v2_check.totalChiSquared(), (int)v2_check.degreesOfFreedom());
                                            if (CL2 > 0.01) {
                                                sv++;
                                                fourmuons_x = v2_check.position().x();
                                                fourmuons_y = v2_check.position().y();
                                                x_diff3 =  v2_check.position().x() - primaryVertex.position().x();
                                                y_diff3 =  v2_check.position().y() - primaryVertex.position().y();
                                                r_xy3 = sqrt(x_diff3 * x_diff3 + y_diff3 * y_diff3);
                                                Tree->Fill();
                                                Tree2->Fill();
                                                Tree3->Fill();
                                                ctau = 0;
                                                ctau2 = 0;
                                                tks2_save1.clear();
                                                tks2_save2.clear();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
}
}

void 
testJpsi_AOD::beginJob()
{

edm::Service<TFileService> fs;

//variabili vertice

}

// ------------ method called once each job just after ending the event loop  ------------
void 
testJpsi_AOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
testJpsi_AOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
testJpsi_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
testJpsi_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
testJpsi_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
testJpsi_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE (testJpsi_AOD);

