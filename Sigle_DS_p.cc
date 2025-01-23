// system include files
#include <memory>
#include <iomanip>
#include <iostream>
#include <string>
#include <set>
#include <utility>
#include <algorithm>  // std::sort, std::swap
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TTree.h"
#include "TProfile.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Version/interface/GetReleaseVersion.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <DataFormats/BeamSpot/interface/BeamSpot.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/Common/interface/Ref.h>
#include <DataFormats/Math/interface/deltaPhi.h>
#include <DataFormats/Common/interface/Ptr.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "Math/VectorUtil.h"

#include "PhysicsTools/PatAlgos/plugins/PATSingleVertexSelector.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include <TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h>
#include <TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h>
#include <TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h>

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/ConfigurableTrimmedVertexFinder.h"

#include "RecoTracker/TransientTrackingRecHit/interface/ProjectedRecHit2D.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit1DMomConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit2DPosConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiPixelRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripMatchedRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit1D.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

using pat::PATSingleVertexSelector; //AOD

using namespace std;
using namespace edm;
using namespace reco;




class testDstar_AOD : public edm::EDAnalyzer {
public:
    explicit testDstar_AOD(const edm::ParameterSet&);
    ~testDstar_AOD();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();
   
   //new for AOD 

   edm::InputTag assocTags_;
   std::map<std::string,TH1F*> histContainer_;

   std::map<std::string,TH2F*> histContainer2_;

   std::map<std::string,TProfile*> tprofile_;

  //Token for AOD


   virtual void beginRun(edm::Run const&, edm::EventSetup const&);
   virtual void endRun(edm::Run const&, edm::EventSetup const&);
   virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
   virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

   edm::EDGetToken tracks_;
   edm::EDGetToken hVtx_;
   edm::EDGetToken beamSpotHandle_;
   edm::EDGetToken trigger;


   TTree * Tree;
   TTree * Tree2;
    TTree * Tree3;
   TTree * treemc;


   

/*
   std::vector<unsigned int>  event_number;
   std::vector<unsigned int>  run_number;

   std::vector<double> PVx, PVy, PVz;

   std::vector<double> trackSlow_pt, trackSlow_eta, trackSlow_phi, trackSlow_charge;
   std::vector<double> track1_pt, track1_eta, track1_phi, track1_charge;
   std::vector<double> track2_pt, track2_eta, track2_phi, track2_charge;
   std::vector<double> track3_pt, track3_eta, track3_phi, track3_charge;
   std::vector<double> track4_pt, track4_eta, track4_phi, track4_charge;
   std::vector<double> D0_vtx_prob;
   std::vector<double> D0_mass;
   std::vector<double> D0_pt;
   std::vector<double> D0_eta;
   std::vector<double> D0_phi;
   std::vector<double> D_cosalpha;
   std::vector<double> D_cosalphaV;
   std::vector<double> D_Lsig;



*/

   //from D*

   int i;
   int event_number = 0, run_number = 0, lumi = 0;

   double Pis_pt, Pis_eta, Pis_phi;
   double Pi_pt, Pi_eta, Pi_phi;
   double K_pt, K_eta, K_phi;
   double DS_pt, DS_eta, DS_phi, DS_mass_reco;
   double D0_pt, D0_eta, D0_phi, D0_mass_reco;
   double CL_vertex, sum_ptVertex2, sum_ptVertex, Delta_mass;
   double x_p, y_p, z_p;
   double x_s, y_s, z_s;
   double L_abs, L_sigma;
   double cos_phi, pt_tracks;
   int Kiter, Piiter, Pisiter;
   int vertex_size, num_tracks;
   bool passTrigger;
   double trackSlow_pt, trackSlow_eta, trackSlow_phi, trackSlow_charge;
   double track1_pt, track1_eta, track1_phi, track1_charge;
   double track2_pt, track2_eta, track2_phi, track2_charge;
   double track3_pt, track3_eta, track3_phi, track3_charge;
   double track4_pt, track4_eta, track4_phi, track4_charge;

};


 testDstar_AOD:: testDstar_AOD(const edm::ParameterSet& iConfig)
{
   edm::Service<TFileService> fs;

   Tree = fs->make<TTree>("tree", "tree");
   Tree2 = fs->make<TTree>("tree2", "tree2");
   Tree3 = fs->make<TTree>("treemc", "treemc");


   Tree->Branch("event_number", &event_number);
   Tree->Branch("run_number", &run_number);
   Tree->Branch("lumi", &lumi, "lumi/I");

   Tree->Branch("Pis_pt",&Pis_pt,"Pis_pt/D");
   Tree->Branch("Pis_eta",&Pis_eta,"Pis_eta/D");
   Tree->Branch("Pis_phi",&Pis_phi,"Pis_phi/D");

   Tree->Branch("Pi_pt",&Pi_pt,"Pi_pt/D");
   Tree->Branch("Pi_eta",&Pi_eta,"Pi_eta/D");
   Tree->Branch("Pi_phi",&Pi_phi,"Pi_phi/D");

   Tree->Branch("K_pt",&K_pt,"K_pt/D");
   Tree->Branch("K_eta",&K_eta,"K_eta/D");
   Tree->Branch("K_phi",&K_phi,"K_phi/D");

   Tree->Branch("D0_pt",&D0_pt,"D0_pt/D");
   Tree->Branch("D0_eta",&D0_eta,"D0_eta/D");
   Tree->Branch("D0_phi",&D0_phi,"D0_phi/D");
   Tree->Branch("D0_mass_reco",&D0_mass_reco,"D0_mass_reco/D");

   Tree->Branch("DS_pt",&DS_pt,"DS_pt/D");
   Tree->Branch("DS_eta",&DS_eta,"DS_eta/D");
   Tree->Branch("DS_phi",&DS_phi,"DS_phi/D");
   Tree->Branch("DS_mass_reco",&DS_mass_reco,"DS_mass_reco/D");

   Tree->Branch("Delta_mass",&Delta_mass,"Delta_mass/D");

   Tree->Branch("Kiter", &Kiter, "Kiter/I");
   Tree->Branch("Piiter", &Piiter, "Piiter/I");
   Tree->Branch("Pisiter", &Pisiter, "Pisiter/I");

   Tree3->Branch("CL_vertex",&CL_vertex,"CL_vertex/D");

   Tree->Branch("L_abs",&L_abs,"L_abs/D");
   Tree->Branch("L_sigma",&L_sigma,"L_sigma/D");
   Tree->Branch("cos_phi",&cos_phi,"cos_phi/D");
   Tree->Branch("vertex_size",&vertex_size,"vertex_size/I");
   Tree->Branch("passTrigger", &passTrigger, "passTrigger/B"); 
   Tree->Branch("trackSlow_pt", &trackSlow_pt, "trackSlow_pt/D");

   Tree->Branch("trackSlow_eta", &trackSlow_eta, "trackSlow_eta/D");
   Tree->Branch("trackSlow_phi", &trackSlow_phi, "trackSlow_phi/D");
   Tree->Branch("trackSlow_charge", &trackSlow_charge, "trackSlow_charge/D"); 

   Tree->Branch("track1_pt", &track1_pt, "track1_pt/D");
   Tree->Branch("track1_eta", &track1_eta, "track1_eta/D");
   Tree->Branch("track1_phi", &track1_phi, "track1_phi/D");
   Tree->Branch("track1_charge", &track1_charge, "track1_charge/D");

   Tree->Branch("track2_pt", &track2_pt, "track2_pt/D");
   Tree->Branch("track2_eta", &track2_eta, "track2_eta/D");
   Tree->Branch("track2_phi", &track2_phi, "track2_phi/D");
   Tree->Branch("track2_charge", &track2_charge, "track2_charge/D");

   Tree2->Branch("track3_pt", &track3_pt, "track3_pt/D");
   Tree2->Branch("track3_eta", &track3_eta, "track3_eta/D");
   Tree2->Branch("track3_phi", &track3_phi, "track3_phi/D");
   Tree2->Branch("track3_charge", &track3_charge, "track3_charge/D");

   Tree->Branch("track4_pt", &track4_pt, "track4_pt/D");
   Tree->Branch("track4_eta", &track4_eta, "track4_eta/D");
   Tree->Branch("track4_phi", &track4_phi, "track4_phi/D");
   Tree->Branch("track4_charge", &track4_charge, "track4_charge/D");


   tracks_= consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
   hVtx_= consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
   beamSpotHandle_= consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
   trigger = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));

}

 testDstar_AOD::~testDstar_AOD()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}



void testDstar_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
 

  //new for AOD

 edm::Handle<reco::TrackCollection> tracks;
 iEvent.getByToken(tracks_ , tracks );

 reco::BeamSpot beamSpot;
 edm::Handle<reco::BeamSpot> beamSpotHandle;
 iEvent.getByToken(beamSpotHandle_ , beamSpotHandle);

 edm::ESHandle<TransientTrackBuilder> theB;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

 edm::ESHandle<TransientTrackBuilder> theB1;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB1);

  

 //form D*

  //new for AOD

  // *******
// TRIGGER
// *******

 edm::Handle<edm::TriggerResults> trigResults;
 iEvent.getByToken(trigger, trigResults);

  const edm::TriggerNames& names = iEvent.triggerNames(*trigResults);
 bool passTrig=false;


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


/* for (unsigned int i=0; i<names.size(); i++){
   if(names.triggerName(i) == "HLT_ZeroBias_v5"
  || names.triggerName(i) == "HLT_ZeroBias_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part0_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part1_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part2_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part3_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part4_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part5_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part6_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part7_v6"
  || names.triggerName(i) == "HLT_ZeroBias_part0_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part1_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part2_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part3_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part4_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part5_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part6_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part7_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part8_v2"
  || names.triggerName(i) == "HLT_ZeroBias_part9_v2" ){
   if (trigResults->accept(names.triggerIndex(names.triggerName(i)))){
    passTrig=true;
    }
   }
 }*/

 edm::Handle<std::vector<reco::Vertex>> hVtx;
 iEvent.getByToken(hVtx_ , hVtx); 
 reco::Vertex primVertex;
 if(hVtx->size() > 0){

   reco::TransientTrack  Ktrans;
   reco::TransientTrack  pitrans;
   std::vector<reco::TransientTrack> tks; 


   vector<reco::Track> myvector;


   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trigger, trigResults);

   reco::Track tr1;
   reco::Track tr2;

   reco::Track goodTrack1;
   reco::Track goodTrack2;


 
  TLorentzVector Kvector;
  TLorentzVector Pivector;
  TLorentzVector D0vector;
  TLorentzVector Pisvector;
  TLorentzVector DSvector;
  TransientVertex v;
  reco::Vertex primaryVertex;
  reco::Vertex primVertex_tmp;
  bool goodvertex = false;

  if(hVtx->size() > 0){
  double bsx = 0, bsy = 0, bsz = 0;
    if ( beamSpotHandle.isValid() ){ beamSpot = *beamSpotHandle;
    bsx = beamSpot.x0();
    bsy = beamSpot.y0();
    bsz = beamSpot.z0();
    }
    GlobalPoint BeamSpotGP(bsx, bsy, bsz);
    primaryVertex = hVtx->at(0);
    
    for(unsigned int t = 0; t<hVtx->size(); t++){
    primVertex_tmp = hVtx->at(t);
    if(!primVertex_tmp.isFake() && primVertex_tmp.isValid() && primVertex_tmp.ndof() > 4 && fabs(primVertex_tmp.z()-bsz) < 10 ){
      goodvertex = true;
      primaryVertex = primVertex_tmp;
      break;
    } 
    }


 

  

   
   double pion_mass = 139.57039/1000.;
   double kaon_mass = 493.677/1000.;
   double D0_mass = 1.864841;

   if(goodvertex){ //already checked
   double somma_ptVertex = 0;
      for (std::vector<TrackBaseRef >::const_iterator tracks = primaryVertex.tracks_begin(); tracks != primaryVertex.tracks_end(); ++tracks) {
    const reco::Track *track = tracks->get();
    somma_ptVertex += track->pt();
   }

   int sameVertex = 0, ntracks = 0, pttracks = 0;
   double inv_mass1 = 0;
   double inv_mass2 = 0, sigma_L = 0, CL = 0;
   bool goodtracks=false, goodDstarcandidate=false;
   double deriv[3];
   double cov_sv[3][3];
   double cov_pv[3][3];
  // int K_iteration = 0, Pi_iteration = 0, Pis_iteration = 0;
   double Lsusigma = 0, xprim = 0, yprim = 0, zprim = 0, xsec = 0, ysec = 0, zsec = 0, dx = 0, dy = 0, dz = 0, px = 0, py = 0, pz = 0, pi = 0, cosphi = 0, L = 0;
   double diff_mass = 100;
   //double diff_mass_min = 0.0;
   double diff_mass_max = 0.16;

    for (reco::TrackCollection::const_iterator track = tracks->begin();  track != tracks->end();  ++track){
    ntracks ++;
    pttracks = pttracks + track->pt();
   if((track->pt() > 0.7) && ( (track->chi2())/(track->ndof()) <= 2.5) && (track->numberOfValidHits() >= 5) && (track->hitPattern().numberOfValidPixelHits() >= 2) && fabs(track->dxy(primaryVertex.position())) < 0.1 && fabs(track->dz(primaryVertex.position())) < 1 && (track->quality(Track::highPurity))) {
    myvector.push_back(*track);
     }
   }
   if (myvector.size() > 1){
    for (reco::TrackCollection::const_iterator track = myvector.begin();  track != myvector.end();  ++track) {
     for (reco::TrackCollection::const_iterator track1 = track + 1; track1 != myvector.end(); ++track1){
      if (fabs(track->vz() - track1->vz()) < 0.8 && track->charge() + track1->charge() == 0){
       goodtracks=false;
       inv_mass1 = 0;
       tr1 = *track;
       tr2 = *track1;
       if ((track->charge()==1) && (track1->charge()==-1)){        
            //D0
            Kvector.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), kaon_mass);
            Pivector.SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), pion_mass);
            D0vector = Kvector + Pivector;
            inv_mass1 = D0vector.M(); 
            goodtracks = true;
           // K_iteration = track1->originalAlgo();
           // Pi_iteration = track->originalAlgo();
            goodtracks = true;
            goodTrack1 = *track;
            goodTrack2 = *track1;
           Ktrans = (*theB).build(tr2);
           pitrans = (*theB).build(tr1);
           }
          else if ((track->charge()==-1) && (track1->charge()==1)){
            //Dbar0
            Kvector.SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), kaon_mass);
            Pivector.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), pion_mass);
           // K_iteration = track->originalAlgo();
            //Pi_iteration = track1->originalAlgo();
            D0vector = Kvector + Pivector;
            inv_mass1 = D0vector.M();
            goodtracks = true;
            goodTrack1 = *track;
            goodTrack2 = *track1;
            Ktrans  = (*theB).build(tr1);
            pitrans = (*theB).build(tr2);
           }


    if (fabs(inv_mass1 - D0_mass) < 0.025){
            CL = 0;
            sameVertex = 0;   
            // ######### vertex fit ##########
            tks.clear();
            tks.push_back(Ktrans);
            tks.push_back(pitrans);
            if (tks.size() > 1){
             KalmanVertexFitter kalman(true);
             if(goodtracks){ v = kalman.vertex(tks);
              if(v.isValid()){sameVertex++;
               CL = 0;
               CL = TMath::Prob(v.totalChiSquared(),(int)v.degreesOfFreedom());
                CL_vertex = CL;
                Tree3->Fill();
              }  
             } 
            }
         
            if(sameVertex > 0 && CL > 0.01 && goodtracks == true){


             xprim = primaryVertex.position().x();
             yprim = primaryVertex.position().y();
             zprim = primaryVertex.position().z();

             xsec = v.position().x();
             ysec = v.position().y();
             zsec = v.position().z();
  
             dx = xsec - xprim;
             dy = ysec - yprim;
             dz = zsec - zprim;
             L_abs = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

             px = D0vector.Px();
             py = D0vector.Py();
             pz = D0vector.Pz();
             pi = D0vector.P();
             L = 0;
             cosphi = (px*dx + py*dy + pz*dz)/(L_abs*pi);
             cov_sv[0][0] = v.positionError().cxx();
             cov_sv[1][0] = v.positionError().cyx();
             cov_sv[2][0] = v.positionError().czx();
             cov_sv[0][1] = cov_sv[1][0];
             cov_sv[1][1] = v.positionError().cyy();
             cov_sv[2][1] = v.positionError().czy();
             cov_sv[0][2] = cov_sv[2][0];
             cov_sv[1][2] = cov_sv[2][1];
             cov_sv[2][2] = v.positionError().czz();

             cov_pv[0][0] = primaryVertex.covariance(0,0);
             cov_pv[1][0] = primaryVertex.covariance(1,0);
             cov_pv[2][0] = primaryVertex.covariance(2,0);
             cov_pv[0][1] = cov_pv[1][0];
             cov_pv[1][1] = primaryVertex.covariance(1,1);
             cov_pv[2][1] = primaryVertex.covariance(2,1);
             cov_pv[0][2] = cov_pv[2][0];
             cov_pv[1][2] = cov_pv[2][1];
             cov_pv[2][2] = primaryVertex.covariance(2,2);
             deriv[0] = dx/L_abs;
             deriv[1] = dy/L_abs;
             deriv[2] = dz/L_abs;
 
             if(cosphi > 0.99){
              L = (px*dx + py*dy + pz*dz)/pi;
              sigma_L = 0;
              for (int m = 0; m < 3; ++m){
               for (int n = 0; n < 3; ++n ){
                sigma_L += deriv[m]*deriv[n]*(cov_pv[m][n] + cov_sv[m][n]);
               }
              }
              sigma_L = sqrt(sigma_L);
              Lsusigma = L/sigma_L;
              if(Lsusigma > 3.){

               if(D0vector.Pt() > 3.){

             for(reco::TrackCollection::const_iterator track3= tracks->begin();  track3!= tracks->end();  ++track3){
               
              if((track3->pt() > 0.3) && (track3->charge()==1) && ( (track3->chi2())/(track3->ndof()) <= 3) && (track3->numberOfValidHits() >= 2) && fabs(track3->dxy(primaryVertex.position())/track3->dxyError()) < 3 && fabs(track3->dz(primaryVertex.position())/track3->dzError()) < 3){
                if(track3->pt() != Pivector.Pt() && track3->eta() != Pivector.Eta()){
                 
                    Pisvector.SetPtEtaPhiM(track3->pt(), track3->eta(), track3->phi(), pion_mass);
                    //Pis_iteration = track3->originalAlgo();
                    DSvector = Pisvector + Kvector + Pivector;
                    inv_mass2= DSvector.M();

         
                    if((DSvector.Pt() > 4.)){


                     if (fabs(DSvector.M() - D0vector.M()) < diff_mass_max){
                      diff_mass = fabs(DSvector.M() - D0vector.M());

                      Pis_pt=Pisvector.Pt();
                      Pis_eta=Pisvector.Eta();
                      Pis_phi=Pisvector.Phi();

                      Pi_pt=Pivector.Pt();
                      Pi_eta=Pivector.Eta();
                      Pi_phi=Pivector.Phi();

                      K_pt=Kvector.Pt();
                      K_eta=Kvector.Eta();
                      K_phi=Kvector.Phi();

                      D0_pt=D0vector.Pt();
                      D0_eta=D0vector.Eta();
                      D0_phi=D0vector.Phi();
                      //D0_mass_reco=D0vector.M();
                      D0_mass_reco=inv_mass1;
                      //std::cout<<"D0.M: "<<D0vector.M()<<" D0 inv_mass1: "<<inv_mass1<<std::endl;

                      DS_pt=DSvector.Pt();
                      DS_eta=DSvector.Eta();
                      DS_phi=DSvector.Phi();
                      DS_mass_reco= inv_mass2;

                      Delta_mass=diff_mass;
                     // std::cout<<""<<Delta_mass<<std::endl;
                      CL_vertex = CL;
                      vertex_size = hVtx->size();
                      cos_phi = (px*dx + py*dy + pz*dz)/(L_abs*pi);

                      x_p = xprim;
                      y_p = yprim;
                      z_p = zprim;

                      x_s = xsec;
                      y_s = ysec;
                      z_s = zsec;
                      L_sigma = L/sigma_L;
                      //event_number = iEvent.id().event();
                      run_number = iEvent.id().run();
                      lumi = iEvent.id().luminosityBlock();
                      sum_ptVertex = somma_ptVertex;
                      sum_ptVertex2 = pow(sum_ptVertex,2);   
                      passTrigger = passTrig;
                      goodDstarcandidate = true;
                   }
               }
              }//terza traccia != altre
             }//if terza traccia
            }//terza traccia
           }//pt D0
          }
         }//cos phi positivo
        }//CL SV
       }//mass
      }//total charge 0
     }//fortrack1
    }//for track
    if (diff_mass < 0.160 && goodDstarcandidate) Tree->Fill();
   }//myvector 
  }//vertex good
 }//hVtx>0 

}
}

void 
testDstar_AOD::beginJob()
{

edm::Service<TFileService> fs;

//variabili vertice

}

// ------------ method called once each job just after ending the event loop  ------------
void 
testDstar_AOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
testDstar_AOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
testDstar_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
testDstar_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
testDstar_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
testDstar_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE (testDstar_AOD);

