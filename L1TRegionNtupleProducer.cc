/*
 * \file L1TRegionNtupleProducer.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Run3Ntuplizer/interface/L1TRegionNtupleProducer.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"


#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>

using namespace edm;
using std::cout;
using std::endl;
using std::vector;



struct cluster{
  float pt;
  float eta;
  float phi;
  
  float topN_pt;
  float topN_eta;
  float topN_phi;
  //....


  //left neighbor
  float leftN_pt;
  float leftN_eta;
  float leftN_phi;

  float id;

  //tells us if the centra seed for this cluster is the maximum in the region
  float maxCluster;
};

// Tunable Distance Parameter
const int R = 1;
// Loop over all regions, no self-comparison, no duplicates
void L1TRegionNtupleProducer::antikt(){
  struct delta {
    float dist;
    float bist;
    float etai;
    float etaj;
    float phii;
    float phij;
    bool cluster;
  };

  int max =100;
    for (int i=0; i<max; i++) {
      for (int j=0; j!=i && j>i; j++) {
	// Region Distance Metric
	float dist = std::min(pow(vRegionEt[i],-2.0),pow(vRegionEt[j],-2.0))*sqrt(pow((vRegionEta[i]-vRegionEta[j]),2.0)+pow((vRegionPhi[i]-vRegionPhi[j]),2.0))/pow(R,2);
	// Beam Distance Metric
	float bist = std::min(pow(vRegionEt[i],-2.0),pow(vRegionEt[j],-2.0));
	// Set vector of parameters of each pair
	std::vector<delta> Delta;
	
	Delta[i].dist = dist;
	Delta[i].bist = bist;
	Delta[i].cluster = dist < bist;
	Delta[i].etai = vRegionEta[i];
	Delta[i].etaj = vRegionEta[j];
	Delta[i].phii = vRegionPhi[i];
	Delta[i].phij = vRegionPhi[j];
    }
  }
};
  
 
L1TRegionNtupleProducer::L1TRegionNtupleProducer( const ParameterSet & cfg ) :
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion")))
  {


    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    //recoPt_              = cfg.getParameter<double>("recoPtCut");
    folder               = tfs_->mkdir(folderName_);
    //folder->cd();
    regionTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree");
    regionTree->Branch("run",        &run,     "run/I");
    regionTree->Branch("lumi",       &lumi,    "lumi/I");
    regionTree->Branch("event",      &event,   "event/I");
    regionTree->Branch("nvtx",       &nvtx,    "nvtx/D");
    regionTree->Branch("vRegionEt",  &vRegionEt  );
    regionTree->Branch("vRegionEta", &vRegionEta );
    regionTree->Branch("vRegionPhi", &vRegionPhi );
    regionTree->Branch("vRegionEG",  &vRegionEG  );
    regionTree->Branch("vRegionTau", &vRegionTau );

    nEvents       = folder.make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );

    regionHitEta  = folder.make<TH1F>( "regionHit_eta"  , "eta", 16, 1, 16. );
    regionHitPhi  = folder.make<TH1F>( "regionHit_phi"  , "phi", 16, 1, 16. );
    regionTotal   = folder.make<TH1F>( "regionHit_total"  , "fullmap", 16, 1, 16. );

    regionEta     = folder.make<TH1F>( "region_eta"  , "eta", 22, 1, 22. );
    regionPhi     = folder.make<TH1F>( "region_phi"  , "phi", 72, 1, 72. );
    regionPt      = folder.make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );

    regionEtaFine   = folder.make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
    regionPhiFine   = folder.make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );

  }

void L1TRegionNtupleProducer::beginJob( const EventSetup & es) {
   std::cout<<"begin job..."<<std::endl;
}

void L1TRegionNtupleProducer::analyze( const Event& evt, const EventSetup& es )
 {
   std::cout<<"Analyzing..."<<std::endl;
   nEvents->Fill(1);
   
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   Handle<L1CaloRegionCollection> regions;
   

   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;
   
   
  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  std::cout<<"making regions"<<std::endl;
  vRegionEt.clear();
  vRegionEta.clear();
  vRegionPhi.clear();
  vRegionTau.clear();
  vRegionEG.clear();

  vector<cluster> clusters;

  UCTGeometry g;
  //************* Get Regions and make region plots


  int i = 0;
  if(!evt.getByToken(regionSource_,regions)){
    std::cout<<"ERROR GETTING THE REGIONS!!!"<<std::endl;}
  else{
    for(vector<L1CaloRegion>::const_iterator testRegion = regions->begin(); testRegion != regions->end(); ++testRegion){
      //UCTRegionProcess uctRegion(*region);
       uint16_t test_raw = testRegion->raw();
       uint32_t test_et = testRegion->et();
       //testRegionTotET += test_et;
       uint32_t test_rEta = testRegion->id().ieta();
       uint32_t test_rPhi = testRegion->id().iphi();
       UCTRegionIndex test_rIndex = g.getUCTRegionIndexFromL1CaloRegion(test_rEta, test_rPhi);
       UCTTowerIndex test_tIndex = g.getUCTTowerIndexFromL1CaloRegion(test_rIndex, test_raw);
       int test_cEta = test_tIndex.first;
       int test_cPhi = test_tIndex.second;

       if(testRegion->et()>0 && fabs(test_cEta)<28 && test_cPhi < 72){
	 float pt = test_et;
	 float eta = 0;
	 if(fabs(test_cEta)<28){
	   eta = towerEtaMap[(int)fabs(test_cEta)];
	   eta = eta*fabs(test_cEta)/test_cEta;
	 }

	 float phi = towerPhiMap[test_cPhi];
	 float isEgammaLike = 0;
	 float isTauLike    = 0;
	 if(!((l1tcalo::RegionEGVeto & testRegion->raw()) == l1tcalo::RegionEGVeto))
	   isEgammaLike = 1;

	 if(!((l1tcalo::RegionTauVeto & testRegion->raw()) == l1tcalo::RegionTauVeto))
	   isTauLike = 1;
	 //std::cout<<"region eta,phi: "<<eta<<" , "<<phi<<std::endl;
	 vRegionEt.push_back(pt);
	 vRegionEta.push_back(eta);
	 vRegionPhi.push_back(phi);
	 vRegionEG.push_back(isEgammaLike);
	 vRegionEG.push_back(isTauLike);
	 

	 // create cluster seeds
	 cluster temp;
	 temp.pt = pt;
	 temp.eta = eta;
	 temp.phi = phi;
	 temp.id = i;
	 temp.maxCluster = 1;

	 regionEta->Fill(eta);
	 regionPhi->Fill(phi);
	 regionPt->Fill(pt);

	 clusters.push_back(temp);
	 i++;
      }
    }
    regionTree->Fill();
  }
  
  //Each region is 0.087 x 4 in eta and 0.087 x 4 in phi
  for(auto seed : clusters){
    for(auto cluster_iter : clusters){
      /// test different sizes for distance parameter in the future
      if((fabs(seed.eta - cluster_iter.eta)< 0.087*4)&&fabs(seed.phi-cluster_iter.phi)<0.0874){
	///this is a neighboring cluster
	if((seed.eta - cluster_iter.eta) > 0){
	  seed.leftN_pt = cluster_iter.pt;
	  seed.leftN_pt = cluster_iter.eta;
	  seed.leftN_pt = cluster_iter.phi;

	  if(cluster_iter.pt > seed.pt)
	    cluster_iter.maxCluster = 0;
	}
	// do this for all neighbors
	
      }
    }
  }
  



}


  
/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */
  
void L1TRegionNtupleProducer::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  
  //std::cout << "ECAL TPGS" << std::endl;
  for (size_t i = 0; i < ecal->size(); ++i) {
    int cal_ieta = (*ecal)[i].id().ieta();
    int cal_iphi = (*ecal)[i].id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecal)[i].compressedEt()*LSB;
    //if(et>0)std::cout<<"et "<< et<<std::endl;
    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    //if(et>0)
    //cout<<"Before filling eTower"
    //<<"ECAL ieta:"<<ieta<<" cal_ieta:"<< cal_ieta<<" iphi:"<<iphi<<" et:"<<et<<endl;

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

 void L1TRegionNtupleProducer::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
					 const ESHandle<L1CaloHcalScale>  hcalScale, 
					 double hTowerETMap[73][57], bool testMode){
  for (size_t i = 0; i < hcal->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcal)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double energy = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(energy>0)std::cout<<"energy "<< energy<<std::endl;

    if(testMode && iphi == 34 && ieta == 12){
      energy = 40;
    }

    if (iphi >= 0 && iphi <= 71 &&
	ieta >= 0 && ieta <= 55) {

      //(*hcal)[i].SOI_compressedEt(), absieta, zside)*LSB; //*LSB
      //if(energy>0)
      //std::cout<<"hcal iphi "<<iphi<<" ieta "<<ieta<<" energy "<<energy<<std::endl;
      hTowerETMap[iphi][ieta] = energy;
      //TPGSum_ +=energy;
      //TPGH_ += energy;
      //double alpha_h = TPGSFp_[cal_ieta]; //v3
      //hCorrTowerETMap[cal_iphi][cal_ieta] = alpha_h*energy;
      //cTPGH_ += alpha_h*energy;
      //if (energy > 0) {
      //std::cout << "hcal eta/phi=" << ieta << "/" << iphi
      //<< " = (" << getEtaTPG(ieta) << "/" << getPhiTPG(iphi) << ") "
      //<< " et=" << (*hcal)[i].SOI_compressedEt()
      //<< " energy=" << energy
      //<< " rctEta="<< twrEta2RegionEta(cal_ieta) << " rctPhi=" << twrPhi2RegionPhi(cal_iphi)
      //<< " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
      //}
      //if (energy>maxTPGHPt){
      //maxTPGHPt=energy;
      //maxTPGHPt_phi = cal_iphi; //this one starts at 0-72
      //maxTPGHPt_eta = cal_ieta; //this one is 0-54
      //} 
    }
    //else
      //std::cout<<"HCAL failed checks iphi "<<iphi<<" ieta "<<ieta<<std::endl;
  }//end HCAL TPG
}

void L1TRegionNtupleProducer::endJob() {
}

L1TRegionNtupleProducer::~L1TRegionNtupleProducer(){

}

DEFINE_FWK_MODULE(L1TRegionNtupleProducer);
