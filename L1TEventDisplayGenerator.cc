/*
 *  \file L1TEventDisplayGenerator.cc
 * 
 *  \author I. Ojalvo
 *  Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Run3Ntuplizer/interface/L1TEventDisplayGenerator.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

L1TEventDisplayGenerator::L1TEventDisplayGenerator( const ParameterSet & cfg ) :
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion"))),
  ecalCaloSrc_(consumes<vector <reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("ecalCaloClusters")))
  {
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Efficiency Tree");

    //efficiencyTree->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    //efficiencyTree->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    //efficiencyTree->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    //efficiencyTree->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    //efficiencyTree->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    //efficiencyTree->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    //efficiencyTree->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    //efficiencyTree->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    //efficiencyTree->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    //efficiencyTree->Branch("rlxTaus", "vector<TLorentzVector>", &rlxTaus, 32000, 0); 
    //efficiencyTree->Branch("isoTaus", "vector<TLorentzVector>", &isoTaus, 32000, 0); 
    //efficiencyTree->Branch("recoTaus", "vector<TLorentzVector>", &recoTaus, 32000, 0); 
    efficiencyTree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0); 
    efficiencyTree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    //efficiencyTree->Branch("signalPFCands", "vector<TLorentzVector>", &signalPFCands, 32000, 0); 
    //efficiencyTree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0); 
    //efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0); 
    //efficiencyTree->Branch("recoJetsDR", "vector<double>", &recoJetsDr, 32000, 0); 
    efficiencyTree->Branch("caloClusters", "vector<TLorentzVector>", &caloClusters, 32000, 0); 

    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    //efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
    //
    //efficiencyTree->Branch("tauEtaEcalEnt", &tauEtaEcalEnt,"tauEtaEcalEnt/D");
    //efficiencyTree->Branch("tauPhiEcalEnt", &tauPhiEcalEnt,"tauPhiEcalEnt/D");

    //efficiencyTree->Branch("recoPt",        &recoPt,   "recoPt/D");
    //efficiencyTree->Branch("isoTauPt",      &isoTauPt, "isoTauPt/D");
    //efficiencyTree->Branch("rlxTauPt",      &rlxTauPt, "rlxTauPt/D");
    //
    //efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
    //efficiencyTree->Branch("isoTauEta",     &isoTauEta, "isoTauEta/D");
    //efficiencyTree->Branch("rlxTauEta",     &rlxTauEta, "rlxTauEta/D");
    //
    //efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");
    //efficiencyTree->Branch("isoTauPhi",     &isoTauPhi, "isoTauPhi/D");
    //efficiencyTree->Branch("rlxTauPhi",     &rlxTauPhi, "rlxTauPhi/D");

    //efficiencyTree->Branch("l1IsoMatched",  &l1IsoMatched, "l1IsoMatched/I");
    //efficiencyTree->Branch("l1RlxMatched",  &l1RlxMatched, "l1RlxMatched/I");
    
  }

void L1TEventDisplayGenerator::beginJob( const EventSetup & es) {
}

void L1TEventDisplayGenerator::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle<reco::VertexCollection> vertices;
  if(evt.getByToken(vtxLabel_, vertices)){
    nvtx = (int) vertices->size();
    std::cout<<"nVertices "<<nvtx<<std::endl;
  }
  
  Handle<L1CaloRegionCollection> regions;
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  edm::Handle < vector<reco::CaloCluster> > recoCaloClusters;

  allRegions->clear(); 
  allEcalTPGs->clear(); 
  allHcalTPGs->clear(); 
  caloClusters->clear(); 

  if(evt.getByToken(ecalCaloSrc_, recoCaloClusters)){
    for( vector<reco::CaloCluster>::const_iterator caloCluster = recoCaloClusters->begin(); 
	 caloCluster != recoCaloClusters->end(); 
	 caloCluster++ ) {
      //fill vector
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(caloCluster->energy(),caloCluster->eta(),caloCluster->phi(),caloCluster->energy());
      caloClusters->push_back(temp);
    }
  }

  UCTGeometry g;
  if(!evt.getByToken(regionSource_,regions)){
    std::cout<<"ERROR GETTING THE REGIONS!!!"<<std::endl;}
  else{
    for(vector<L1CaloRegion>::const_iterator testRegion = regions->begin(); testRegion != regions->end(); ++testRegion){
      //UCTRegionProcess uctRegion(*region);
      uint16_t test_raw = testRegion->raw();
      uint32_t test_et = testRegion->et();
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
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(pt,eta,phi,pt);
      allRegions->push_back(temp);
      }
    }
  }

  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < ecalTPGs->size(); ++i) {
      int cal_ieta = (*ecalTPGs)[i].id().ieta();
      int cal_iphi = (*ecalTPGs)[i].id().iphi();
      if(cal_iphi==0)
	std::cout<<"cal_phi is 0"<<std::endl;
      if(cal_ieta<-28)
	continue;
      if(cal_ieta>28)
	continue;
      int ieta = TPGEtaRange(cal_ieta);
      short zside = (*ecalTPGs)[i].id().zside();
      // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
      // TPG ieta ideal goes from 0-55.
      double LSB = 0.5;
      double et= (*ecalTPGs)[i].compressedEt()*LSB;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);
      //if(et>0)
      //std::cout<<"et "<<et<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      //if(et>5)
      //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
      allEcalTPGs->push_back(temp);
    }

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }

  efficiencyTree->Fill();
 }

int L1TEventDisplayGenerator::get5x5TPGs(const int maxTPGPt_eta, 
				     const int maxTPGPt_phi, 
				     const double eTowerETMap[73][57], 
				     const double hTowerETMap[73][57], 
				     std::vector<double>* hcalTpgs_pt, 
				     std::vector<double>* hcalTpgs_eta, 
				     std::vector<double>* hcalTpgs_phi, 
				     std::vector<double>* ecalTpgs_pt, 
				     std::vector<double>* ecalTpgs_eta, 
				     std::vector<double>* ecalTpgs_phi,
				     std::vector<double>* sumTpgs_pt, 
				     std::vector<double>* sumTpgs_eta, 
				     std::vector<double>* sumTpgs_phi){
  for (int j = -5; j < 6; ++j) {//phi
    for (int k = -5; k < 6; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      hcalTpgs_pt->push_back(hTowerETMap[tpgsquarephi][tpgsquareeta]);
      hcalTpgs_eta->push_back(towerEtaMap[k]);
      hcalTpgs_phi->push_back(towerPhiMap[j]);

      ecalTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]);
      ecalTpgs_eta->push_back(towerEtaMap[tpgsquareeta]);
      ecalTpgs_phi->push_back(towerPhiMap[tpgsquarephi]);

      sumTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]+eTowerETMap[tpgsquarephi][tpgsquareeta]);
      sumTpgs_eta->push_back(towerEtaMap[tpgsquareeta]);
      sumTpgs_phi->push_back(towerPhiMap[tpgsquarephi]);
    }
  }
  //std::cout<<"TPGe5x5_ "<<TPGe5x5_<<" TPGh5x5_ "<<TPGh5x5_<<std::endl;
  //return (TPGe5x5_ + TPGh5x5_);
  return 1;
}

/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */

void L1TEventDisplayGenerator::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
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

    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

void L1TEventDisplayGenerator::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
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

void L1TEventDisplayGenerator::endJob() {
}

L1TEventDisplayGenerator::~L1TEventDisplayGenerator(){
}

DEFINE_FWK_MODULE(L1TEventDisplayGenerator);
