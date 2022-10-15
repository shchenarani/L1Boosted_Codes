// system include files
#include <memory>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTTower.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

#include "L1Trigger/L1TCaloLayer1/src/UCTObject.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometryExtended.hh"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace l1extra;
using namespace std;
using namespace l1tcalo;

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

float getRecoEtaNew(int caloEta){
  float eta = -999.;
  static bool first = true;
  static double twrEtaValues[42];
  if(first) {
    twrEtaValues[0] = 0;
    for(unsigned int i = 0; i < 20; i++) {
      twrEtaValues[i + 1] = 0.0436 + i * 0.0872;
    }
    twrEtaValues[21] = 1.785;
    twrEtaValues[22] = 1.880;
    twrEtaValues[23] = 1.9865;
    twrEtaValues[24] = 2.1075;
    twrEtaValues[25] = 2.247;
    twrEtaValues[26] = 2.411;
    twrEtaValues[27] = 2.575;
    twrEtaValues[28] = 2.825;
    twrEtaValues[29] = 999.;
    twrEtaValues[30] = (3.15+2.98)/2.;
    twrEtaValues[31] = (3.33+3.15)/2.;
    twrEtaValues[32] = (3.50+3.33)/2.;
    twrEtaValues[33] = (3.68+3.50)/2.;
    twrEtaValues[34] = (3.68+3.85)/2.;
    twrEtaValues[35] = (3.85+4.03)/2.;
    twrEtaValues[36] = (4.03+4.20)/2.;
    twrEtaValues[37] = (4.20+4.38)/2.;
    twrEtaValues[38] = (4.74+4.38*3)/4.;
    twrEtaValues[39] = (4.38+4.74*3)/4.;
    twrEtaValues[40] = (5.21+4.74*3)/4.;
    twrEtaValues[41] = (4.74+5.21*3)/4.;
    first = false;
  }
  uint32_t absCaloEta = abs(caloEta);
  if(absCaloEta <= 41) {
    if(caloEta < 0)
      eta =  -twrEtaValues[absCaloEta];
    else
      eta = +twrEtaValues[absCaloEta];
  }
  return eta;
};

float getRecoPhiNew(int caloPhi){
  float phi = -999.;
  if(caloPhi > 72) phi = +999.;
  uint32_t absCaloPhi = std::abs(caloPhi) - 1;
  if(absCaloPhi < 36)
    phi = (((double) absCaloPhi + 0.5) * 0.0872);
  else
    phi = (-(71.5 - (double) absCaloPhi) * 0.0872);
  return phi;
};


// class declaration
class BoostedJetStudies : public edm::EDAnalyzer {
public:
  explicit BoostedJetStudies(const edm::ParameterSet&);
  ~BoostedJetStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void zeroOutAllVariables();

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
  virtual void beginJob() override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<reco::CaloJet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;
  edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
  edm::EDGetTokenT<l1t::JetBxCollection> stage2JetToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> stage2EGToken_;
  edm::EDGetTokenT<l1t::TauBxCollection> stage2TauToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> stage2EtSumToken_;
  edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1BoostedToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPToken_;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPToken_;
  edm::EDGetTokenT<vector<pat::Tau> > tauSrc_;
  TH1F* nEvents;
  TH1F* Hist_tauVeto;
  TH1F* Hist_tauVetoEt25;
  TH1F* Hist_tauVetoEt50;
  TH1F* Hist_tauVetoEt100;
  TH1F* Hist_pt;
  TH1F* Hist_eta;
  TH1F* Hist_phi;
  TH1F* Hist_dR;
  TH1F* Hist_regionEta;
  TH1F* Hist_regionPhi;
  TH1F* Hist_pt_higgs;
TH1F* Hist_pt_res;
TH1F* boostedJetpt;

int nevent0=0;
int nevent1=0;
int nevent2=0;
int nevent3=0;
int nevent4=0;
int nevent5=0;
int nevent6=0;

int nevent10=0;
int nevent11=0;
int nevent12=0;
int nevent13=0;
int nevent14=0;
int nevent15=0;
int nevent16=0;

int nevent20=0;
int nevent21=0;
int nevent22=0;
int nevent23=0;
int nevent24=0;
int nevent25=0;
int nevent26=0;

int nevent30=0;
int nevent31=0;
int nevent32=0;
int nevent33=0;
int nevent34=0;
int nevent35=0;
int nevent36=0;

int nevent40=0;
int nevent41=0;
int nevent42=0;
int nevent43=0;
int nevent44=0;
int nevent45=0;
int nevent46=0;

int nevent50=0;
int nevent51=0;
int nevent52=0;
int nevent53=0;
int nevent54=0;
int nevent55=0;
int nevent56=0;

int nevent60=0;
int nevent61=0;
int nevent62=0;
int nevent63=0;
int nevent64=0;
int nevent65=0;
int nevent66=0;

int nevent000=0;


int tauveto_0=0;
int tauveto_1=0;
int tauveto_2=0;
int tauveto_3=0;
int tauveto_4=0;
int tauveto_5=0;
int tauveto_6=0;
int tauveto_7=0;
int tauveto_8=0;

  int run, lumi, event, nPV;

  double genPt_1, genEta_1, genPhi_1, genM_1, genDR;
  int genId, genMother;
  double recoPt_1, recoEta_1, recoPhi_1;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  double seedPt_1, seedEta_1, seedPhi_1;

  int l1NthJet_1;
  int recoNthJet_1;
  int seedNthJet_1;

  std::vector<string>  regionEta, regionPhi;
  double recoPt_;
  std::vector<int> nSubJets, nBHadrons, HFlav;
  std::vector<std::vector<int>> subJetHFlav;
  std::vector<float> tau1, tau2, tau3;

  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *seed180  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *egseed = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *tauseed  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *htseed = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ak8Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *subJets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ecalTPGs  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *hcalTPGs  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *boostedTaus = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *L1boostedTaus = new std::vector<TLorentzVector>;

  bool passL1SingleJet180, passL1HTT360;

  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  edm::Service<TFileService> tfs_;
  uint32_t nPumBins;

  std::vector< std::vector< std::vector < uint32_t > > > pumLUT;

  double caloScaleFactor;

  uint32_t jetSeed;
  uint32_t tauSeed;
  float tauIsolationFactor;
  uint32_t eGammaSeed;
  double eGammaIsolationFactor;
  double boostedJetPtFactor;

  bool verbose;
  int fwVersion;

  edm::EDGetTokenT<L1CaloRegionCollection> regionToken;

  UCTLayer1 *layer1;
  UCTSummaryCard *summaryCard;
};

BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  jetSrc_(    consumes<vector<reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  genSrc_( consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticles"))),
  stage2JetToken_(consumes<l1t::JetBxCollection>( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  stage2EGToken_(consumes<l1t::EGammaBxCollection>( edm::InputTag("caloStage2Digis","EGamma","RECO"))),
  stage2TauToken_(consumes<l1t::TauBxCollection>( edm::InputTag("caloStage2Digis","Tau","RECO"))),
  stage2EtSumToken_(consumes<l1t::EtSumBxCollection>( edm::InputTag("caloStage2Digis","EtSum","RECO"))),
  l1BoostedToken_(consumes<vector<l1extra::L1JetParticle>>( edm::InputTag("uct2016EmulatorDigis","Boosted",""))),
  vtxToken_(consumes<reco::VertexCollection>( edm::InputTag("offlineSlimmedPrimaryVertices"))),
  trigobjectsMINIAODToken_(consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("slimmedPatTrigger"))),
  trgresultsToken_(consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT"))),
  ecalTPToken_(consumes<EcalTrigPrimDigiCollection>( edm::InputTag("l1tCaloLayer1Digis"))),
  hcalTPToken_(consumes<HcalTrigPrimDigiCollection>( edm::InputTag("l1tCaloLayer1Digis"))),
  tauSrc_(consumes<vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("boostedTaus"))),

  nPumBins(iConfig.getParameter<unsigned int>("nPumBins")),
  pumLUT(nPumBins, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(13))),
  caloScaleFactor(iConfig.getParameter<double>("caloScaleFactor")),
  jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  tauIsolationFactor(iConfig.getParameter<double>("tauIsolationFactor")),
  eGammaSeed(iConfig.getParameter<unsigned int>("eGammaSeed")),
  eGammaIsolationFactor(iConfig.getParameter<double>("eGammaIsolationFactor")),
  boostedJetPtFactor(iConfig.getParameter<double>("boostedJetPtFactor")),
  verbose(iConfig.getParameter<bool>("verbose")),
  fwVersion(iConfig.getParameter<int>("firmwareVersion")),
  regionToken(consumes<L1CaloRegionCollection>( edm::InputTag("simCaloStage2Layer1Digis") ))

{
  // Initialize the Tree
  std::vector<double> pumLUTData;
  char pumLUTString[10];
  for(uint32_t pumBin = 0; pumBin < nPumBins; pumBin++) {
    for(uint32_t side = 0; side < 2; side++) {
      if(side == 0) sprintf(pumLUTString, "pumLUT%2.2dp", pumBin);
      else sprintf(pumLUTString, "pumLUT%2.2dn", pumBin);
      pumLUTData = iConfig.getParameter<std::vector < double > >(pumLUTString);
      for(uint32_t iEta = 0; iEta < std::max((uint32_t) pumLUTData.size(), MaxUCTRegionsEta); iEta++) {
  pumLUT[pumBin][side][iEta] = (uint32_t) round(pumLUTData[iEta] / caloScaleFactor);
      }
      if(pumLUTData.size() != (MaxUCTRegionsEta))
  edm::LogError("L1TCaloSummary") << "PUM LUT Data size integrity check failed; Expected size = " << MaxUCTRegionsEta
      << "; Provided size = " << pumLUTData.size()
      << "; Will use what is provided :(" << std::endl;
  }
   }

  summaryCard = new UCTSummaryCard(&pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);
  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  Hist_tauVeto = tfs_->make<TH1F>("Hist_tauVeto","Hist_tauVeto",10,0,10);
  Hist_tauVetoEt25 = tfs_->make<TH1F>("Hist_tauVetoEt25","Hist_tauVetoEt25",10,0,10);
  Hist_tauVetoEt50 = tfs_->make<TH1F>("Hist_tauVetoEt50","Hist_tauVetoEt50",10,0,10);
  Hist_tauVetoEt100 = tfs_->make<TH1F>("Hist_tauVetoEt100","Hist_tauVetoEt100",10,0,10);
  Hist_pt = tfs_->make<TH1F>("Hist_pt","Hist_pt",20,0,1000);
  Hist_eta = tfs_->make<TH1F>("Hist_eta","Hist_eta",30,-3.0,3.0);
  Hist_phi = tfs_->make<TH1F>("Hist_phi","Hist_phi",40,-4.0,4.0);
  Hist_dR = tfs_->make<TH1F>("Hist_dR","Hist_dR",40,0,2.0);
  Hist_regionEta = tfs_->make<TH1F>("Hist_regionEta","Hist_regionEta",10,0,10);
  Hist_regionPhi = tfs_->make<TH1F>("Hist_regionPhi","Hist_regionPhi",10,0,10);
  Hist_pt_higgs = tfs_->make<TH1F>("Hist_Higgs_pt","Higgs_Hist_pt",20,0,1000);
  boostedJetpt= tfs_->make<TH1F>("boostedJetpt","boostedJetpt",20,0,1000);
  Hist_pt_res=tfs_->make<TH1F>("Hist_pt_res","Hist_pt_res",20,-1,1);

  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Gen Matched Jet Tree ");
  createBranches(efficiencyTree);
}

BoostedJetStudies::~BoostedJetStudies() {
  if(summaryCard != 0) delete summaryCard;
}

//
// member functions
// ------------ method called to produce the data  ------------



void BoostedJetStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  using namespace edm;
  nEvents->Fill(1);
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  nPV = 0;
  edm::Handle<reco::VertexCollection> vtxHandle;
  if(!evt.getByToken(vtxToken_, vtxHandle)) cout<<"ERROR GETTING PRIMARY VERTICES"<<endl;
  evt.getByToken(vtxToken_, vtxHandle);
  nPV = vtxHandle->size();

  std::vector<reco::CaloJet> goodJets;
  std::vector<pat::Jet> goodJetsAK8;
  std::vector<l1t::Jet> seeds;
  std::vector<pat::Tau> goodTaus;

  l1Jets->clear();
  seed180->clear();
  egseed->clear();
  tauseed->clear();
  htseed->clear();
  ak8Jets->clear();
  subJets->clear();
  ecalTPGs->clear();
  hcalTPGs->clear();
  nSubJets.clear();
  nBHadrons.clear();
  subJetHFlav.clear();
  tau1.clear();
  tau2.clear();
  tau3.clear();
  boostedTaus->clear();
  regionEta.clear();
  regionPhi.clear();
  L1boostedTaus->clear();

  std::unique_ptr<L1JetParticleCollection> TauCands(new L1JetParticleCollection);

  UCTGeometry g;

  // Here we read region data from the region collection created by L1TCaloLayer1 instead of
  // // independently creating regions from TPGs for processing by the summary card. This results
  // // in a single region vector of size 252 whereas from independent creation we had 3*6 vectors
  // // of size 7*2. Indices are mapped in UCTSummaryCard accordingly.

  summaryCard->clearRegions();
    std::vector<UCTRegion*> inputRegions;
    inputRegions.clear();
    edm::Handle<std::vector<L1CaloRegion>> regionCollection;
    if(!evt.getByToken(regionToken, regionCollection)) edm::LogError("L1TCaloSummary") << "UCT: Failed to get regions from region collection!" ;
    evt.getByToken(regionToken, regionCollection);
    for (const L1CaloRegion &i : *regionCollection) {
      UCTRegionIndex r = g.getUCTRegionIndexFromL1CaloRegion(i.gctEta(), i.gctPhi());
      UCTTowerIndex t = g.getUCTTowerIndexFromL1CaloRegion(r, i.raw());
      uint32_t absCaloEta = std::abs(t.first);
      uint32_t absCaloPhi = std::abs(t.second);
      bool negativeEta = false;
      if (t.first < 0)
        negativeEta = true;
      uint32_t crate = g.getCrate(t.first, t.second);
      uint32_t card = g.getCard(t.first, t.second);
      uint32_t region = g.getRegion(absCaloEta, absCaloPhi);
      UCTRegion* test = new UCTRegion(crate, card, negativeEta, region, fwVersion);
      test->setRegionSummary(i.raw());
      inputRegions.push_back(test);
    }
    summaryCard->setRegionData(inputRegions);

    if(!summaryCard->process()) {
      edm::LogError("L1TCaloSummary") << "UCT: Failed to process summary card" << std::endl;
      exit(1);
    }

double pt = -999;
double eta = -999.;
double phi = -999.;
double mass= 0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////  Started by Abdollah
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//// We start the code from here
double pt_matched=-100;
double eta_matched=-100;
double phi_matched=-100;
//float min=100;
float genH_pt = -100;
float min_Pt = 1000;
float dR_matched=-100;

std::list<UCTObject*> boostedJetObjs = summaryCard->getBoostedJetObjs();
edm::Handle<reco::GenParticleCollection> genParticles;
if(evt.getByToken(genSrc_, genParticles)){
    // loop over gen particle and select only Higgs boson
for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
    if   (genparticle->pdgId()!=25) continue;
    genH_pt=genparticle->pt();
    Hist_pt_higgs->Fill(genH_pt);

    for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor * boostedJetPtFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    boostedJetpt->Fill(pt);

    float  DeltaR = deltaR(eta,phi ,genparticle->eta(), genparticle->phi());
    if (DeltaR > 1.0) continue;
    if (fabs(pt - genH_pt)< min_Pt){
    min_Pt= fabs(pt - genH_pt);
    eta_matched=eta;
    phi_matched=phi;
    pt_matched=pt;
    dR_matched=DeltaR;
//  cout<<"min_Pt = " <<min_Pt <<"  L1Pt="<<pt <<"  genPt=" <<genH_pt<< "  DeltaR="<<DeltaR<< "  genparticle->status="<<genparticle->status() << "\n";
            }
        }
     break;
//                cout<<"genparticle->pdgId() = " <<genparticle->pdgId() <<"  genparticle->status()="<<genparticle->status() << "  genparticle->pt()=" <<genparticle->pt()<<"  genparticle->eta()=" <<genparticle->eta()<< "  genparticle->phi()=" <<genparticle->phi()<< "\n";
    }
}



//Hist_pt_higgs->Fill(genH_pt);
Hist_pt->Fill(pt_matched);
Hist_eta->Fill(eta_matched);
Hist_phi->Fill(phi_matched);
Hist_dR->Fill(dR_matched);
Hist_pt_res->Fill((pt_matched-genH_pt)/genH_pt);

//std::list<UCTObject*> boostedJetObjs = summaryCard->getBoostedJetObjs();
for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
const UCTObject* object = *i;
pt = ((double) object->et()) * caloScaleFactor * boostedJetPtFactor;
eta = g.getUCTTowerEta(object->iEta());
phi = g.getUCTTowerPhi(object->iPhi());


    if (  pt == pt_matched &&  eta ==eta_matched &&  phi == phi_matched) {
    // loop over all 9 regions and check how many regions pass TauVeto
        int TauVetoIndex=0;
        int TauVetoIndexEt25=0;
        int TauVetoIndexEt50=0;
        int TauVetoIndexEt100=0;
        for(uint32_t iEta = 0; iEta < 3; iEta++){
            for(uint32_t iPhi = 0; iPhi < 3; iPhi++){
if (object->boostedJetRegionTauVeto()[3*iEta+iPhi] == 1 ) TauVetoIndex++;
if ( (object->et() * caloScaleFactor * boostedJetPtFactor) > 25 && object->boostedJetRegionTauVeto()[3*iEta+iPhi] == 1 ) TauVetoIndexEt25++;
if ( (object->et() * caloScaleFactor * boostedJetPtFactor) > 50 && object->boostedJetRegionTauVeto()[3*iEta+iPhi] == 1 ) TauVetoIndexEt50++;
if ( (object->et() * caloScaleFactor * boostedJetPtFactor) > 100 && object->boostedJetRegionTauVeto()[3*iEta+iPhi] == 1 ) TauVetoIndexEt100++;
            }
        }
std::cout<<"caloScaleFactor * boostedJetPtFactor "<<caloScaleFactor <<"   " <<boostedJetPtFactor<<"\n";
Hist_tauVeto->Fill(TauVetoIndex);
Hist_tauVetoEt25->Fill(TauVetoIndexEt25);
Hist_tauVetoEt50->Fill(TauVetoIndexEt50);
Hist_tauVetoEt100->Fill(TauVetoIndexEt100);

bitset<3> activeRegionEtaPattern = 0;
for(uint32_t iEta = 0; iEta < 3; iEta++){
bool activeStrip = false;
for(uint32_t iPhi = 0; iPhi < 3; iPhi++){
if(object->boostedJetRegionET()[3*iEta+iPhi] > 30 && object->boostedJetRegionET()[3*iEta+iPhi] > object->et()*0.0625  && object->boostedJetRegionTauVeto()[3*iEta+iPhi] ==0) activeStrip = true;

            }
            if(activeStrip) activeRegionEtaPattern |= (0x1 << iEta);
        }

bitset<3> activeRegionPhiPattern = 0;
for(uint32_t iPhi = 0; iPhi < 3; iPhi++){
bool activeStrip = false;

for(uint32_t iEta = 0; iEta < 3; iEta++){
  if (object->boostedJetRegionET()[3*iEta+iPhi] > 30 && object->boostedJetRegionET()[3*iEta+iPhi] > object->et()*0.0625  && object->boostedJetRegionTauVeto()[3*iEta+iPhi]==0 ) activeStrip = true;
            }
  if(activeStrip) activeRegionPhiPattern |= (0x1 << iPhi);
        }

string regionEta1 = activeRegionEtaPattern.to_string<char,std::string::traits_type,std::string::allocator_type>();
string regionPhi1 = activeRegionPhiPattern.to_string<char,std::string::traits_type,std::string::allocator_type>();

regionEta.push_back(activeRegionEtaPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());
regionPhi.push_back(activeRegionPhiPattern.to_string<char,std::string::traits_type,std::string::allocator_type>());
//cout<<"regionEta1="<<regionEta1<<endl;

        if ( (regionEta1=="100" && (regionPhi1 =="100" ||  regionPhi1== "001" || regionPhi1== "010" || regionPhi1== "110" )) ||
            (regionEta1=="001" && (regionPhi1 =="100" ||  regionPhi1== "001" || regionPhi1== "010" || regionPhi1== "110" ) ) ||
            (regionEta1=="010" && (regionPhi1 =="100" ||  regionPhi1== "001" || regionPhi1== "010" || regionPhi1== "110" ) )
        ){
            TauCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));
        }
    }
}

/*cout<<"nevent000= "<<nevent000<<endl;
cout<<"nevent0= "<<nevent0<<endl;
cout<<"nevent1= "<<nevent1<<endl;
cout<<"nevent2= "<<nevent2<<endl;
cout<<"nevent3= "<<nevent3<<endl;
cout<<"nevent4= "<<nevent4<<endl;
cout<<"nevent5= "<<nevent5<<endl;
cout<<"nevent6= "<<nevent6<<endl;

cout<<"nevent10= "<<nevent10<<endl;
cout<<"nevent11= "<<nevent11<<endl;
cout<<"nevent12= "<<nevent12<<endl;
cout<<"nevent13= "<<nevent13<<endl;
cout<<"nevent14= "<<nevent14<<endl;
cout<<"nevent15= "<<nevent15<<endl;
cout<<"nevent16= "<<nevent16<<endl;

cout<<"nevent20= "<<nevent20<<endl;
cout<<"nevent21= "<<nevent21<<endl;
cout<<"nevent22= "<<nevent22<<endl;
cout<<"nevent23= "<<nevent23<<endl;
cout<<"nevent24= "<<nevent24<<endl;
cout<<"nevent25= "<<nevent25<<endl;
cout<<"nevent26= "<<nevent26<<endl;

cout<<"nevent30= "<<nevent30<<endl;
cout<<"nevent31= "<<nevent31<<endl;
cout<<"nevent32= "<<nevent32<<endl;
cout<<"nevent33= "<<nevent33<<endl;
cout<<"nevent34= "<<nevent34<<endl;
cout<<"nevent35= "<<nevent35<<endl;
cout<<"nevent36= "<<nevent36<<endl;

cout<<"nevent40= "<<nevent40<<endl;
cout<<"nevent41= "<<nevent41<<endl;
cout<<"nevent42= "<<nevent42<<endl;
cout<<"nevent43= "<<nevent43<<endl;
cout<<"nevent44= "<<nevent44<<endl;
cout<<"nevent45= "<<nevent45<<endl;
cout<<"nevent46= "<<nevent46<<endl;

cout<<"nevent50= "<<nevent50<<endl;
cout<<"nevent51= "<<nevent51<<endl;
cout<<"nevent52= "<<nevent52<<endl;
cout<<"nevent53= "<<nevent53<<endl;
cout<<"nevent54= "<<nevent54<<endl;
cout<<"nevent55= "<<nevent55<<endl;
cout<<"nevent56= "<<nevent56<<endl;

cout<<"nevent60= "<<nevent60<<endl;
cout<<"nevent61= "<<nevent61<<endl;
cout<<"nevent62= "<<nevent62<<endl;
cout<<"nevent63= "<<nevent63<<endl;
cout<<"nevent64= "<<nevent64<<endl;
cout<<"nevent65= "<<nevent65<<endl;
cout<<"nevent66= "<<nevent66<<endl;
cout<<"All= "<<nevent0+nevent1+nevent2+nevent3+nevent4+nevent5+nevent6+
nevent10+nevent11+nevent12+nevent13+nevent14+nevent15+nevent16+
nevent20+nevent21+nevent22+nevent23+nevent24+nevent25+nevent26+
nevent30+nevent31+nevent32+nevent33+nevent34+nevent35+nevent36+
nevent40+nevent41+nevent42+nevent43+nevent44+nevent45+nevent46+
nevent50+nevent51+nevent52+nevent53+nevent54+nevent55+nevent56+
nevent60+nevent61+nevent62+nevent63+nevent64+nevent65+nevent66<<endl;
*/

    for (size_t j = 0; j < TauCands->size(); j++){
        float pt = TauCands->at(j).pt();
        float eta = TauCands->at(j).eta();
        float phi= TauCands->at(j).phi();
        float et = TauCands->at(j).et();
        TLorentzVector temp ;
        temp.SetPtEtaPhiE(pt,eta,phi,et);
        L1boostedTaus->push_back(temp);

    }




    edm::Handle<vector<pat::Tau> > taus;
    if(evt.getByToken(tauSrc_, taus)){
      for (const pat::Tau &itau : *taus) {
        if(itau.pt() > recoPt_ ) {
          goodTaus.push_back(itau);
          TLorentzVector temp ;
          temp.SetPtEtaPhiE(itau.pt(),itau.eta(),itau.phi(),itau.et());
          boostedTaus->push_back(temp);
        }
      }
    }


  edm::Handle<EcalTrigPrimDigiCollection> ecalTPs;
  evt.getByToken(ecalTPToken_, ecalTPs);
  for ( const auto& ecalTp : *ecalTPs ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int et = ecalTp.compressedEt();
    if(et != 0) {
      float eta = getRecoEtaNew(caloEta);
      float phi = getRecoPhiNew(caloPhi);
      TLorentzVector temp;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      ecalTPGs->push_back(temp);
    }
  }

  edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
  evt.getByToken(hcalTPToken_, hcalTPs);
  for ( const auto& hcalTp : *hcalTPs ) {
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    if(absCaloEta == 29) {
      continue;
    }
    else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
      continue;
    }
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      if(caloPhi <= 72) {
	int et = hcalTp.SOI_compressedEt();
	if(et != 0) {
          float eta = getRecoEtaNew(caloEta);
          float phi = getRecoPhiNew(caloPhi);
          TLorentzVector temp;
          temp.SetPtEtaPhiE(et,eta,phi,et);
          hcalTPGs->push_back(temp);
	}
      }
      else {
	std::cerr << "Illegal Tower: caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;
      }
    }
    else {
      std::cerr << "Illegal Tower: caloEta = " << caloEta << std::endl;
    }
  }

  passL1SingleJet180 = false; passL1HTT360 = false;

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  if(!evt.getByToken(trigobjectsMINIAODToken_, triggerObjects)) cout<<"ERROR GETTING TRIGGER OBJECTS"<<endl;
  evt.getByToken(trigobjectsMINIAODToken_, triggerObjects);

  edm::Handle<edm::TriggerResults> trigResults;
  if(!evt.getByToken(trgresultsToken_, trigResults)) cout<<"ERROR GETTING TRIGGER RESULTS"<<endl;
  evt.getByToken(trgresultsToken_, trigResults);

  const edm::TriggerNames &names = evt.triggerNames(*trigResults);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackFilterLabels(evt,*trigResults);
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      string myfillabl=obj.filterLabels()[h];
      if(myfillabl=="hltL1sSingleJet180") { passL1SingleJet180 = true; }
      if(myfillabl=="hltL1sHTT380erIorHTT320er") { passL1HTT360 = true; }
    }
  }

  // Accessing existing L1 seed stored in MINIAOD
  edm::Handle<l1t::JetBxCollection> stage2Jets;
  if(!evt.getByToken(stage2JetToken_, stage2Jets)) cout<<"ERROR GETTING THE STAGE 2 JETS"<<endl;
  evt.getByToken(stage2JetToken_, stage2Jets);
  for (l1t::JetBxCollection::const_iterator obj = stage2Jets->begin(0); obj != stage2Jets->end(0); obj++) {
        seeds.push_back(*obj);
        TLorentzVector temp;
        temp.SetPtEtaPhiE(obj->pt(), obj->eta(), obj->phi(), obj->et());
        seed180->push_back(temp);
  }

  edm::Handle<l1t::EGammaBxCollection> stage2EGs;
  if(!evt.getByToken(stage2EGToken_, stage2EGs)) cout<<"ERROR GETTING THE STAGE 2 EGAMMAS"<<endl;
  evt.getByToken(stage2EGToken_, stage2EGs);
  for (l1t::EGammaBxCollection::const_iterator obj = stage2EGs->begin(0); obj != stage2EGs->end(0); obj++) {
        TLorentzVector temp;
        temp.SetPtEtaPhiE(obj->pt(), obj->eta(), obj->phi(), obj->et());
        egseed->push_back(temp);
  }

  edm::Handle<l1t::TauBxCollection> stage2Taus;
  if(!evt.getByToken(stage2TauToken_, stage2Taus)) cout<<"ERROR GETTING THE STAGE 2 TAUS"<<endl;
  evt.getByToken(stage2TauToken_, stage2Taus);
  for (l1t::TauBxCollection::const_iterator obj = stage2Taus->begin(0); obj != stage2Taus->end(0); obj++) {
        TLorentzVector temp;
        temp.SetPtEtaPhiE(obj->pt(), obj->eta(), obj->phi(), obj->et());
        tauseed->push_back(temp);
  }

  edm::Handle<l1t::EtSumBxCollection> stage2EtSum;
  if(!evt.getByToken(stage2EtSumToken_, stage2EtSum)) cout<<"ERROR GETTING THE STAGE 2 ETSUM"<<endl;
  evt.getByToken(stage2EtSumToken_, stage2EtSum);
  for (l1t::EtSumBxCollection::const_iterator obj = stage2EtSum->begin(0); obj != stage2EtSum->end(0); obj++) {
    if(obj->getType() == l1t::EtSum::kTotalHt){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(obj->pt(), obj->eta(), obj->phi(), obj->et());
      htseed->push_back(temp);
    }
  }

  // Accessing L1boosted collection
  edm::Handle<vector<l1extra::L1JetParticle>> l1Boosted;
  if(!evt.getByToken(l1BoostedToken_, l1Boosted)) cout<<"ERROR GETTING THE L1BOOSTED JETS"<<endl;
  evt.getByToken(l1BoostedToken_, l1Boosted);
  const vector<l1extra::L1JetParticle> &l1B = *l1Boosted;
  for(auto obj : l1B) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.et());
    l1Jets->push_back(temp);
  }

  // Start Running Analysis
  Handle<vector<reco::CaloJet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
    for (const reco::CaloJet &jet : *jets) {
      if(jet.pt() > recoPt_ ) {
	goodJets.push_back(jet);
      }
    }
  }
  else
    int uu=0;
//    cout<<"Error getting calo jets"<<endl;

  Handle<vector<pat::Jet> > jetsAK8;

  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting AK8 Jets
    for (const pat::Jet &jetAK8 : *jetsAK8) {
      if(jetAK8.pt() > recoPt_ ) {
        nSubJets.push_back(jetAK8.subjets("SoftDropPuppi").size());
        nBHadrons.push_back(jetAK8.jetFlavourInfo().getbHadrons().size());
        TLorentzVector temp ;
        temp.SetPtEtaPhiE(jetAK8.pt(),jetAK8.eta(),jetAK8.phi(),jetAK8.et());
        ak8Jets->push_back(temp);
        if(jetAK8.subjets("SoftDropPuppi").size() ==  2 /*&& jetAK8.jetFlavourInfo().getbHadrons().size() > 1*/){
          goodJetsAK8.push_back(jetAK8);
        }
      }
    }
  }
  else
    cout<<"Error getting AK8 jets"<<endl;






  zeroOutAllVariables();
//  if(boostedTaus->size()>0){
  if(goodJetsAK8.size()>0){

    for(auto jet:goodJetsAK8){
      tau1.push_back(jet.userFloat("NjettinessAK8Puppi:tau1"));
      tau2.push_back(jet.userFloat("NjettinessAK8Puppi:tau2"));
      tau3.push_back(jet.userFloat("NjettinessAK8Puppi:tau3"));
      HFlav.clear();
      for(unsigned int isub=0; isub<((jet.subjets("SoftDropPuppi")).size()); isub++){
        HFlav.push_back(jet.subjets("SoftDropPuppi")[isub]->hadronFlavour());
        TLorentzVector temp;
        temp.SetPtEtaPhiE(jet.subjets("SoftDropPuppi")[isub]->pt(),jet.subjets("SoftDropPuppi")[isub]->eta(),jet.subjets("SoftDropPuppi")[isub]->phi(),jet.subjets("SoftDropPuppi")[isub]->et());
        subJets->push_back(temp);
      }
      subJetHFlav.push_back(HFlav);
      //take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    }
//}
    //Match to boosted jets and see if we can match subjettiness functions...
    vector<l1extra::L1JetParticle> l1JetsSorted;
    for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1Boosted->begin(); l1Jet != l1Boosted->end(); l1Jet++ ){
      l1JetsSorted.push_back(*l1Jet);
    }
    if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}

    pat::Jet recoJet_1;

    recoPt_1  = goodJetsAK8.at(0).pt();
    recoEta_1 = goodJetsAK8.at(0).eta();
    recoPhi_1 = goodJetsAK8.at(0).phi();
    recoJet_1 = goodJetsAK8.at(0);

  /*  pat::Tau recoJet_1;

    recoPt_1  = boostedTaus->at(0).Pt();
    recoEta_1 = boostedTaus->at(0).Eta();
    recoPhi_1 = boostedTaus->at(0).Phi();
    recoJet_1 = goodTaus.at(0);
*/

    int i = 0;
    int foundL1Jet_1 = 0;
    l1extra::L1JetParticle l1Jet_1;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        if(reco::deltaR(jet, recoJet_1)<0.4 && foundL1Jet_1 == 0 ){
          l1Jet_1 = jet;
          l1Pt_1  = jet.pt();
          l1Eta_1 = jet.eta();
          l1Phi_1 = jet.phi();
          l1NthJet_1 = i;
          foundL1Jet_1 = 1;
        }
        i++;
      }
    }

    int j = 0;
    int foundSeed_1 = 0;
    if(seeds.size() > 0){
      for(auto seed : seeds){
        if(reco::deltaR(seed, recoJet_1)<0.4 && foundSeed_1 == 0 ){
          seedPt_1  = seed.pt();
          seedEta_1 = seed.eta();
          seedPhi_1 = seed.phi();
          seedNthJet_1 = j;
          foundSeed_1 = 1;
        }
        j++;
      }
    }

  }

//  edm::Handle<reco::GenParticleCollection> genParticles;
  if(evt.getByToken(genSrc_, genParticles)){//Begin Getting Gen Particles
    for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      double DR = reco::deltaR(recoEta_1, recoPhi_1, genparticle->eta(), genparticle->phi());
      if (DR < genDR && genparticle->status() > 21 && genparticle->status() < 41){
        genDR = DR;
        genId = genparticle->pdgId();
        genMother = genparticle->motherRef(0)->pdgId();
        genPt_1 = genparticle->pt();
        genEta_1 = genparticle->eta();
        genPhi_1 = genparticle->phi();
        genM_1 = genparticle->mass();
      }
    }
  }



  efficiencyTree->Fill();
}

void BoostedJetStudies::zeroOutAllVariables(){
  genPt_1=-99; genEta_1=-99; genPhi_1=-99; genM_1=-99; genDR=99; genId=-99; genMother=-99;
  seedPt_1=-99; seedEta_1=-99; seedPhi_1=-99; seedNthJet_1=-99;
  recoPt_1=-99; recoEta_1=-99; recoPhi_1=-99; recoNthJet_1=-99;
  l1Pt_1=-99; l1Eta_1=-99; l1Phi_1=-99; l1NthJet_1=-99;
  recoPt_=-99;
}

void BoostedJetStudies::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,    "lumi/I");
    tree->Branch("event",   &event,   "event/I");
    tree->Branch("nPV",     &nPV,     "nPV/I");
    tree->Branch("genPt_1",       &genPt_1,     "genPt_1/D");
    tree->Branch("genEta_1",      &genEta_1,    "genEta_1/D");
    tree->Branch("genPhi_1",      &genPhi_1,    "genPhi_1/D");
    tree->Branch("genM_1",        &genM_1,      "genM_1/D");
    tree->Branch("genDR",         &genDR,       "genDR/D");
    tree->Branch("genId",         &genId,       "genId/I");
    tree->Branch("genMother",     &genMother,   "genMother/I");
    tree->Branch("seedPt_1",      &seedPt_1,     "seedPt_1/D");
    tree->Branch("seedEta_1",     &seedEta_1,    "seedEta_1/D");
    tree->Branch("seedPhi_1",     &seedPhi_1,    "seedPhi_1/D");
    tree->Branch("seedNthJet_1",  &seedNthJet_1, "seedNthJet_1/I");
    tree->Branch("recoPt_1",      &recoPt_1,     "recoPt_1/D");
    tree->Branch("recoEta_1",     &recoEta_1,    "recoEta_1/D");
    tree->Branch("recoPhi_1",     &recoPhi_1,    "recoPhi_1/D");
    tree->Branch("recoNthJet_1",  &recoNthJet_1, "recoNthJet_1/I");
    tree->Branch("l1Pt_1",        &l1Pt_1,       "l1Pt_1/D");
    tree->Branch("l1Eta_1",       &l1Eta_1,      "l1Eta_1/D");
    tree->Branch("l1Phi_1",       &l1Phi_1,      "l1Phi_1/D");
    tree->Branch("l1NthJet_1",    &l1NthJet_1,   "l1NthJet_1/I");
    tree->Branch("tau1",          &tau1);
    tree->Branch("tau2",          &tau2);
    tree->Branch("tau3",          &tau3);
    tree->Branch("nSubJets",      &nSubJets);
    tree->Branch("subJetHFlav",   &subJetHFlav);
    tree->Branch("nBHadrons",     &nBHadrons);
    tree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0);
    tree->Branch("seed180", "vector<TLorentzVector>", &seed180, 32000, 0);
    tree->Branch("egseed", "vector<TLorentzVector>", &egseed, 32000, 0);
    tree->Branch("tauseed", "vector<TLorentzVector>", &tauseed, 32000, 0);
    tree->Branch("htseed", "vector<TLorentzVector>", &htseed, 32000, 0);
    tree->Branch("ak8Jets", "vector<TLorentzVector>", &ak8Jets, 32000, 0);
    tree->Branch("subJets", "vector<TLorentzVector>", &subJets, 32000, 0);
    tree->Branch("hcalTPGs", "vector<TLorentzVector>", &hcalTPGs, 32000, 0);
    tree->Branch("ecalTPGs", "vector<TLorentzVector>", &ecalTPGs, 32000, 0);
    tree->Branch("passL1HTT360",  &passL1HTT360,  "passL1HTT360/B");
    tree->Branch("passL1SingleJet180",  &passL1SingleJet180,  "passL1SingleJet180/B");
    tree->Branch("boostedTaus", "vector<TLorentzVector>", &boostedTaus, 32000, 0);
    tree->Branch("L1boostedTaus", "vector<TLorentzVector>", &L1boostedTaus, 32000, 0);
    tree->Branch("regionEta",       &regionEta);
    tree->Branch("regionPhi",       &regionPhi);



  }


// ------------ method called once each job just before starting event loop  ------------
void
BoostedJetStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
BoostedJetStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

// ------------ method called when ending the processing of a run  ------------
/*
  void
  BoostedJetStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BoostedJetStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BoostedJetStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedJetStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedJetStudies);
