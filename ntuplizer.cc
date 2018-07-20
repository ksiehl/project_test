// -*- C++ -*-
//
// Package:    Analysis/ntuplizer
// Class:      ntuplizer
// 
/**\class ntuplizer ntuplizer.cc Analysis/ntuplizer/plugins/ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  kevin siehl
//         Created:  Fri, 23 Jun 2017 19:05:21 GMT
//
//

//#define MULTI_WEIGHTS 
// system include files
#include <memory>
#include <iostream>
#include <fstream> //not auto-gen
#include <algorithm>
#include <bitset>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

////////////////////////

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen particle
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Electron
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"

// Pileup
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// LHE weights
//#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Utilities
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

////////////////////////////////////////////////////////

//DataFormats
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//basic root stuff
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <TRandom3.h>

//and something else
#include "Analysis/ntuplizer/metZcalc/METzCalculator.cc"
//and for n-subjettiness
//#include "/uscms_data/d3/ksiehl/CMSSW_8_0_26/src/Analysis/ntuplizer/jet_tools/JetTools.cc"
//#include "RecoJets/JetProducers/interface/NjettinessAdder.h"

struct particle
{
  Double_t Px;
  Double_t Py;
  Double_t Pz;
  Double_t E;
  Double_t M;
  Double_t Pt;
  Double_t Eta;
  Double_t Theta;
  Double_t Phi;
  Double_t Y;
};




Double_t d_costheta1 = 0.0;
Double_t d_costheta2 = 0.0;
Double_t d_phi = 0.0;
Double_t d_costhetastar = 0.0;
Double_t d_phi1 = 0.0;
Double_t d_phi2 = 0.0;


Double_t d_leptons_in_lep_px = 0.0;
Double_t d_leptons_in_lep_py = 0.0;
Double_t d_leptons_in_lep_pz = 0.0;
Double_t d_partons_in_lep_px = 0.0;
Double_t d_partons_in_lep_py = 0.0;
Double_t d_partons_in_lep_pz = 0.0;
Double_t d_parton1_in_lep_px = 0.0;
Double_t d_parton2_in_lep_px = 0.0;
Double_t d_parton1_in_lep_py = 0.0;
Double_t d_parton2_in_lep_py = 0.0;
Double_t d_parton1_in_lep_pz = 0.0;
Double_t d_parton2_in_lep_pz = 0.0;
Double_t d_lepton1_in_lep_px = 0.0;
Double_t d_lepton1_in_lep_py = 0.0;
Double_t d_lepton1_in_lep_pz = 0.0;
Double_t d_lepton1_dotted_x = 0.0;
Double_t d_lepton1_dotted_y = 0.0;
Double_t d_lepton1_dotted_z = 0.0;
Double_t d_leptons_in_had_px = 0.0;
Double_t d_leptons_in_had_py = 0.0;
Double_t d_leptons_in_had_pz = 0.0;
Double_t d_lepton1_in_had_px = 0.0;
Double_t d_lepton1_in_had_py = 0.0;
Double_t d_lepton1_in_had_pz = 0.0;
Double_t d_lepton2_in_had_px = 0.0;
Double_t d_lepton2_in_had_py = 0.0;
Double_t d_lepton2_in_had_pz = 0.0;
Double_t d_parton1_in_had_px = 0.0;
Double_t d_parton1_in_had_py = 0.0;
Double_t d_parton1_in_had_pz = 0.0;
Double_t d_parton1_dotted_x = 0.0;
Double_t d_parton1_dotted_y = 0.0;
Double_t d_parton1_dotted_z = 0.0;
Double_t d_complicated1_px = 0.0;
Double_t d_complicated1_py = 0.0;
Double_t d_complicated1_pz = 0.0;
Double_t d_complicated2_px = 0.0;
Double_t d_complicated2_py = 0.0;
Double_t d_complicated2_pz = 0.0;
Double_t d_lepton1WWframe_X = 0.0;
Double_t d_lepton1WWframe_Y = 0.0;
Double_t d_lepton1WWframe_Z = 0.0;
Double_t d_lepton_sumWWframe_X = 0.0;
Double_t d_lepton_sumWWframe_Y = 0.0;
Double_t d_lepton_sumWWframe_Z = 0.0;
Double_t d_parton1WWframe_X = 0.0;
Double_t d_parton1WWframe_Y = 0.0;
Double_t d_parton1WWframe_Z = 0.0;
Double_t d_parton_sumWWframe_X = 0.0;
Double_t d_parton_sumWWframe_Y = 0.0;
Double_t d_parton_sumWWframe_Z = 0.0;
Double_t d_boostWWframe_X = 0.0;
Double_t d_boostWWframe_Y = 0.0;
Double_t d_boostWWframe_Z = 0.0;
Double_t d_boostWlep_X = 0.0;
Double_t d_boostWlep_Y = 0.0;
Double_t d_boostWlep_Z = 0.0;
Double_t d_boostWhad_X = 0.0;
Double_t d_boostWhad_Y = 0.0;
Double_t d_boostWhad_Z = 0.0;


int extra_muon_count = 0;
int extra_electron_count = 0;
int no_lepton_found_count = 0;
int loose_muon_count = 0;
int loose_electron_count = 0;
int no_jet_count = 0;
int merged_combo_count = 0;
int extra_merged_combo_count = 0;
int missing_et_cutoff_count = 0;

int nonmerged_count = 0;
int inverted_count = 0;
int neither_count = 0;

int skip_counter = 0;
int backup_skip_counter = 0;
int merged_counter = 0;

int merged_event_count = 0;

bool talking = false;

const Double_t Wmass = 80.4; //other file has a 4bodymass function here, as well as namespace lhapdf functions
const Double_t Zmass = 90.2;

const Double_t el_mass = 0.00051099891;
const Double_t mu_mass = 0.105658367;

void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Double_t& costheta1, Double_t& costheta2, Double_t& phi, Double_t& costhetastar, Double_t& phi1, Double_t& phi2);

void intermediate_steps(TLorentzVector lepton1, TLorentzVector lepton2, TLorentzVector parton1, TLorentzVector parton2, Double_t& leptons_in_lep_px, Double_t& leptons_in_lep_py, Double_t& leptons_in_lep_pz, Double_t& partons_in_lep_px, Double_t& partons_in_lep_py, Double_t& partons_in_lep_pz, Double_t& parton1_in_lep_px, Double_t& parton2_in_lep_px, Double_t& parton1_in_lep_py, Double_t& parton2_in_lep_py, Double_t& parton1_in_lep_pz, Double_t& parton2_in_lep_pz, Double_t& lepton1_in_lep_px, Double_t& lepton1_in_lep_py, Double_t& lepton1_in_lep_pz, Double_t& lepton1_dotted_x, Double_t& lepton1_dotted_y, Double_t& lepton1_dotted_z, Double_t& leptons_in_had_px, Double_t& leptons_in_had_py, Double_t& leptons_in_had_pz, Double_t& lepton1_in_had_px, Double_t& lepton1_in_had_py, Double_t& lepton1_in_had_pz, Double_t& lepton2_in_had_px, Double_t& lepton2_in_had_py, Double_t& lepton2_in_had_pz, Double_t& parton1_in_had_px, Double_t& parton1_in_had_py, Double_t& parton1_in_had_pz, Double_t& parton1_dotted_x, Double_t& parton1_dotted_y, Double_t& parton1_dotted_z, Double_t& complicated1_px, Double_t& complicated1_py, Double_t& complicated1_pz, Double_t& complicated2_px, Double_t& complicated2_py, Double_t& complicated2_pz, Double_t& lepton_sumWWframe_X, Double_t& lepton_sumWWframe_Y, Double_t& lepton_sumWWframe_Z, Double_t& lepton1WWframe_X, Double_t& lepton1WWframe_Y, Double_t& lepton1WWframe_Z, Double_t& parton_sumWWframe_X, Double_t& parton_sumWWframe_Y, Double_t& parton_sumWWframe_Z, Double_t& parton1WWframe_X, Double_t& parton1WWframe_Y, Double_t& parton1WWframe_Z, Double_t& costhetastar, Double_t& costheta1, Double_t& phi, Double_t& costheta2, Double_t& phi1, Double_t& phi2, Double_t& boostWWframe_X, Double_t& boostWWframe_Y, Double_t& boostWWframe_Z, Double_t& boostWlep_X, Double_t& boostWlep_Y, Double_t& boostWlep_Z, Double_t& boostWhad_X, Double_t& boostWhad_Y, Double_t& boostWhad_Z);

double jetCharge(const reco::PFJet &jet, const double kappa)
{
  double charge = 0;
  const unsigned int nPFCands = jet.nConstituents ();
  for(unsigned int ipf = 0; ipf < nPFCands; ipf++)
    {
      const reco::Candidate* cand = jet.getJetConstituentsQuick ()[ipf];
      const pat::PackedCandidate *packCand = dynamic_cast<const pat::PackedCandidate *>(cand);
      if(packCand != 0)
	{ 
	  if(packCand->charge() == 0) continue;
	  charge += pow(packCand->pt(),kappa)*packCand->charge();
	}
      else
	{ 
	  const reco::PFCandidatePtr pfcand    = jet.getPFConstituents().at(ipf);
	  const reco::TrackRef       track  = pfcand->trackRef();
	  if(track.isNull()) continue;
	  charge += pow(pfcand->pt(),kappa)*track->charge();
	}
    }
  double denom = (kappa==1) ? jet.pt() : pow(jet.pt(),kappa);
  return charge/denom;
}
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit ntuplizer(const edm::ParameterSet&);
  ~ntuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::ofstream outfile;


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  
  edm::EDGetTokenT<pat::JetCollection> ak4jetToken_;
  edm::EDGetTokenT<pat::JetCollection> ak8jetToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetToken_;
  edm::EDGetTokenT<pat::JetCollection> ak8CHSSoftDropSubjetsToken_;
  //edm::EDGetTokenT<pat::JetCollection> ak8PuppiSoftDropSubjetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak4genjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak8genjetToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<double> rhoToken_;  
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  //edm::EDGetTokenT<edm::TriggerResults> triggerResultsMETFilterToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  //edm::EDGetTokenT<bool> badMuonFilterToken_;
  //edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  //edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
  //edm::EDGetTokenT<LHEEventProduct> theSrc_;
  //edm::EDGetTokenT<GenEventInfoProduct> pdfToken_;
  //edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;                 
  //edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

  bool isRun2016F_;
  
  //std::vector<std::string>  jecPayloadsAK4chs_;
  //std::vector<std::string>  jecPayloadsAK8chs_;
  //std::vector<std::string>  jecPayloadsAK4pup_;
  //std::vector<std::string>  jecPayloadsAK8pup_;

  std::vector<std::string>  jecPayloadsAK4chsSecondary_;
  std::vector<std::string>  jecPayloadsAK8chsSecondary_;
  std::vector<std::string>  jecPayloadsAK4pupSecondary_;
  std::vector<std::string>  jecPayloadsAK8pupSecondary_;

  //std::string jertextAK4_;   // jet resolution AK4 jets
  //std::string jertextAK8_;   // jet resolution AK8 jets
  //std::string jerSFtext_ ;   // jer SF
  
  boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4chs;
  //boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8chs;
  //boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4pup;
  //boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8pup;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4chs;
  //boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8chs;
  //boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4pup;
  //boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8pup;

  std::string jet_scale_fac_textfile_;
  std::string ak4_textfile_;
  std::string ak8_textfile_;
  
  TH1D* wp_mass; //a histogram is necesarry to cause the output to be in the right format;
  TH1D* wm_mass; //a histogram is necesarry to cause the output to be in the right format;
  //TH1D* pt_smear; //difference between MET pt and neutrino_gen pt
  TH1D *comp_hist0, *comp_hist1, *comp_hist2, *comp_hist3, *comp_hist4, *comp_hist5, *comp_hist6, *comp_hist7;
  TH1D *ww_mass, *ww_pt, *ww_rap;
  

  TH1D *wp_pt, *wm_pt, *whad_pt, *wlep_pt, *neu_pt, *elr_pt, *mur_pt, *elg_pt, *mug_pt, *met_pt;
  TH1D *wp_rap, *wm_rap, *whad_rap, *wlep_rap, *neu_rap, *elr_rap, *mur_rap, *elg_rap, *mug_rap, *met_rap;
  
  TH1D *whad_pt_smear, *wlep_pt_smear, *el_pt_smear, *mu_pt_smear, *met_pt_smear, *parton0_pt_smear, *parton1_pt_smear;
  TH1D *whad_phi_smear, *wlep_phi_smear, *el_phi_smear, *mu_phi_smear, *met_phi_smear, *parton0_phi_smear, *parton1_phi_smear;
  TH1D *whad_eta_smear, *wlep_eta_smear, *el_eta_smear, *mu_eta_smear, *met_eta_smear, *parton0_eta_smear, *parton1_eta_smear;
  
  TH1D *wminus_theta, *wplus_theta, *lepton_theta, *neutrino_theta, *quark_theta, *antiquark_theta;
  TH1D *wminus_phi, *wplus_phi, *lepton_phi, *neutrino_phi, *quark_phi, *antiquark_phi;
  
  //other file has many more bools, string vectors, shared pointers, string, tfile, th1ds, ints


  //MULTIPLE WEIGHTS
  edm::EDGetTokenT<Double_t> mgreweight1_;
  edm::EDGetTokenT<Double_t> mgreweight2_;
  edm::EDGetTokenT<Double_t> mgreweight3_;
  edm::EDGetTokenT<Double_t> mgreweight4_;
  edm::EDGetTokenT<Double_t> mgreweight5_;
  edm::EDGetTokenT<Double_t> mgreweight6_;
  edm::EDGetTokenT<Double_t> mgreweight7_;
  /*edm::EDGetTokenT<Double_t> mgreweight8_;
  edm::EDGetTokenT<Double_t> mgreweight9_;
  edm::EDGetTokenT<Double_t> mgreweight10_;
  edm::EDGetTokenT<Double_t> mgreweight11_;
  edm::EDGetTokenT<Double_t> mgreweight12_;
  edm::EDGetTokenT<Double_t> mgreweight13_;
  edm::EDGetTokenT<Double_t> mgreweight14_;
  edm::EDGetTokenT<Double_t> mgreweight15_;
  edm::EDGetTokenT<Double_t> mgreweight16_;
  edm::EDGetTokenT<Double_t> mgreweight17_;
  edm::EDGetTokenT<Double_t> mgreweight18_;
  edm::EDGetTokenT<Double_t> mgreweight19_;
  edm::EDGetTokenT<Double_t> mgreweight20_;
  / *edm::EDGetTokenT<Double_t> mgreweight21_;
  edm::EDGetTokenT<Double_t> mgreweight22_;
  edm::EDGetTokenT<Double_t> mgreweight23_;
  edm::EDGetTokenT<Double_t> mgreweight24_;
  edm::EDGetTokenT<Double_t> mgreweight25_;
  edm::EDGetTokenT<Double_t> mgreweight26_;
  edm::EDGetTokenT<Double_t> mgreweight27_;
  edm::EDGetTokenT<Double_t> mgreweight28_;
  edm::EDGetTokenT<Double_t> mgreweight29_;
  edm::EDGetTokenT<Double_t> mgreweight30_;
  edm::EDGetTokenT<Double_t> mgreweight31_;
  edm::EDGetTokenT<Double_t> mgreweight32_;
  edm::EDGetTokenT<Double_t> mgreweight33_;
  edm::EDGetTokenT<Double_t> mgreweight34_;
  edm::EDGetTokenT<Double_t> mgreweight35_;
  edm::EDGetTokenT<Double_t> mgreweight36_;
  edm::EDGetTokenT<Double_t> mgreweight37_;
  edm::EDGetTokenT<Double_t> mgreweight38_;
  edm::EDGetTokenT<Double_t> mgreweight39_;
  edm::EDGetTokenT<Double_t> mgreweight40_;
  edm::EDGetTokenT<Double_t> mgreweight41_;
  edm::EDGetTokenT<Double_t> mgreweight42_;
  edm::EDGetTokenT<Double_t> mgreweight43_;
  edm::EDGetTokenT<Double_t> mgreweight44_;
  edm::EDGetTokenT<Double_t> mgreweight45_;
  edm::EDGetTokenT<Double_t> mgreweight46_;
  edm::EDGetTokenT<Double_t> mgreweight47_;
  edm::EDGetTokenT<Double_t> mgreweight48_;
  edm::EDGetTokenT<Double_t> mgreweight49_;
  edm::EDGetTokenT<Double_t> mgreweight50_;*/

  TTree *MyTree; //other file had some vectors and a string before listings

  std::vector<std::string> *Trig_Names = new std::vector<std::string>;
  std::vector<int> *Trig_Prescales = new std::vector<int>;
  std::vector<bool> *Trig_Pass = new std::vector<bool>;

  std::string Trig_Accpt_Bits;

  Int_t good_events = 0;

  Double_t min_mass_diff = -99.9;
  Double_t min_submass_diff = -99.9;
  Double_t hadronic_gen_2body_mass = -99.9;
  Double_t leptonic_2body_mass = -99.9;
  Double_t leptonic_gen_2body_mass = -99.9;

  Double_t min_quark_delta = -99.9;
  Double_t min_antiquark_delta = -99.9;

  Double_t gen_konstant = -99.9;
  Double_t reco_konstant = -99.9;

  particle Wpls_gen;

  particle Wmns_gen;

  particle Zntrl_gen;

  particle qrk_gen;
 
  particle antiqrk_gen;

  particle lptn_gen;

  particle ntrno_gen;

  Double_t lepton_deltaR = -99.9;
  Double_t lepton_deltaPhi = -99.9;

  Double_t neutrino_deltaR = -99.9;
  Double_t other_neutrino_deltaR = -99.9;
  Double_t truth_neutrino_deltaR = -99.9;
  Double_t false_neutrino_deltaR = -99.9;

  Double_t neutrino_deltaPhi = -99.9;
  Double_t other_neutrino_deltaPhi = -99.9;
  Double_t neutrino_deltaPt = -99.9;

  Double_t jet0_deltaR = -99.9;
  Double_t jet1_deltaR = -99.9;

  Double_t subjet0_deltaR = -99.9;
  Double_t subjet1_deltaR = -99.9;
  
  Double_t costheta1 = -99.9;
  Double_t costheta2 = -99.9;
  Double_t costhetastar = -99.9;
  Double_t phi1 = -99.9;
  Double_t phi2 = -99.9;
  Double_t phi = -99.9;

  Double_t costheta1_gen = -99.9;
  Double_t costheta2_gen = -99.9;
  Double_t costhetastar_gen = -99.9;
  Double_t phi1_gen = -99.9;
  Double_t phi2_gen = -99.9;
  Double_t phi_gen = -99.9;

  Double_t other_costheta1 = -99.9;
  Double_t other_costheta2 = -99.9;
  Double_t other_costhetastar = -99.9;
  Double_t other_phi1 = -99.9;
  Double_t other_phi2 = -99.9;
  Double_t other_phi = -99.9;

  Double_t truth_costheta1 = -99.9;
  Double_t truth_costheta2 = -99.9;
  Double_t truth_costhetastar = -99.9;
  Double_t truth_phi1 = -99.9;
  Double_t truth_phi2 = -99.9;
  Double_t truth_phi = -99.9;

  Double_t false_costheta1 = -99.9;
  Double_t false_costheta2 = -99.9;
  Double_t false_costhetastar = -99.9;
  Double_t false_phi1 = -99.9;
  Double_t false_phi2 = -99.9;
  Double_t false_phi = -99.9;

  Double_t Tau21 = -99.9; Double_t Tau32 = -99.9; Double_t PTau21 = -99.9; Double_t PTau32 = -99.9;

  //LONG CHAIN
  particle Whad_reco, Wlep_reco, ntrno_reco, other_ntrno_reco, truth_ntrno_reco, false_ntrno_reco, muon_reco, electrn_reco, lepton_reco, jt0, jt1, subjt0, subjt1;

  Int_t extra_muon = 0;
  Int_t extra_electron = 0;
  Int_t no_leptons = 0;
  Int_t loose_muon = 0;
  Int_t loose_electron = 0;
  Int_t no_jets = 0;
  Int_t lonesome_jet = 0;
  Int_t loose_jet = 0;
  Int_t extra_merged_jet = 0;
  Int_t met_small = 0;
  Int_t bad_dijet = 0;
  Int_t undefined_met = 0;
  Int_t no_merged_jet = 0;
  Int_t reco_lep_generation = 0;
  Int_t gen_lep_generation = 0;

  Int_t mu_tight = -1;

  Int_t merged_event = -1;
  Int_t mu_event = -1;
  Int_t el_event = -1;

  //MULTIPLE WEIGTHS
  Double_t anom_weight1;
  Double_t anom_weight2;
  Double_t anom_weight3;
  Double_t anom_weight4;
  Double_t anom_weight5;
  Double_t anom_weight6;
  Double_t anom_weight7;
  /*Double_t anom_weight8;
  Double_t anom_weight9;
  Double_t anom_weight10;
  Double_t anom_weight11;
  Double_t anom_weight12;
  Double_t anom_weight13;
  Double_t anom_weight14;
  Double_t anom_weight15;
  Double_t anom_weight16;
  Double_t anom_weight17;
  Double_t anom_weight18;
  Double_t anom_weight19;
  Double_t anom_weight20;
  / *Double_t anom_weight21;
  Double_t anom_weight22;
  Double_t anom_weight23;
  Double_t anom_weight24;
  Double_t anom_weight25;
  Double_t anom_weight26;
  Double_t anom_weight27;
  Double_t anom_weight28;
  Double_t anom_weight29;
  Double_t anom_weight30;
  Double_t anom_weight31;
  Double_t anom_weight32;
  Double_t anom_weight33;
  Double_t anom_weight34;
  Double_t anom_weight35;
  Double_t anom_weight36;
  Double_t anom_weight37;
  Double_t anom_weight38;
  Double_t anom_weight39;
  Double_t anom_weight40;
  Double_t anom_weight41;
  Double_t anom_weight42;
  Double_t anom_weight43;
  Double_t anom_weight44;
  Double_t anom_weight45;
  Double_t anom_weight46;
  Double_t anom_weight47;
  Double_t anom_weight48;
  Double_t anom_weight49;
  Double_t anom_weight50;*/

  //INTERMEDIATE STEPS VARIABLES
  Double_t leptons_in_lep_px = -99.9;
  Double_t leptons_in_lep_py = -99.9;
  Double_t leptons_in_lep_pz = -99.9;
  Double_t partons_in_lep_px = -99.9;
  Double_t partons_in_lep_py = -99.9;
  Double_t partons_in_lep_pz = -99.9;
  Double_t parton1_in_lep_px = -99.9;
  Double_t parton2_in_lep_px = -99.9;
  Double_t parton1_in_lep_py = -99.9;
  Double_t parton2_in_lep_py = -99.9;
  Double_t parton1_in_lep_pz = -99.9;
  Double_t parton2_in_lep_pz = -99.9;
  Double_t lepton1_in_lep_px = -99.9;
  Double_t lepton1_in_lep_py = -99.9;
  Double_t lepton1_in_lep_pz = -99.9;
  Double_t lepton1_dotted_x = -99.9;
  Double_t lepton1_dotted_y = -99.9;
  Double_t lepton1_dotted_z = -99.9;
  Double_t leptons_in_had_px = -99.9;
  Double_t leptons_in_had_py = -99.9;
  Double_t leptons_in_had_pz = -99.9;
  Double_t lepton1_in_had_px = -99.9;
  Double_t lepton1_in_had_py = -99.9;
  Double_t lepton1_in_had_pz = -99.9;
  Double_t lepton2_in_had_px = -99.9;
  Double_t lepton2_in_had_py = -99.9;
  Double_t lepton2_in_had_pz = -99.9;
  Double_t parton1_in_had_px = -99.9;
  Double_t parton1_in_had_py = -99.9;
  Double_t parton1_in_had_pz = -99.9;
  Double_t parton1_dotted_x = -99.9;
  Double_t parton1_dotted_y = -99.9;
  Double_t parton1_dotted_z = -99.9;
  Double_t complicated1_px = -99.9;
  Double_t complicated1_py = -99.9;
  Double_t complicated1_pz = -99.9;
  Double_t complicated2_px = -99.9;
  Double_t complicated2_py = -99.9;
  Double_t complicated2_pz = -99.9;
  Double_t lepton1WWframe_X = -99.9;
  Double_t lepton1WWframe_Y = -99.9;
  Double_t lepton1WWframe_Z = -99.9;
  Double_t lepton_sumWWframe_X = -99.9;
  Double_t lepton_sumWWframe_Y = -99.9;
  Double_t lepton_sumWWframe_Z = -99.9;
  Double_t parton1WWframe_X = -99.9;
  Double_t parton1WWframe_Y = -99.9;
  Double_t parton1WWframe_Z = -99.9;
  Double_t parton_sumWWframe_X = -99.9;
  Double_t parton_sumWWframe_Y = -99.9;
  Double_t parton_sumWWframe_Z = -99.9;

  Double_t boostWWframe_X = -99.9;
  Double_t boostWWframe_Y = -99.9;
  Double_t boostWWframe_Z = -99.9;
  Double_t boostWlep_X = -99.9;
  Double_t boostWlep_Y = -99.9;
  Double_t boostWlep_Z = -99.9;
  Double_t boostWhad_X = -99.9;
  Double_t boostWhad_Y = -99.9;
  Double_t boostWhad_Z = -99.9;

  Int_t imaginary_neutrino = -1;

  Int_t Wcount = 0;
  Int_t qcount = 0;
  
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
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig):
  ak4jetToken_(consumes<pat::JetCollection>(edm::InputTag("slimmedJets"))),
  ak8jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8chsInput"))),//edm::InputTag("slimmedJetsAK8"))),
  ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
  ak4genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"))),
  ak8genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"))),
  //rhoToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  vtxToken_(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  //triggerResultsMETFilterToken_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"))),  //"PAT"
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),//"TriggerResults", "", "HLT2"))),
  //triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("patTrigger"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),  //   selectedPatTrigger))),
  //badMuonFilterToken_(consumes<bool>(edm::InputTag("BadPFMuonFilter", ""))),
  //badChargedCandidateFilterToken_(consumes<bool>(edm::InputTag("BadChargedCandidateFilter", ""))),
  muonToken_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
  electronToken_(consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"))),
  metToken_(consumes<pat::METCollection>(edm::InputTag("slimmedMETs"))),//#multiweight
  //pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  //theSrc_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("theSrc"))),
  //pdfToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  //beamSpotToken_(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
  //conversionsToken_(consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma:reducedConversions"))),
  //jecPayloadsAK4chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4chs")),
  //jecPayloadsAK8chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8chs")),
  //jecPayloadsAK4pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4pup")),
  //jecPayloadsAK8pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8pup")),
  //jet_scale_fac_textfile_ (iConfig.getParameter<std::string> ("jet_scale_fac_textfile")),
  //ak4_textfile_ (iConfig.getParameter<std::string> ("ak4_textfile"))
  //ak8_textfile_ (iConfig.getParameter<std::string> ("ak8_textfile"))

  mgreweight1_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight1", "HLT"))),
  mgreweight2_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight2", "HLT"))),
  mgreweight3_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight3", "HLT"))),
  mgreweight4_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight4", "HLT"))),
  mgreweight5_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight5", "HLT"))),
  mgreweight6_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight6", "HLT"))),
  mgreweight7_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight7", "HLT")))/*,
  mgreweight8_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight8", "HLT"))),
  mgreweight9_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight9", "HLT"))),
  mgreweight10_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight10", "HLT"))),
  mgreweight11_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight11", "HLT"))),
  mgreweight12_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight12", "HLT"))),
  mgreweight13_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight13", "HLT"))),
  mgreweight14_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight14", "HLT"))),
  mgreweight15_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight15", "HLT"))),
  mgreweight16_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight16", "HLT"))),
  mgreweight17_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight17", "HLT"))),
  mgreweight18_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight18", "HLT"))),
  mgreweight19_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight19", "HLT"))),
  mgreweight20_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight20", "HLT")))//,*/
  /*mgreweight21_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight21", "HLT"))),
  mgreweight22_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight22", "HLT"))),
  mgreweight23_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight23", "HLT"))),
  mgreweight24_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight24", "HLT"))),
  mgreweight25_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight25", "HLT"))),
  mgreweight26_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight26", "HLT"))),
  mgreweight27_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight27", "HLT"))),
  mgreweight28_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight28", "HLT"))),
  mgreweight29_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight29", "HLT"))),
  mgreweight30_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight30", "HLT"))),
  mgreweight31_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight31", "HLT"))),
  mgreweight32_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight32", "HLT"))),
  mgreweight33_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight33", "HLT"))),
  mgreweight34_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight34", "HLT"))),
  mgreweight35_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight35", "HLT"))),
  mgreweight36_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight36", "HLT"))),
  mgreweight37_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight37", "HLT"))),
  mgreweight38_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight38", "HLT"))),
  mgreweight39_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight39", "HLT"))),
  mgreweight40_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight40", "HLT"))),
  mgreweight41_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight41", "HLT"))),
  mgreweight42_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight42", "HLT"))),
  mgreweight43_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight43", "HLT"))),
  mgreweight44_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight44", "HLT"))),
  mgreweight45_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight45", "HLT"))),
  mgreweight46_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight46", "HLT"))),
  mgreweight47_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight47", "HLT"))),
  mgreweight48_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight48", "HLT"))),
  mgreweight49_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight49", "HLT"))),
  mgreweight50_(consumes<Double_t>(edm::InputTag("allWeights", "mgreweight50", "HLT")))*/
{
  //before useresource, other file has lhapdf::initpdfset something
  
  //now do what ever initialization is needed
  usesResource("TFileService");

  edm::Service<TFileService> fileservice;
  //= fileservice->make<TH1D>("","", 200, ,);

  wp_mass = fileservice->make<TH1D>("wp_mass", "gen w plus mass", 200, 60.0, 100.0);
  wm_mass = fileservice->make<TH1D>("wm_mass", "gen w minus mass", 200, 60.0, 100.0);
  //pt_smear = fileservice->make<TH1D>("pt_smear", "pt smearing", 200, -100.0, 100.0);

  comp_hist0 = fileservice->make<TH1D>("comp_hist0", "gen pz minus method 0", 200, -300.0, 300.0);
  comp_hist1 = fileservice->make<TH1D>("comp_hist1", "gen pz minus method 1", 200, -300.0, 300.0);
  comp_hist2 = fileservice->make<TH1D>("comp_hist2", "gen pz minus method 2", 200, -300.0, 300.0);
  comp_hist3 = fileservice->make<TH1D>("comp_hist3", "gen pz minus method 3", 200, -300.0, 300.0);
  comp_hist4 = fileservice->make<TH1D>("comp_hist4", "gen pz minus method 4", 200, -300.0, 300.0);
  comp_hist5 = fileservice->make<TH1D>("comp_hist5", "gen pz minus method 5", 200, -300.0, 300.0);
  comp_hist6 = fileservice->make<TH1D>("comp_hist6", "gen pz minus method 6", 200, -300.0, 300.0);
  comp_hist7 = fileservice->make<TH1D>("comp_hist7", "gen pz minus method 7", 200, -300.0, 300.0);

  ww_mass = fileservice->make<TH1D>("ww_mass", "gen ww mass", 200, 200.0, 1500.0);
  ww_pt = fileservice->make<TH1D>("ww_pt", "gen ww pt", 200, 120.0, 1500.0);
  ww_rap = fileservice->make<TH1D>("ww_rap", "gen ww rap", 200, -2.5, 2.5);

  wp_pt = fileservice->make<TH1D>("wp_pt", "gen wplus pt", 200, 0.0, 999.0);
  wm_pt = fileservice->make<TH1D>("wm_pt", "gen wminus pt", 200, 0.0, 999.0);
  whad_pt = fileservice->make<TH1D>("whad_pt", "reco whad pt", 200, 0.0, 999.0);
  wlep_pt = fileservice->make<TH1D>("wlep_pt", "reco wlep pt", 200, 0.0, 999.0);
  neu_pt = fileservice->make<TH1D>("neu_pt", "gen neutrino pt", 200, 0.0, 500.0);
  met_pt = fileservice->make<TH1D>("met_pt", "reco MET pt", 200, 0.0, 500.0);
  elr_pt = fileservice->make<TH1D>("elr_pt", "reco electron pt", 200, 0.0, 500.0);
  mur_pt = fileservice->make<TH1D>("mur_pt", "reco muon pt", 200, 0.0, 500.0);
  elg_pt = fileservice->make<TH1D>("elg_pt", "gen electron pt", 200, 0.0, 500.0);
  mug_pt = fileservice->make<TH1D>("mug_pt", "gen muon pt", 200, 0.0, 500.0);

  whad_pt_smear = fileservice->make<TH1D>("whad_pt_smear", "whad pt smearing", 200, -100.0, 100.0);
  wlep_pt_smear = fileservice->make<TH1D>("wlep_pt_smear", "wlep pt smearing", 200, -100.0, 100.0);
  met_pt_smear = fileservice->make<TH1D>("met_pt_smear", "MET pt smearing", 200, -100.0, 100.0);
  el_pt_smear = fileservice->make<TH1D>("el_pt_smear", "electron pt smearing", 200, -10.0, 10.0);
  mu_pt_smear = fileservice->make<TH1D>("mu_pt_smear", "muon pt smearing", 200, -20.0, 20.0);
  parton0_pt_smear = fileservice->make<TH1D>("parton0_pt_smear", "parton0 pt smearing", 200, -80.0, 80.0);
  parton1_pt_smear = fileservice->make<TH1D>("parton1_pt_smear", "parton1 pt smearing", 200, -80.0, 80.0);

  whad_phi_smear = fileservice->make<TH1D>("whad_phi_smear", "whad phi smearing", 200, -0.2, 0.2);
  wlep_phi_smear = fileservice->make<TH1D>("wlep_phi_smear", "wlep phi smearing", 200, -0.4, 0.4);
  met_phi_smear = fileservice->make<TH1D>("met_phi_smear", "MET phi smearing", 200, -1.0, 1.0);
  el_phi_smear = fileservice->make<TH1D>("el_phi_smear", "electron phi smearing", 200, -0.015, 0.015);
  mu_phi_smear = fileservice->make<TH1D>("mu_phi_smear", "muon phi smearing", 200, -0.005, 0.005);
  parton0_phi_smear = fileservice->make<TH1D>("parton0_phi_smear", "parton0 phi smearing", 200, -0.2, 0.2);
  parton1_phi_smear = fileservice->make<TH1D>("parton1_phi_smear", "parton1 phi smearing", 200, -0.2, 0.2);

  whad_eta_smear = fileservice->make<TH1D>("whad_eta_smear", "whad eta smearing", 200, -1.0, 1.0);
  wlep_eta_smear = fileservice->make<TH1D>("wlep_eta_smear", "wlep eta smearing", 200, -1.5, 1.5);
  met_eta_smear = fileservice->make<TH1D>("met_eta_smear", "MET eta smearing", 200, -4.0, 4.0);
  el_eta_smear = fileservice->make<TH1D>("el_eta_smear", "electron eta smearing", 200, -0.005, 0.005);
  mu_eta_smear = fileservice->make<TH1D>("mu_eta_smear", "muon eta smearing", 200, -0.002, 0.002);
  parton0_eta_smear = fileservice->make<TH1D>("parton0_eta_smear", "parton0 eta smearing", 200, -0.2, 0.2);
  parton1_eta_smear = fileservice->make<TH1D>("parton1_eta_smear", "parton1 eta smearing", 200, -0.2, 0.2);
  
  wp_rap = fileservice->make<TH1D>("wp_rap", "gen wplus rap", 200, -3.0, +3.0);
  wm_rap = fileservice->make<TH1D>("wm_rap", "gen wminus rap", 200, -3.0, +3.0);
  whad_rap = fileservice->make<TH1D>("whad_rap", "reco whad rap", 200, -2.5, +2.5);
  wlep_rap = fileservice->make<TH1D>("wlep_rap", "reco wlep rap", 200, -3.0, +3.0);
  neu_rap = fileservice->make<TH1D>("neu_rap", "gen neutrino rap", 200, -4.0, +4.0);
  met_rap = fileservice->make<TH1D>("met_rap", "reco MET rap", 200, -4.0, +4.0);
  elr_rap = fileservice->make<TH1D>("elr_rap", "reco electron rap", 200, -2.5, +2.5);
  mur_rap = fileservice->make<TH1D>("mur_rap", "reco muon rap", 200, -2.5, +2.5);
  elg_rap = fileservice->make<TH1D>("elg_rap", "gen electron rap", 200, -2.5, +2.5);
  mug_rap = fileservice->make<TH1D>("mug_rap", "gen muon rap", 200, -2.5, +2.5);

  wplus_theta = fileservice->make<TH1D>("wplus_theta", "gen wplus theta", 200, 0.0, 3.2);
  wminus_theta = fileservice->make<TH1D>("wminus_theta", "gen wminus theta", 200, 0.0, 3.2);
  neutrino_theta = fileservice->make<TH1D>("neutrino_theta", "gen neutrino theta", 200, 0.0, 3.2);
  lepton_theta = fileservice->make<TH1D>("lepton_theta", "gen lepton theta", 200, 0.0, 3.2);
  quark_theta = fileservice->make<TH1D>("quark_theta", " gen quark theta", 200, 0.0, 3.2);
  antiquark_theta = fileservice->make<TH1D>("antiquark_theta", "gen antiquark theta", 200, 0.0, 3.2);

  wplus_phi = fileservice->make<TH1D>("wplus_phi", "gen wplus phi", 200, -3.2, 3.2);
  wminus_phi = fileservice->make<TH1D>("wminus_phi", "gen wminus phi", 200, -3.2, 3.2);
  neutrino_phi = fileservice->make<TH1D>("neutrino_phi", "gen neutrino phi", 200, -3.2, 3.2);
  lepton_phi = fileservice->make<TH1D>("lepton_phi", "gen lepton phi", 200, -3.2, 3.2);
  quark_phi = fileservice->make<TH1D>("quark_phi", "gen quark phi", 200, -3.2, 3.2);
  antiquark_phi = fileservice->make<TH1D>("antiquark_phi", "gen antiquark phi", 200, -3.2, 3.2);
 
  
  //initialization of th1d's, pdfset

  //after useresource, other file has edmservicetfileservicefs, and a bunch of fs->make<th1d>s
   

  MyTree = new TTree("MyTree", "MyTree");

  //MyTree->Branch("", &, "/D");
  MyTree->Branch("imaginary_neutrino", &imaginary_neutrino, "imaginary_neutrino/I");
  MyTree->Branch("gen_konstant", &gen_konstant, "gen_konstant/D");
  MyTree->Branch("reco_konstant", &reco_konstant, "reco_konstant/D");

  MyTree->Branch("Wcount", &Wcount, "Wcount/I");
  MyTree->Branch("qcount", &qcount, "qcount/I");

  MyTree->Branch("Trig_Names", "vector<std::string>", &Trig_Names);
  MyTree->Branch("Trig_Prescales", "vector<int>", &Trig_Prescales);
  MyTree->Branch("Trig_Pass", "vector<bool>", &Trig_Pass);

  MyTree->Branch("Trig_Accpt_Bits", &Trig_Accpt_Bits);

  MyTree->Branch("min_mass_diff", &min_mass_diff, "min_mass_diff/D");
  MyTree->Branch("min_submass_diff", &min_submass_diff, "min_submass_diff/D");
  MyTree->Branch("hadronic_gen_2body_mass", &hadronic_gen_2body_mass, "hadronic_gen_2body_mass/D");
  MyTree->Branch("leptonic_2body_mass", &leptonic_2body_mass, "leptonic_2body_mass/D");
  MyTree->Branch("leptonic_gen_2body_mass", &leptonic_gen_2body_mass, "leptonic_gen_2body_mass/D");

  MyTree->Branch("min_quark_delta", &min_quark_delta, "min_quark_delta/D");
  MyTree->Branch("min_antiquark_delta", &min_antiquark_delta, "min_antiquark_delta/D");

  MyTree->Branch("costheta1", &costheta1, "costheta1/D");
  MyTree->Branch("costheta2", &costheta2, "costheta2/D");
  MyTree->Branch("costhetastar", &costhetastar, "costhetastar/D");
  MyTree->Branch("phi1", &phi1, "phi1/D");
  MyTree->Branch("phi2", &phi2, "phi2/D");
  MyTree->Branch("phi", &phi, "phi/D");

  MyTree->Branch("other_costheta1", &other_costheta1, "other_costheta1/D");
  MyTree->Branch("other_costheta2", &other_costheta2, "other_costheta2/D");
  MyTree->Branch("other_costhetastar", &other_costhetastar, "other_costhetastar/D");
  MyTree->Branch("other_phi1", &other_phi1, "other_phi1/D");
  MyTree->Branch("other_phi2", &other_phi2, "other_phi2/D");
  MyTree->Branch("other_phi", &other_phi, "other_phi/D");

  MyTree->Branch("truth_costheta1", &truth_costheta1, "truth_costheta1/D");
  MyTree->Branch("truth_costheta2", &truth_costheta2, "truth_costheta2/D");
  MyTree->Branch("truth_costhetastar", &truth_costhetastar, "truth_costhetastar/D");
  MyTree->Branch("truth_phi1", &truth_phi1, "truth_phi1/D");
  MyTree->Branch("truth_phi2", &truth_phi2, "truth_phi2/D");
  MyTree->Branch("truth_phi", &truth_phi, "truth_phi/D");

  MyTree->Branch("false_costheta1", &false_costheta1, "false_costheta1/D");
  MyTree->Branch("false_costheta2", &false_costheta2, "false_costheta2/D");
  MyTree->Branch("false_costhetastar", &false_costhetastar, "false_costhetastar/D");
  MyTree->Branch("false_phi1", &false_phi1, "false_phi1/D");
  MyTree->Branch("false_phi2", &false_phi2, "false_phi2/D");
  MyTree->Branch("false_phi", &false_phi, "false_phi/D");


  MyTree->Branch("costheta1_gen", &costheta1_gen, "costheta1_gen/D");
  MyTree->Branch("costheta2_gen", &costheta2_gen, "costheta2_gen/D");
  MyTree->Branch("costhetastar_gen", &costhetastar_gen, "costhetastar_gen/D");
  MyTree->Branch("phi1_gen", &phi1_gen, "phi1_gen/D");
  MyTree->Branch("phi2_gen", &phi2_gen, "phi2_gen/D");
  MyTree->Branch("phi_gen", &phi_gen, "phi_gen/D");

  MyTree->Branch("Tau21", &Tau21, "Tau21/D");
  MyTree->Branch("Tau32", &Tau32, "Tau32/D");
  MyTree->Branch("PTau21", &PTau21, "PTau21/D");
  MyTree->Branch("PTau32", &PTau32, "PTau32/D");
  
  
  //PARTICLES START HERE
  //MyTree->Branch("", &, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  
  MyTree->Branch("Wpls_gen", &Wpls_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("Wmns_gen", &Wmns_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("Zntrl_gen", &Zntrl_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("qrk_gen", &qrk_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("antiqrk_gen", &antiqrk_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("lptn_gen", &lptn_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("ntrno_gen", &ntrno_gen, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");

  MyTree->Branch("Whad_reco", &Whad_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("Wlep_reco", &Wlep_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("muon_reco", &muon_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("electrn_reco", &electrn_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("lepton_reco", &lepton_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("jt0", &jt0, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("jt1", &jt1, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("subjt0", &subjt0, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("subjt1", &subjt1, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");

  MyTree->Branch("ntrno_reco", &ntrno_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("other_ntrno_reco", &other_ntrno_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("truth_ntrno_reco", &truth_ntrno_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");
  MyTree->Branch("false_ntrno_reco", &false_ntrno_reco, "Px/D:Py:Pz:E:M:Pt:Eta:Theta:Phi:Y");


  MyTree->Branch("lepton_deltaR", &lepton_deltaR, "lepton_deltaR/D");
  MyTree->Branch("neutrino_deltaR", &neutrino_deltaR, "neutrino_deltaR/D");
  MyTree->Branch("other_neutrino_deltaR", &other_neutrino_deltaR, "other_neutrino_deltaR/D");
  MyTree->Branch("truth_neutrino_deltaR", &truth_neutrino_deltaR, "truth_neutrino_deltaR/D");
  MyTree->Branch("false_neutrino_deltaR", &false_neutrino_deltaR, "false_neutrino_deltaR/D");
  MyTree->Branch("jet0_deltaR", &jet0_deltaR, "jet0_deltaR/D");
  MyTree->Branch("jet1_deltaR", &jet1_deltaR, "jet1_deltaR/D");
  MyTree->Branch("subjet0_deltaR", &subjet0_deltaR, "subjet0_deltaR/D");
  MyTree->Branch("subjet1_deltaR", &subjet1_deltaR, "subjet1_deltaR/D");

  MyTree->Branch("lepton_deltaPhi", &lepton_deltaPhi, "lepton_deltaPhi/D");
  MyTree->Branch("neutrino_deltaPhi", &neutrino_deltaPhi, "neutrino_deltaPhi/D");
  MyTree->Branch("other_neutrino_deltaPhi", &other_neutrino_deltaPhi, "other_neutrino_deltaPhi/D");
  MyTree->Branch("neutrino_deltaPt", &neutrino_deltaPt, "neutrino_deltaPt/D");

  MyTree->Branch("extra_muon", &extra_muon, "extra_muon/I");
  MyTree->Branch("extra_electron", &extra_electron, "extra_electron/I");
  MyTree->Branch("no_leptons", &no_leptons, "no_leptons/I");
  MyTree->Branch("loose_muon", &loose_muon, "loose_muon/I");
  MyTree->Branch("loose_electron", &loose_electron, "loose_electron/I");
  MyTree->Branch("no_jets", &no_jets, "no_jets/I");
  MyTree->Branch("lonesome_jet", &lonesome_jet, "lonesome_jet/I");
  MyTree->Branch("loose_jet", &loose_jet, "loose_jet/I");
  MyTree->Branch("extra_merged_jet", &extra_merged_jet, "extra_merged_jet/I");
  MyTree->Branch("met_small", &met_small, "met_small/I");
  MyTree->Branch("bad_dijet", &bad_dijet, "bad_dijet/I");
  MyTree->Branch("undefined_met", &undefined_met, "undefined_met/I");
  MyTree->Branch("no_merged_jet", &no_merged_jet, "no_merged_jet/I");
  MyTree->Branch("mu_tight", &mu_tight, "mu_tight/I");
  MyTree->Branch("merged_event", &merged_event, "merged_event/I");
  MyTree->Branch("mu_event", &mu_event, "mu_event/I");
  MyTree->Branch("el_event", &el_event, "el_event/I");

  MyTree->Branch("anom_weight1", &anom_weight1, "anom_weight1/D");
  MyTree->Branch("anom_weight2", &anom_weight2, "anom_weight2/D");
  MyTree->Branch("anom_weight3", &anom_weight3, "anom_weight3/D");
  MyTree->Branch("anom_weight4", &anom_weight4, "anom_weight4/D");
  MyTree->Branch("anom_weight5", &anom_weight5, "anom_weight5/D");
  MyTree->Branch("anom_weight6", &anom_weight6, "anom_weight6/D");
  MyTree->Branch("anom_weight7", &anom_weight7, "anom_weight7/D");
  /*MyTree->Branch("anom_weight8", &anom_weight8, "anom_weight8/D");
  MyTree->Branch("anom_weight9", &anom_weight9, "anom_weight9/D");
  MyTree->Branch("anom_weight10", &anom_weight10, "anom_weight10/D");
  MyTree->Branch("anom_weight11", &anom_weight11, "anom_weight11/D");
  MyTree->Branch("anom_weight12", &anom_weight12, "anom_weight12/D");
  MyTree->Branch("anom_weight13", &anom_weight13, "anom_weight13/D");
  MyTree->Branch("anom_weight14", &anom_weight14, "anom_weight14/D");
  MyTree->Branch("anom_weight15", &anom_weight15, "anom_weight15/D");
  MyTree->Branch("anom_weight16", &anom_weight16, "anom_weight16/D");
  MyTree->Branch("anom_weight17", &anom_weight17, "anom_weight17/D");
  MyTree->Branch("anom_weight18", &anom_weight18, "anom_weight18/D");
  MyTree->Branch("anom_weight19", &anom_weight19, "anom_weight19/D");
  MyTree->Branch("anom_weight20", &anom_weight20, "anom_weight20/D");
  / *MyTree->Branch("anom_weight21", &anom_weight21, "anom_weight21/D");
  MyTree->Branch("anom_weight22", &anom_weight22, "anom_weight22/D");
  MyTree->Branch("anom_weight23", &anom_weight23, "anom_weight23/D");
  MyTree->Branch("anom_weight24", &anom_weight24, "anom_weight24/D");
  MyTree->Branch("anom_weight25", &anom_weight25, "anom_weight25/D");
  MyTree->Branch("anom_weight26", &anom_weight26, "anom_weight26/D");
  MyTree->Branch("anom_weight27", &anom_weight27, "anom_weight27/D");
  MyTree->Branch("anom_weight28", &anom_weight28, "anom_weight28/D");
  MyTree->Branch("anom_weight29", &anom_weight29, "anom_weight29/D");
  MyTree->Branch("anom_weight30", &anom_weight30, "anom_weight30/D");
  MyTree->Branch("anom_weight31", &anom_weight31, "anom_weight31/D");
  MyTree->Branch("anom_weight32", &anom_weight32, "anom_weight32/D");
  MyTree->Branch("anom_weight33", &anom_weight33, "anom_weight33/D");
  MyTree->Branch("anom_weight34", &anom_weight34, "anom_weight34/D");
  MyTree->Branch("anom_weight35", &anom_weight35, "anom_weight35/D");
  MyTree->Branch("anom_weight36", &anom_weight36, "anom_weight36/D");
  MyTree->Branch("anom_weight37", &anom_weight37, "anom_weight37/D");
  MyTree->Branch("anom_weight38", &anom_weight38, "anom_weight38/D");
  MyTree->Branch("anom_weight39", &anom_weight39, "anom_weight39/D");
  MyTree->Branch("anom_weight40", &anom_weight40, "anom_weight40/D");
  MyTree->Branch("anom_weight41", &anom_weight41, "anom_weight41/D");
  MyTree->Branch("anom_weight42", &anom_weight42, "anom_weight42/D");
  MyTree->Branch("anom_weight43", &anom_weight43, "anom_weight43/D");
  MyTree->Branch("anom_weight44", &anom_weight44, "anom_weight44/D");
  MyTree->Branch("anom_weight45", &anom_weight45, "anom_weight45/D");
  MyTree->Branch("anom_weight46", &anom_weight46, "anom_weight46/D");
  MyTree->Branch("anom_weight47", &anom_weight47, "anom_weight47/D");
  MyTree->Branch("anom_weight48", &anom_weight48, "anom_weight48/D");
  MyTree->Branch("anom_weight49", &anom_weight49, "anom_weight49/D");
  MyTree->Branch("anom_weight50", &anom_weight50, "anom_weight50/D");*/
  /*
  MyTree->Branch("leptons_in_lep_px", &leptons_in_lep_px, "leptons_in_lep_px/D");
  MyTree->Branch("leptons_in_lep_py", &leptons_in_lep_py, "leptons_in_lep_py/D");
  MyTree->Branch("leptons_in_lep_pz", &leptons_in_lep_pz, "leptons_in_lep_pz/D");
  MyTree->Branch("partons_in_lep_px", &partons_in_lep_px, "partons_in_lep_px/D");
  MyTree->Branch("partons_in_lep_py", &partons_in_lep_py, "partons_in_lep_py/D");
  MyTree->Branch("partons_in_lep_pz", &partons_in_lep_pz, "partons_in_lep_pz/D");
  MyTree->Branch("parton1_in_lep_px", &parton1_in_lep_px, "parton1_in_lep_px/D");
  MyTree->Branch("parton2_in_lep_px", &parton2_in_lep_px, "parton2_in_lep_px/D");
  MyTree->Branch("parton1_in_lep_py", &parton1_in_lep_py, "parton1_in_lep_py/D");
  MyTree->Branch("parton2_in_lep_py", &parton2_in_lep_py, "parton2_in_lep_py/D");
  MyTree->Branch("parton1_in_lep_pz", &parton1_in_lep_pz, "parton1_in_lep_pz/D");
  MyTree->Branch("parton2_in_lep_pz", &parton2_in_lep_pz, "parton2_in_lep_pz/D");
  MyTree->Branch("lepton1_in_lep_px", &lepton1_in_lep_px, "lepton1_in_lep_px/D");
  MyTree->Branch("lepton1_in_lep_py", &lepton1_in_lep_py, "lepton1_in_lep_py/D");
  MyTree->Branch("lepton1_in_lep_pz", &lepton1_in_lep_pz, "lepton1_in_lep_pz/D");
  MyTree->Branch("lepton1_dotted_x", &lepton1_dotted_x, "lepton1_dotted_x/D");
  MyTree->Branch("lepton1_dotted_y", &lepton1_dotted_y, "lepton1_dotted_y/D");
  MyTree->Branch("lepton1_dotted_z", &lepton1_dotted_z, "lepton1_dotted_z/D");
  MyTree->Branch("leptons_in_had_px", &leptons_in_had_px, "leptons_in_had_px/D");
  MyTree->Branch("leptons_in_had_py", &leptons_in_had_py, "leptons_in_had_py/D");
  MyTree->Branch("leptons_in_had_pz", &leptons_in_had_pz, "leptons_in_had_pz/D");
  MyTree->Branch("lepton1_in_had_px", &lepton1_in_had_px, "lepton1_in_had_px/D");
  MyTree->Branch("lepton1_in_had_py", &lepton1_in_had_py, "lepton1_in_had_py/D");
  MyTree->Branch("lepton1_in_had_pz", &lepton1_in_had_pz, "lepton1_in_had_pz/D");
  MyTree->Branch("lepton2_in_had_px", &lepton2_in_had_px, "lepton2_in_had_px/D");
  MyTree->Branch("lepton2_in_had_py", &lepton2_in_had_py, "lepton2_in_had_py/D");
  MyTree->Branch("lepton2_in_had_pz", &lepton2_in_had_pz, "lepton2_in_had_pz/D");
  MyTree->Branch("parton1_in_had_px", &parton1_in_had_px, "parton1_in_had_px/D");
  MyTree->Branch("parton1_in_had_py", &parton1_in_had_py, "parton1_in_had_py/D");
  MyTree->Branch("parton1_in_had_pz", &parton1_in_had_pz, "parton1_in_had_pz/D");
  MyTree->Branch("parton1_dotted_x", &parton1_dotted_x, "parton1_dotted_x/D");
  MyTree->Branch("parton1_dotted_y", &parton1_dotted_y, "parton1_dotted_y/D");
  MyTree->Branch("parton1_dotted_z", &parton1_dotted_z, "parton1_dotted_z/D");
  MyTree->Branch("complicated1_px", &complicated1_px, "complicated1_px/D");
  MyTree->Branch("complicated1_py", &complicated1_py, "complicated1_py/D");
  MyTree->Branch("complicated1_pz", &complicated1_pz, "complicated1_pz/D");
  MyTree->Branch("complicated2_px", &complicated2_px, "complicated2_px/D");
  MyTree->Branch("complicated2_py", &complicated2_py, "complicated2_py/D");
  MyTree->Branch("complicated2_pz", &complicated2_pz, "complicated2_pz/D");
  
  MyTree->Branch("lepton_sumWWframe_X", &lepton_sumWWframe_X, "lepton_sumWWframe_X/D");
  MyTree->Branch("lepton_sumWWframe_Y", &lepton_sumWWframe_Y, "lepton_sumWWframe_Y/D");
  MyTree->Branch("lepton_sumWWframe_Z", &lepton_sumWWframe_Z, "lepton_sumWWframe_Z/D");
  
  MyTree->Branch("parton_sumWWframe_X", &parton_sumWWframe_X, "parton_sumWWframe_X/D");
  MyTree->Branch("parton_sumWWframe_Y", &parton_sumWWframe_Y, "parton_sumWWframe_Y/D");
  MyTree->Branch("parton_sumWWframe_Z", &parton_sumWWframe_Z, "parton_sumWWframe_Z/D");
  
  MyTree->Branch("lepton1WWframe_X", &lepton1WWframe_X, "lepton1WWframe_X/D");
  MyTree->Branch("lepton1WWframe_Y", &lepton1WWframe_Y, "lepton1WWframe_Y/D");
  MyTree->Branch("lepton1WWframe_Z", &lepton1WWframe_Z, "lepton1WWframe_Z/D");
  
  MyTree->Branch("parton1WWframe_X", &parton1WWframe_X, "parton1WWframe_X/D");
  MyTree->Branch("parton1WWframe_Y", &parton1WWframe_Y, "parton1WWframe_Y/D");
  MyTree->Branch("parton1WWframe_Z", &parton1WWframe_Z, "parton1WWframe_Z/D");

  MyTree->Branch("boostWWframe_X", &boostWWframe_X, "boostWWframe_X/D");
  MyTree->Branch("boostWWframe_Y", &boostWWframe_Y, "boostWWframe_Y/D");
  MyTree->Branch("boostWWframe_Z", &boostWWframe_Z, "boostWWframe_Z/D");
  MyTree->Branch("boostWlep_X", &boostWlep_X, "boostWlep_X/D");
  MyTree->Branch("boostWlep_Y", &boostWlep_Y, "boostWlep_Y/D");
  MyTree->Branch("boostWlep_Z", &boostWlep_Z, "boostWlep_Z/D");
  MyTree->Branch("boostWhad_X", &boostWhad_X, "boostWhad_X/D");
  MyTree->Branch("boostWhad_Y", &boostWhad_Y, "boostWhad_Y/D");
  MyTree->Branch("boostWhad_Z", &boostWhad_Z, "boostWhad_Z/D");*/

  if (talking) outfile.open("special_output.txt");
}


ntuplizer::~ntuplizer()
{
  if (talking)
    {
      outfile << "extra_muon_count " << extra_muon_count << std::endl;
      outfile << "extra_electron_count " << extra_electron_count << std::endl; 
      outfile << "no_lepton_found_count " << no_lepton_found_count << std::endl; 
      outfile << "loose_muon_count " << loose_muon_count << std::endl; 
      outfile << "loose_electron_count " << loose_electron_count << std::endl; 
      outfile << "no_jet_count " << no_jet_count << std::endl; 
      outfile << "merged_combo_count " << merged_combo_count << std::endl; 
      outfile << "extra_merged_combo_count " << extra_merged_combo_count << std::endl; 
      outfile << "missing_et_cutoff_count " << missing_et_cutoff_count << std::endl;

      outfile << "nonmerged_count " << nonmerged_count << std::endl;
      outfile << "inverted_count " << inverted_count << std::endl;
      outfile << "neither_count " << neither_count << std::endl;

      outfile << "skip_counter " << skip_counter << std::endl;
      outfile << "backup_skip_counter " << backup_skip_counter << std::endl;
  
      outfile << "here's the sum of counters: merged_counter " << merged_counter << std::endl;
      outfile << "here's the number of merged events: merged_event_count " << merged_event_count << std::endl;
      // do anything here that needs to be done at desctruction time
      // (e.g. close files, deallocate resources etc.)
      outfile.close();
    }
   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  //std::cout << "\n\nbegin event.\n";
  Bool_t bad_event = false;
//#if 0
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  //using namespace LHAPDF;

  reco_lep_generation = 0;
  gen_lep_generation = 0;
  merged_event = -1;
  if (talking) outfile << "===================" << std::endl;
  good_events++;

    
  
  //MARKER HIGH LEVEL TRIGGER
  //  _   _ _       _       _         _   _____                   
  // | | | (_) __ _| |__   | | __   _| | |_   _| __ __ _ _ __ ___ 
  // | |_| | |/ _` | '_ \  | | \ \ / / |   | || '__/ _` | '__/ __|
  // |  _  | | (_| | | | | | |__\ V /| |   | || | | (_| | |  \__ \ |
  // |_| |_|_|\__, |_| |_| |_____\_/ |_|   |_||_|  \__, |_|  |___/
  //          |___/                                |___/          
  //

  if (bad_event) return; //{std::cout << "bad event reached HLT\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
   
  int match_count = 0;
  
  edm::Handle<edm::TriggerResults> triggerBits;
  //edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;


  iEvent.getByToken(triggerBits_, triggerBits);
  //iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  //LATER ON
  //const edm::TriggerNames &list_of_names = iEvent.triggerNames(*triggerBits);
  //int trigger_count = 0;
  //const int triggerBitsSize = triggerBits->size();

  vector<string> triggers_list;

  // HT 
  //triggers_list.push_back("HLT_PFHT300_v");
  //triggers_list.push_back("HLT_PFHT350_v");
  //triggers_list.push_back("HLT_PFHT400_v");
  //triggers_list.push_back("HLT_PFHT475_v");
  //triggers_list.push_back("HLT_PFHT600_v");
  //triggers_list.push_back("HLT_PFHT650_v");
  //triggers_list.push_back("HLT_PFHT800_v");
  //triggers_list.push_back("HLT_PFHT900_v");
  //triggers_list.push_back("HLT_PFHT650_WideJetMJJ900"); //HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v6
  //triggers_list.push_back("HLT_PFHT650_WideJetMJJ950"); //HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v6

  // Single jet
  //triggers_list.push_back("HLT_CaloJet500_NoJetID_v");
  //triggers_list.push_back("HLT_PFJet320_v");
  //triggers_list.push_back("HLT_PFJet400_v");
  //triggers_list.push_back("HLT_PFJet450_v");
  //triggers_list.push_back("HLT_PFJet500_v");
  //triggers_list.push_back("HLT_AK8PFJet450_v");
  //triggers_list.push_back("HLT_AK8PFJet500_v");

  // Substructure
  //triggers_list.push_back("HLT_AK8PFJet360_TrimMass30_v");
  //triggers_list.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  //triggers_list.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");

  // Substructure + b-tag
  // triggers_list.push_back("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  //triggers_list.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");

  // Muon
  //triggers_list.push_back("HLT_Mu45_eta2p1_v");
  triggers_list.push_back("HLT_Mu50_v");//one used
  //triggers_list.push_back("HLT_Mu55_v");
  triggers_list.push_back("HLT_TkMu50_v");
  //triggers_list.push_back("HLT_IsoMu22_eta2p1_v");
  //triggers_list.push_back("HLT_IsoMu24_v");
  //triggers_list.push_back("HLT_IsoMu27_v");

  // Muon + jet
  //triggers_list.push_back("HLT_Mu30_eta2p1_PFJet150_PFJet50_v");
  //triggers_list.push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v");

  // Electron
  //triggers_list.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
  //triggers_list.push_back("HLT_Ele35_WPLoose_Gsf_v");
  triggers_list.push_back("HLT_Ele27_WPLoose_Gsf_v");
  triggers_list.push_back("HLT_Ele30_WPTight_Gsf_v");
  triggers_list.push_back("HLT_Ele45_WPLoose_Gsf_v");
  triggers_list.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v");

  triggers_list.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  triggers_list.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");

  // Electron + jet
  //triggers_list.push_back("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
  //triggers_list.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_v");
  //triggers_list.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v");
  //triggers_list.push_back("HLTriggerFirstPa");

  //junk test
  const int num_of_triggers = triggers_list.size();
  //std::cout << "number of triggers is " << num_of_triggers << std::endl;

  //same thing section
  std::bitset<38> high_lvl_trig_bit; //converted to string for branch, Trig_Accept_Bits
  vector<bool> accepted_trigger; //seems redundant with Trig_Pass

  Trig_Names->clear();//lost
  Trig_Prescales->clear();
  Trig_Pass->clear();
  
  const edm::TriggerNames &list_of_names = iEvent.triggerNames(*triggerBits);
  int trigger_count = 0;
  const int triggerBitsSize = triggerBits->size();
  //std::cout << "trigger bit size is " << triggerBitsSize << std::endl;
  bool some_pass = false;

  for (int i = 0; i < num_of_triggers; i++)
    {
     // std::cout << "trying to find " << triggers_list[i] << std::endl;
      //bool found_trig  = false;
      bool passed = false;
      int num_of_prescales = 0;
      string name_in_list = "";
      int accepted = -1;
      
      for (int j = 0; j < triggerBitsSize; j++)
	{
	  //std::cout << "======================" << std::endl << "starting iteration of for loop with i = " << i << " and j = " << j << std::endl;
	  //std::cout << "n is equal to " << n << std::endl;
	  name_in_list = list_of_names.triggerName(j);
	  //std::cout << "name_in_list is " << name_in_list << std::endl;
	  //std::cout << "current index in truggers_list is " << triggers_list[i] << std::endl;
	  std::size_t finds = name_in_list.find(triggers_list[i]);//LOOK FOR MATCH
	  //std::cout << "we just defined finds.\n";
	  //std::cout << "ending iteration of for loop with i = " << i << " and j = " << j << std::endl << "====================" << std::endl;
	  
	  if (finds != std::string::npos)//IF MATCH WAS FOUND
	    {
	      match_count++;
	      //std::cout << "found a match\n";
	      //found_trig = true;
	      accepted = triggerBits->accept(j);//PASSED OR NOT?
	      //std::cout << "accepted is " << accepted << std::endl;
	      //std::cout << "finds is " << finds << std::endl;
	      if (accepted == 1) passed = true;
	      if (passed) some_pass = true;
	      //if (accepted > 1) {std::cout << "accepted greater than one?\n"; exit(42);}
	      //std::cout << "now to set num_of_prescales" << std::endl;
	      num_of_prescales = triggerPrescales->getPrescaleForIndex(j);
	      //std::cout << "just set num_of_prescales" << std::endl;
	      //string result = passed ? "pass" : "fail";

	      //std::cout << "Found a Trigger " << name_in_list << " number(?) prescale is " << num_of_prescales << " : " << result << std::endl;
	      break;
	    }
	}// end loop over all triggers in event

      accepted_trigger.push_back(passed);
	      
      Trig_Names->push_back(name_in_list);
      Trig_Prescales->push_back(num_of_prescales);
      Trig_Pass->push_back(passed);
	      
      if (passed) high_lvl_trig_bit[trigger_count] = 1;  
      trigger_count++;
    }// end loop over list of triggers to save in tree

 
  //int accepted_trig_size = accepted_trigger.size();

  //std::cout << "accepted_trigger size is " << accepted_trig_size << std::endl;

  //for (int k = 0; k < accepted_trig_size; k++) std::cout << "accepted trigger is " << accepted_trigger[accepted_trig_size - 1 - k] << std::endl;

  Trig_Accpt_Bits = high_lvl_trig_bit.to_string();

  if (talking) outfile << "number of matches is " << match_count << std::endl;

  if (!some_pass) return; //{std::cout << "no triggers.\n"; return;}
  
  //std::cout << std::endl << "Trig Accpt Bits is " << Trig_Accpt_Bits << std::endl;


  // MARKER GEN PARTICLES
  //    .d8888b.  8888888888 888b    888     8888888b.                   888    d8b          888                   
  //   d88P  Y88b 888        8888b   888     888   Y88b                  888    Y8P          888                   
  //   888    888 888        88888b  888     888    888                  888                 888                   
  //   888        8888888    888Y88b 888     888   d88P  8888b.  888d888 888888 888  .d8888b 888  .d88b.  .d8888b  
  //   888  88888 888        888 Y88b888     8888888P"      "88b 888P"   888    888 d88P"    888 d8P  Y8b 88K      
  //   888    888 888        888  Y88888     888        .d888888 888     888    888 888      888 88888888 "Y8888b. 
  //   Y88b  d88P 888        888   Y8888     888        888  888 888     Y88b.  888 Y88b.    888 Y8b.          X88 
  //    "Y8888P88 8888888888 888    Y888     888        "Y888888 888      "Y888 888  "Y8888P 888  "Y8888   88888P' 
  //                                                                                                               

  if (bad_event) return;
  TLorentzVector Wplus_gen;
  TLorentzVector Wminus_gen;
  TLorentzVector Zneutral_gen;
  TLorentzVector quark_gen;
  TLorentzVector antiquark_gen;
  TLorentzVector lepton_gen;
  TLorentzVector neutrino_gen;

  Wplus_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  Wminus_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  Zneutral_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  quark_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  antiquark_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  lepton_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  neutrino_gen.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);

  //lepton.SetPxPyPzE(99.0, 99.0, 99.0, 99.0);

  //Int_t charge = 0;
  //charge = charge;//keep asinine compiler from throwing a fit if I later set it but don't use it.
  
  if (!iEvent.isRealData())// and runGenLoop_) //if event is mc, not data
    {
      Handle<edm::View<reco::GenParticle>> genpart; //create pointer to .... particle array
      iEvent.getByToken(prunedGenToken_, genpart);  //connect said pointer to event

      //std::cout << "SEARCHING: About to enter for loop, looking for status of particles.\n";

      bool lepton_found = false;
      bool neutrino_found = false;

      Wcount = 0;
      qcount = 0;
      //lepcount = 0;
      //neucount = 0;
      //int wflag = -10;
   
      for (size_t i = 0; i < genpart->size(); i++) //i is particle number, in a single event, sets, and resets, a bunch of quantities for the event (entire analyze function applied event-by-event)
	{
	  int id        = (*genpart)[i].pdgId(); //genpart[i] a particle
	  int status    = (*genpart)[i].status();
	  Double_t px     = (*genpart)[i].px();
	  Double_t py     = (*genpart)[i].py();
	  Double_t pz     = (*genpart)[i].pz();
	  Double_t e      = (*genpart)[i].energy();
	  int ndau      = (*genpart)[i].numberOfDaughters();
	  //const reco::Candidate*  mother_    = (*genpart)[i].mother();
	  //int mother = (*mother_).pdgId();
	  /*
	    std::string statsp = " ";
	    std::string daughtsp = " ";
	    std::string idsp = " ";
	    std::string xspace = " ";
	    std::string yspace = " ";
	    std::string zspace = " ";
	  */

	
	  
	  //EWK Gauge bosons

	  if (id ==  24)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (status != 2 && talking) outfile << "status of W is " << status << std::endl;
	      if (ndau != 2 && talking) outfile << "number of daughtors for the W is " << ndau << std::endl;
	      Wcount++;
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      Wplus_gen.SetPxPyPzE(px, py, pz, e);
	      //wflag = +1;
	      //std::cout << "SEARCHING: status of Wplus  is" << statsp << status << " number of Wplus  daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }
	  
	  if (id == -24)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (status != 2 && talking) outfile << "status of W is " << status << std::endl;
	      if (ndau != 2 && talking) outfile << "number of daughtors for the W is " << ndau << std::endl;
	      Wcount++;
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      Wminus_gen.SetPxPyPzE(px, py, pz, e);
	      //wflag = -1;
	      //std::cout << "SEARCHING: status of Wminus is" << statsp << status << " number of Wminus daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }
	    
	  if (id ==  23)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      Zneutral_gen.SetPxPyPzE(px, py, pz, e);
	      //wflag = 0;
	      //std::cout << "SEARCHING: status of Zneu   is" << statsp << status << " number of Zneu   daughters is" << daughtsp << ndau << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  //decay product fermions
	  
	  if (id == 11 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "lepton has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = 1;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }


	  if (id == -11 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "lepton has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = -1;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  if (id == 13 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "lepton has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = 2;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  if (id == -13 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "lepton has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = -2;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  if (id == 15 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "tau has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = 3;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  if (id == -15 && !lepton_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "tau has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      lepton_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id < 0) charge =  1;
	      //if (id > 0) charge = -1;
	      lepton_found = true;
	      gen_lep_generation = -3;
	      //if (wflag == -1) version = negative;
	      //if (wflag == +1) version = positive;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  lepton is" << statsp << status << " number of  lepton daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  
	  if (id == 12 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = -1;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  
	  if (id == -12 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = 1;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  
	  if (id == 14 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = -2;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	  
	  if (id == -14 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = 2;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	    
	  if (id == 16 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = -3;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }

	    
	  if (id == -16 && !neutrino_found)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      if (ndau != 0 && talking) outfile << "neutrino has " << ndau << " daughters.\n";
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace= " ";
	      //if (pz >= 0) zspace = " +";
	      neutrino_gen.SetPxPyPzE(px, py, pz, e);
	      //if (id > 0) charge =  1;
	      //if (id < 0) charge = -1;
	      neutrino_found = true;
	      //gen_lep_generation = 3;
	      //std::cout << "SEARCHING: status of neutrno is" << statsp << status << " number of neutrno daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }
	    

	  if (id >=  1 && id <=  6 && status == 23)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      qcount++;
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      quark_gen.SetPxPyPzE(px, py, pz, e);
	      //if (wflag == -1) version = positive;
	      //if (wflag == +1) version = negative;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of   quark is" << statsp << status << " number of   quark daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }
	      
	  if (id <= -1 && id >= -6 && status == 23)
	    { //std::cout << "id is " << id << " and mother is " << mother << std::endl;
	      qcount++;
	      //if (status < 10) statsp = "  ";
	      //if (status >= 10) statsp = " ";
	      //if (ndau < 10) daughtsp = "  ";
	      //if (ndau >= 10) daughtsp = " ";
	      //if (id <= -10) idsp = " ";
	      //if (id > -10 && id < 0) idsp = "  ";
	      //if (id < 10 && id >= 0) idsp = "  +";
	      //if (id >= 10) idsp = " +";
	      //if (px < 0) xspace = " ";
	      //if (px >= 0) xspace = " +";
	      //if (py < 0) yspace = " ";
	      //if (py >= 0) yspace = " +";
	      //if (pz < 0) zspace = " ";
	      //if (pz >= 0) zspace = " +";
	      antiquark_gen.SetPxPyPzE(px, py, pz, e);
	      //if (wflag == -1) version = positive;
	      //if (wflag == +1) version = negative;
	      //if (wflag == 0) version = neutral;
	      //std::cout << "SEARCHING: status of  aquark is" << statsp << status << " number of  aquark daughters is" << daughtsp << ndau << " id is" << idsp << id << " px is" << xspace << px << " py is" << yspace << py << " pz is" << zspace << pz << " e is " << e << std::endl;
	    }
	
	} // end genParticle loop
      
      //std::cout << "SEARCHING: Just exited for loop, looking for status of particles.\n";
     
    
    }

       
      
  // MARKER VERTICIES
  //
  //
  //  __     _______ ____ _____ ___ ____ ___ _____ ____  
  //  \ \   / / ____|  _ \_   _|_ _/ ___|_ _| ____/ ___| 
  //   \ \ / /|  _| | |_) || |  | | |    | ||  _| \___ \ |
  //    \ V / | |___|  _ < | |  | | |___ | || |___ ___) |
  //     \_/  |_____|_| \_\|_| |___\____|___|_____|____/ 
  //
  //

  if (bad_event) return; //{std::cout << "bad event reached verticies.\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
  
  edm::Handle<std::vector<reco::Vertex>> verticies;
  iEvent.getByToken(vtxToken_, verticies);
  /////int num_verticies = verticies->size();
  if (verticies->empty()) return; //{std::cout << "no primary vertex found.\n"; return;} // skip the event if no PV found
  const reco::Vertex &primary_vertex = verticies->front();  // save PV for tight muon ID

  int num_good_verticies = 0;//all from here down just c-out statements
  bool fake_vertex = true;
  for (std::vector<reco::Vertex>::const_iterator vertex_iterator = verticies->begin(); vertex_iterator != verticies->end(); ++vertex_iterator)
    {
      fake_vertex = (vertex_iterator->chi2() == 0 && vertex_iterator->ndof() == 0);   //// bool fake_vertex = vertex_iterator->fake_vertex();  // for AOD
      if (!fake_vertex && vertex_iterator->ndof() >= 4.0 && vertex_iterator->position().Rho() <= 2.0 && fabs(vertex_iterator->position().Z()) <= 24.0) num_good_verticies++;
    }
  //std::cout << "number of verticies " << num_verticies << " number of good verticies " << num_good_verticies << std::endl;


  // MARKER RECO LEPTONS
  // ____  _____ ____ ___    _     _____ ____ _____ ___  _   _ ____  
  //|  _ \| ____/ ___/ _ \  | |   | ____|  _ \_   _/ _ \| \ | / ___| 
  //| |_) |  _|| |  | | | | | |   |  _| | |_) || || | | |  \| \___ \ |
  //|  _ <| |__| |__| |_| | | |___| |___|  __/ | || |_| | |\  |___) |
  //|_| \_\_____\____\___/  |_____|_____|_|    |_| \___/|_| \_|____/ 
  //

  if (bad_event) return; //{std::cout << "bad event reached reco leptons"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
  
  Int_t electron_value_set = 0;
  Int_t muon_value_set = 0;
  
  Int_t electron_out_of_bounds = 0;
  Int_t muon_out_of_bounds = 0;
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  TLorentzVector muon;
  TLorentzVector electron;
  TLorentzVector lepton;

  
  muon.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
  electron.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
  lepton.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

  el_event = -1;
  mu_event = -1;
 

  Int_t el_charge = 0;
  Int_t mu_charge = 0;
  Int_t counter = -1;

  Int_t mu_mark = -1;
  Int_t el_mark = -1;

  Bool_t lepton_found = false;
  
  
  //int met_counter = 0;
  //int non_met_counter = 0;

  extra_muon = 0;
  extra_electron = 0;
  loose_muon = 0;
  loose_electron = 0;
 
  no_leptons = 0;
 
 
  const Double_t tight_muon_cutoff_pt = 53.0; //was 25.0
  const Double_t tight_electron_cutoff_pt = 50.0; //was 30.0
  const Double_t loose_muon_cutoff_pt = 20.0; //was ?? channel
  const Double_t loose_electron_cutoff_pt = 20.0; //35.0; //was 20.0 channel
  const Double_t muon_rapidity = 2.4; //was 2.1 (2.5 loose)
  const Double_t electron_rapidity = 2.5; //was 2.5 
  const Double_t mu_channel_missing_et_cutoff = 40.0; //was 25.0
  const Double_t el_channel_missing_et_cutoff = 80.0; //was 25.0
 
    
  Double_t loose_lepton_cutoff_pt = -1.0;
  Double_t missing_et_cutoff = -1.0;

  //electron eta cannot be between 1.44 and 1.57 (is this absolute value?)
  //std::cout << "COUNTER\n";

  for (const pat::Muon &mu_index : *muons) //tight muons
    {
      //if (bad_event) break;
      counter++;
      if (mu_index.pt() > tight_muon_cutoff_pt && fabs(mu_index.eta()) < muon_rapidity)
	{
	  if (lepton_found) {bad_event = true; extra_muon = 1; extra_muon_count++;}// std::cout << "BAD EVENT: more than one muon.\n";}
	  else
	    {
	      lepton_found = true;
	      muon.SetPxPyPzE(mu_index.px(), mu_index.py(), mu_index.pz(), sqrt(mu_index.px()*mu_index.px() + mu_index.py()*mu_index.py() + mu_index.pz()*mu_index.pz() + mu_mass*mu_mass));
	      muon_value_set = 1;
	      mu_event = 1;
	      el_event = 0;
	      mu_mark = counter;
	      mu_charge = mu_index.charge();
	      if (mu_index.isTightMuon(primary_vertex)) mu_tight = 1;
	      else mu_tight = 0;
	    }
	}
      else muon_out_of_bounds = 1;
    }
  //std::cout << "MUON COUNTER IS: " << counter << std::endl; 
  counter = -1;

  for (const pat::Electron &el_index : *electrons) //tight electrons
    {
      //if (bad_event) break;
      counter++;
      if (el_index.et() > tight_electron_cutoff_pt && fabs(el_index.eta()) < electron_rapidity)
	{
	  if (lepton_found) {bad_event = true; extra_electron = 1; extra_electron_count++;} //std::cout << "BAD EVENT: more than one electron, or electron and muon.\n";}
	  else
	    {
	      lepton_found = true;
	      electron.SetPxPyPzE(el_index.px(), el_index.py(), el_index.pz(), sqrt(el_index.px()*el_index.px() + el_index.py()*el_index.py() + el_index.pz()*el_index.pz() + el_mass*el_mass));
	      electron_value_set = 1;
	      el_event = 1;
	      mu_event = 0;
	      el_mark = counter;
	      el_charge = el_index.charge();
	    }
	}
	else electron_out_of_bounds = 1;
    }
  //std::cout << "ELECTRON COUNTER IS: " << counter << std::endl;
  if (!lepton_found) {bad_event = true; no_lepton_found_count++;} //std::cout << "BAD EVENT: no leptons found.\n"; no_leptons = 1;}

  if (el_mark != -1 && el_charge < 0)
    {
      if (mu_mark != -1) {std::cout << "both el mark and mu mark can't be right.\n"; exit(3);}
      reco_lep_generation = 1;
      loose_lepton_cutoff_pt = loose_electron_cutoff_pt;
      lepton = electron;
      missing_et_cutoff = el_channel_missing_et_cutoff;
      if (el_event != 1 || mu_event != 0) {std::cout << "lepton misclassification.\n"; exit(3);}
    }
  if (el_mark != -1 && el_charge > 0)
    {
      if (mu_mark != -1) {std::cout << "both el mark and mu mark can't be right.\n"; exit(3);}
      reco_lep_generation = -1;
      loose_lepton_cutoff_pt = loose_electron_cutoff_pt;
      lepton = electron;
      missing_et_cutoff = el_channel_missing_et_cutoff;
      if (el_event != 1 || mu_event != 0) {std::cout << "lepton misclassification.\n"; exit(3);}
    }

  if (mu_mark != -1 && mu_charge < 0)
    {
      if (el_mark != -1) {std::cout << "both el mark and mu mark can't be right.\n"; exit(3);}
      reco_lep_generation = 2;
      loose_lepton_cutoff_pt = loose_muon_cutoff_pt;
      lepton = muon;
      missing_et_cutoff = mu_channel_missing_et_cutoff;
      if (mu_event != 1 || el_event != 0) {std::cout << "lepton misclassification.\n"; exit(3);}
    }
  if (mu_mark != -1 && mu_charge > 0)
    {
      //std::cout << "we're in a conditional where mu_mark isn't -1, it's " << mu_mark << std::endl;
      if (el_mark != -1) {std::cout << "both el mark and mu mark can't be right.\n"; exit(3);}
      reco_lep_generation = -2;
      loose_lepton_cutoff_pt = loose_muon_cutoff_pt;
      lepton = muon;
      missing_et_cutoff = mu_channel_missing_et_cutoff;
      if (mu_event != 1 || el_event != 0) {std::cout << "lepton misclassification.\n"; exit(3);}
    }

  if (mu_event == 1 && el_event == 0)
    {
      if (lepton.Px() == electron.Px() && lepton.Py() == electron.Py() && lepton.Pz() == electron.Pz() && lepton.E() == electron.E())
	{
	  std::cout << "What should be a muon looks like an electron.\n";
	  exit(4);
	}
      if (lepton.Px() != muon.Px() || lepton.Py() != muon.Py() || lepton.Pz() != muon.Pz() || lepton.E() != muon.E())
	{
	  std::cout << "It should be a muon, but it's not.\n";
	  exit(4);
	}
    }

 if (el_event == 1 && mu_event == 0)
    {
      if (lepton.Px() == muon.Px() && lepton.Py() == muon.Py() && lepton.Pz() == muon.Pz() && lepton.E() == muon.E())
	{
	  std::cout << "What should be an electron looks like a muon.\n";
	  exit(4);
	}
      if (lepton.Px() != electron.Px() || lepton.Py() != electron.Py() || lepton.Pz() != electron.Pz() || lepton.E() != electron.E())
	{
	  std::cout << "It should be an electron, but it's not.\n";
	  exit(4);
	}
    }

 if (muon_value_set == 1 && electron_value_set == 1) {std::cout << "Both lepton types set. Abort.\n"; exit(3);}
    
  counter = -1;

  for (const pat::Muon &mu_index : *muons) //loose muons
    {
      //if (bad_event) break;
      counter++;
      if (counter == mu_mark) continue; //avoid counting the tight muon as a loose muon
      if (mu_index.pt() > loose_lepton_cutoff_pt) {bad_event = true; loose_muon = 1; loose_muon_count++;}// std::cout << "BAD EVENT: loose muon.\n";}	  
    }
  //std::cout << "LOOSE MUON COUNTER IS: " << counter << std::endl;
  counter = -1;

  for (const pat::Electron &el_index : *electrons) //loose electrons
    {
      //if (bad_event) break;
      counter++;
      if (counter == el_mark) continue; //avoid counting the tight electron as a loose electron
      if (el_index.pt() > loose_lepton_cutoff_pt) {bad_event = true; loose_electron = 1; loose_electron_count++;} // std::cout << "BAD EVENT: loose electron.\n";}
    }
  //std::cout << "LOOSE ELECTRON COUNTER IS: " << counter << std::endl;

  //MARKER NSUBJETTINESS

  //   _   _      ____        _     _      _   _   _                     
  //  | \ | |    / ___| _   _| |__ (_) ___| |_| |_(_)_ __   ___  ___ ___ 
  //  |  \| |____\___ \| | | | '_ \| |/ _ \ __| __| | '_ \ / _ \/ __/ __|
  //  | |\  |_____|__) | |_| | |_) | |  __/ |_| |_| | | | |  __/\__ \__ \ |
  //  |_| \_|    |____/ \__,_|_.__// |\___|\__|\__|_|_| |_|\___||___/___/
  //                             |__/                                    
  /*
  edm::Handle<pat::JetCollection> AK8CHS;
  iEvent.getByToken(ak8jetToken_, AK8CHS);

  for (const pat::Jet &ijet: *AK8CHS)
    {
      double tau1 = ijet.userFloat("NjettinessAK8:tau1");
      double tau2 = ijet.userFloat("NjettinessAK8:tau2");
      double tau3 = ijet.userFloat("NjettinessAK8:tau3");

      double tau21 = 99.0;
      double tau32 = 99.0;

      if (tau1 != 0.0) tau21 = tau2/tau1;
      if (tau2 != 0.0) tau32 = tau3/tau2;

      std::cout << "tau21 is " << tau21 << std::endl;
      std::cout << "tau32 is " << tau32 << std::endl;

      double puppi_tau1 = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
      double puppi_tau2 = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
      double puppi_tau3 = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");

      double puppi_tau21 = 99.0;
      double puppi_tau32 = 99.0;

      if (puppi_tau1 != 0.0) puppi_tau21 = puppi_tau2/puppi_tau1;
      if (puppi_tau2 != 0.0) puppi_tau32 = puppi_tau3/puppi_tau2;

      std::cout << "puppi_tau21 is " << puppi_tau21 << std::endl;
      std::cout << "puppi_tau32 is " << puppi_tau32 << std::endl;
    }
  */
#if 0
  edm::EDGetTokenT<pat::JetCollection>     fTokPatJetName;
  fTokPatJetName    = iC.consumes<pat::JetCollection>   (fJetName);
  // Get jet collection
  edm::Handle<pat::JetCollection> hJetProduct;
  iEvent.getByToken(fTokPatJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const pat::JetCollection *jetCol = hJetProduct.product();//NEEDED FOR N-SUBJETTINESS TO ITERATE OVER --> jet_iterator_AOD_fill ->jetCol-begin to end;
//jet_iterator_AOD_fill input to addjet from fill, to get vartype_pfjet_input_to_addjet
 for ( ) 
  edm::EDGetTokenT<edm::ValueMap<float> >    fTokQGLSubJets      ;
  edm::EDGetTokenT<reco::PFJetCollection>    fTokSubJets         ;
  edm::EDGetTokenT<reco::JetTagCollection>   fTokCSVbtagSubJetName;

  //Get Quark Gluon Likelihood on subjets
  edm::Handle<edm::ValueMap<float> > hQGLikelihoodSubJets; //NEEDED FOR N-SUBJETTINESS
  iEvent.getByToken(fTokQGLSubJets,hQGLikelihoodSubJets);
  assert(hQGLikelihoodSubJets.isValid());

    // Get sub-jet collection  
  edm::Handle<reco::PFJetCollection> hSubJetProduct;
  iEvent.getByToken(fTokSubJets,hSubJetProduct);
  assert(hSubJetProduct.isValid());
  const reco::PFJetCollection *subJetCol = hSubJetProduct.product(); //NEEDED FOR N-SUBJETTINESS

  // Get b sub-jets
  edm::Handle<reco::JetTagCollection> hCSVbtagsSubJets; //NEEDED FOR N-SUBJETTINESS
  iEvent.getByToken(fTokCSVbtagSubJetName, hCSVbtagsSubJets);
  assert(hCSVbtagsSubJets.isValid());
  
  const double fConeSize = 0.4;
  const reco::PFJet &vartype_pfjet_input_to_addjet;
  //==================================================================================START
  // find/sort up to 4 hardest subjets
  const reco::PFJet *subjet1=0, *subjet2=0, *subjet3=0, *subjet4=0;
  double csv1=-2, csv2=-2, csv3=-2, csv4=-2;
  double qgid1=-2, qgid2=-2, qgid3=-2, qgid4=-2;
  double q1=-100, q2=-100, q3=-100, q4=-100;
  for(reco::PFJetCollection::const_iterator subjet_iterator = subJetCol->begin(); subjet_iterator!=subJetCol->end(); ++subjet_iterator)
    {
      if(reco::deltaR(vartype_pfjet_input_to_addjet.eta(),vartype_pfjet_input_to_addjet.phi(),subjet_iterator->eta(),subjet_iterator->phi())>fConeSize) continue;  // (!) get associated subjets by dR...is there a better way???

      reco::PFJetRef subjetRef(hSubJetProduct, subjet_iterator - subJetCol->begin());
      reco::JetBaseRef subjetBaseRef(subjetRef);

      if(!subjet1 || subjet_iterator->pt() > subjet1->pt())
	{
	  subjet4 = subjet3;
	  csv4    = csv3;
	  qgid4   = qgid3;
	  q4      = q3;
      
	  subjet3 = subjet2;
	  csv3    = csv2;
	  qgid3   = qgid2;
	  q3      = q2;
      
	  subjet2 = subjet1;
	  csv2    = csv1;
	  qgid2   = qgid1;
	  q2      = q1;
      
	  subjet1 = &(*subjet_iterator);      
	  csv1    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];      
	  qgid1   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
	  q1      = jetCharge(*subjet_iterator);
      
	}
      else if(!subjet2 || subjet_iterator->pt() > subjet2->pt())
	{
	  subjet4 = subjet3;
	  csv4    = csv3;
	  qgid4   = qgid3;
	  q4      = q3;

	  subjet3 = subjet2;
	  csv3    = csv2;
	  qgid3   = qgid2;
	  q3      = q2;

	  subjet2 = &(*subjet_iterator);      
	  csv2    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef]; 
	  qgid2   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
	  q2      = jetCharge(*subjet_iterator);

	}
      else if(!subjet3 || subjet_iterator->pt() > subjet3->pt())
	{
	  subjet4 = subjet3;
	  csv4    = csv3;
	  qgid4   = qgid3;
	  q4      = q3;

	  subjet3 = &(*subjet_iterator);
	  csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
	  qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
	  q3      = jetCharge(*subjet_iterator);
      
	}
      else if(!subjet4 || subjet_iterator->pt() > subjet4->pt())
	{
	  subjet4 = &(*subjet_iterator);
	  csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
	  qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
	  q4      = jetCharge(*subjet_iterator);
	}
    }
  //=============================================================FINISH
#endif

  // MARKER RESOLVED JETS
  //   ____  _____ ____   ___  _ __     _______ ____        _ _____ _____ ____  
  //  |  _ \| ____/ ___| / _ \| |\ \   / / ____|  _ \      | | ____|_   _/ ___| 
  //  | |_) |  _| \___ \| | | | | \ \ / /|  _| | | | |  _  | |  _|   | | \___ \ |
  //  |  _ <| |___ ___) | |_| | |__\ V / | |___| |_| | | |_| | |___  | |  ___) |
  //  |_| \_\_____|____/ \___/|_____\_/  |_____|____/   \___/|_____| |_| |____/ 

  if (bad_event) return; //{std::cout << "bad event reached resolved jets\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";

  //start jet algorithm here
  ///////////////////////////////////GOOD JETS

  edm::Handle<pat::JetCollection> ak8_jets;// max = 4; mode = 2; mean = 1.7;
  iEvent.getByToken(ak8jetToken_, ak8_jets);

  edm::Handle<reco::GenJetCollection> gen_ak8_jets;// max = 5; mode = 2; mean = 2.0;
  iEvent.getByToken(ak8genjetToken_, gen_ak8_jets);

  TLorentzVector jet0, jet1;

  jet0.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
  jet1.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

  Double_t tau21(99.9), tau32(99.9), puppi_tau21(99.9), puppi_tau32(99.9);
  /////////////////////////////////////////////////////////////////////////////////////////////
  Int_t jet0_mark = -1;
  Int_t jet1_mark = -1;

  no_jets = 0;
  lonesome_jet = 0;
  extra_merged_jet = 0;
  bad_dijet = 0;
  loose_jet = 0;

  const Double_t jet8_rapidity = 2.4; //was infinity
  //const Double_t jet4_rapidity = 2.4;
  const Double_t goodjet_pt_cutoff = 200.0; //was 40.0 (35.0 second)
  //const Double_t badjet_pt_cutoff = 30.0; //was 30.0
  const Double_t Tau21max = 0.55;
  
  counter = -1;
  Double_t leading_jet_pt = goodjet_pt_cutoff; //set initial pt to cutoff value, increase as more jets come in

  for (const pat::Jet &jet_index : *ak8_jets) //search for leading jet
    {
      //if (bad_event) break;
      counter++;
      if (jet_index.pt() > leading_jet_pt && fabs(jet_index.eta()) < jet8_rapidity)
	{
	  leading_jet_pt = jet_index.pt();
	  jet0.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	  jet0_mark = counter;
	  backup_skip_counter += counter;
          Double_t tau1 = jet_index.userFloat("NjettinessAK8:tau1");
	  Double_t tau2 = jet_index.userFloat("NjettinessAK8:tau2");
	  Double_t tau3 = jet_index.userFloat("NjettinessAK8:tau3");

	  if (tau1 != 0.0) tau21 = tau2/tau1;
	  if (tau2 != 0.0) tau32 = tau3/tau2;

	  Double_t puppi_tau1 = jet_index.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
	  Double_t puppi_tau2 = jet_index.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
	  Double_t puppi_tau3 = jet_index.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");

	  if (puppi_tau1 != 0.0) puppi_tau21 = puppi_tau2/puppi_tau1;
	  if (puppi_tau2 != 0.0) puppi_tau32 = puppi_tau3/puppi_tau2;
	}
    }
  //std::cout << "JET COUNTER IS: " << counter << std::endl;
  if (tau21 >= Tau21max) {bad_event = true;}// std::cout << "njettiness rejection.\n";}
   
  counter = -1;
  Double_t sub_leading_jet_pt = goodjet_pt_cutoff; //again, set to cutoff value, skipping jet already marked
    
  if (jet0_mark == -1) {bad_event = true; no_jet_count++;} //std::cout << "BAD EVENT: UNMERGED: no jets found.\n"; n o jets = 1;}

  for (const pat::Jet &jet_index : *ak8_jets) //search for subleading jet
    {
      //if (bad_event) break;
      counter++;
      if (counter == jet0_mark) {skip_counter += counter; continue;} //avoid treating leading jet as subleading jet
      if (jet_index.pt() > sub_leading_jet_pt && fabs(jet_index.eta()) < jet8_rapidity)
	{	   
	  sub_leading_jet_pt = jet_index.pt();
	  jet1.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	  jet1_mark = counter;
	  merged_counter += counter;
	  //outfile << "jet1_mark is " << jet1_mark << std::endl;
	}
    }
  //std::cout << "SECOND JET COUNTER IS: " << counter << std::endl;
  //////////////////////////////BAD JETS

  counter = -1;
  //Double_t background_jet_pt_max = badjet_pt_cutoff;
  //if (talking) outfile << "merged event is " << merged_event << std::endl;
  //const int zeroth_tester = merged_event;

  if (jet1_mark == -1) {lonesome_jet = 1; merged_event = 1;}// outfile << "jet1_mark is " << jet1_mark << std::endl;}//bad_event = true;// std::cout << "BAD EVENT: UNMERGED: only one jet.\n"; lone some jet = 1;}
  else merged_event = 0;

  //if (talking) outfile << "merged event is " << merged_event << std::endl;
  //const int first_tester = merged_event;
  
  //std::cout << "BAD JET COUNTER IS: " << counter << std::endl;
  //std::cout << "COUNTER\n";
  TLorentzVector dijet_system;
  dijet_system.SetPxPyPzE(jet0.Px() + jet1.Px(), jet0.Py() + jet1.Py(), jet0.Pz() + jet1.Pz(), jet0.E() + jet1.E());
  
  if (dijet_system.Pt() <= 70.0) bad_dijet = 1;

  // MARKER MERGED JETS
  //
  //  m    m mmmmmm mmmmm    mmm  mmmmmm mmmm          mmm  mmmmmm mmmmmmm  mmmm 
  //  ##  ## #      #   "# m"   " #      #   "m          #  #         #    #"   "
  //  # ## # #mmmmm #mmmm" #   mm #mmmmm #    #          #  #mmmmm    #    "#mmm 
  //  # "" # #      #   "m #    # #      #    #          #  #         #        "#
  //  #    # #mmmmm #    "  "mmm" #mmmmm #mmm"       "mmm"  #mmmmm    #    "mmm#"
  //
	
  ///////////////////////////////////////////////////////////////MERGED STUFF

  if (bad_event) return; //{std::cout << "bad event reached merged jets\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";

  edm::Handle<pat::JetCollection> sub_jets;// max = 7; mode = 3; mean = 2.5;
  iEvent.getByToken(ak8CHSSoftDropSubjetsToken_, sub_jets);

  TLorentzVector subjet0, subjet1;
  Int_t good_jet = 0;

  subjet0.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);
  subjet1.SetPxPyPzE(-99.9, -99.9, -99.9, -99.9);

  Int_t count0 = -1;
  Int_t count1 = -1;

  min_submass_diff = 100.0;

  for (const pat::Jet subjet_index0 : *sub_jets)
    {
      count0++;
      TLorentzVector dummy0;
      dummy0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
      if (dummy0.DeltaR(jet0) > 1.0) continue; //FIXME: NOT SURE IF THIS IS CORRECT; PROBABLY NOT
      for (const pat::Jet subjet_index1 : *sub_jets)
	{
	  count1++;
	  TLorentzVector dummy1;
	  dummy1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
	  if (dummy1.DeltaR(jet0) > 1.0) continue; //FIXME: NOT SURE IF THIS IS CORRECT; PROBABLY NOT
	  if (count1 == count0) continue;
	  TLorentzVector dummy_sum = dummy0 + dummy1;
	  Double_t two_body_mass = std::sqrt(dummy_sum.Dot(dummy_sum));
	  Double_t mass_diff = fabs(two_body_mass - Wmass);
	  if (mass_diff < min_submass_diff)
	    {
	      min_submass_diff = mass_diff;
	      subjet0 = dummy0;
	      subjet1 = dummy1;
	      good_jet = 1;
	    }
	}
      count1 = -1;
    }
  if (min_submass_diff > 20.0) bad_event = true;
  
  /*
  Int_t subjet0_mark = -1;
  //Int_t subjet1_mark = -1;

  Int_t merged_mark = -1;
  Int_t subcounter = -1;
 
  //if boosted (merged?), reconstructed w_pt must be over 200.0
  Double_t merged_jet_pt_min = 200.0;
  counter = -1;
  const Double_t merged_otherjet_pt_max = 150;//80.0;//seems to be a killer
	
  for (const pat::Jet &jet_index : *ak8_jets)
    {
      //if (bad_event) break;
      counter++;
      if (jet_index.pt() > merged_jet_pt_min)
	{
	  merged_jet_pt_min = jet_index.pt();
	  //merged_jet.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
	  merged_mark = counter;
	}
    }

  if (merged_mark == -1 && !bad_event) {no_merged_jet = 1; std::cout << "no merged jet.\n"; exit(2);}
  if (merged_mark == -1 && merged_event == 1) {bad_event = true; merged_combo_count++;}
  if (talking) outfile << "merged event is " << merged_event << std::endl;
  const int second_tester = merged_event;
  counter = -1;

  for (const pat::Jet &jet_index : *ak8_jets)
    {
      //if (bad_event) break;
      counter++;
      if (counter == merged_mark) continue; //skip over already found jet
      if (jet_index.pt() > merged_otherjet_pt_max) extra_merged_jet = 1;
    }

  if (extra_merged_jet == 1 && merged_event == 1) {bad_event = true; extra_merged_combo_count++;}//doesn't add up
  const int third_tester = merged_event;
  if (extra_merged_jet != 1 && merged_event == 1) nonmerged_count++;
  const int fourth_tester = merged_event;
  if (extra_merged_jet == 1 && merged_event != 1) inverted_count++;
  const int fifth_tester = merged_event;
  if (extra_merged_jet != 1 && merged_event != 1) neither_count++;
  const int sixth_tester = merged_event;
  if (merged_event) merged_event_count++;
  const int seventh_tester = merged_event;
  counter = -1;
	
  //if (zeroth_tester != first_tester) {std::cout << "zeroth and first don't match. zeroth is " << zeroth_tester << " first is " << first_tester << std::endl;}
  if (first_tester != second_tester) {std::cout << "first and second don't match. first is " << first_tester << " second is " << second_tester << std::endl; exit(43);}
  if (second_tester != third_tester) {std::cout << "second and third don't match. second is " << second_tester << " third is " << third_tester << std::endl; exit(43);}
  if (third_tester != fourth_tester) {std::cout << "third and fourth don't match. third is " << third_tester << " fourth is " << fourth_tester << std::endl; exit(43);}
  if (fourth_tester != fifth_tester) {std::cout << "fourth and fifth don't match. fourth is " << fourth_tester << " fifth is " << fifth_tester << std::endl; exit(43);}
  if (fifth_tester != sixth_tester) {std::cout << "fifth and sixth don't match. fifth is " << fifth_tester << " sixth is " << sixth_tester << std::endl; exit(43);}
  if (sixth_tester != seventh_tester) {std::cout << "sixth and seventh don't match. sixth is " << sixth_tester << " seventh is " << seventh_tester << std::endl; exit(43);}
  //now for the subjets
  Double_t subjet_pt_max = 0.0;
  Double_t subjet_pt_next_to_max = 0.0;
  for (const pat::Jet &jet_index : *ak8_jets) //HERE'S THE TWO SUB-LOOPS FOR THE SUBJETS 0 AND 1
    {
      if (jet_index.pt() == 0) { } //keep compiler from warning about unused index
      //if (bad_event) break;
      counter++;
      if (counter != merged_mark) continue; //this time we want only the merged jet we found
	    
      for (const pat::Jet &subjet_index : *sub_jets)
	{
	  subcounter++;
	  if (subjet_index.pt() > subjet_pt_max)
	    {
	      subjet_pt_max = subjet_index.pt();
	      subjet0.SetPxPyPzE(subjet_index.px(), subjet_index.py(), subjet_index.pz(), subjet_index.energy());
	      subjet0_mark = subcounter;
	    }
	}

      subcounter = -1;
	    
      for (const pat::Jet &subjet_index : *sub_jets)
	{
	  subcounter++;
	  if (subcounter == subjet0_mark) continue; //don't want to count the same subjet twice
	  if (subjet_index.pt() > subjet_pt_next_to_max)
	    {
	      subjet_pt_next_to_max = subjet_index.pt();
	      subjet1.SetPxPyPzE(subjet_index.px(), subjet_index.py(), subjet_index.pz(), subjet_index.energy());
	      //subjet1_mark = subcounter;
	    }
	}
    }
  //*/ 
  ///////////////////////////////////////////////////////////////END MERGED STUFF

 
  
  //MARKER MET
  //
  // #    #  ######   #####
  // ##  ##  #          #
  // # ## #  #####      #
  // #    #  #          #
  // #    #  #          #
  // #    #  ######     #
  //

  if (bad_event) return; //{std::cout << "bad event reached MET.\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
  
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();

  //MISSING ET

  bool angle_test = false;
	
  TLorentzVector met4vector, other_met4vector;
  Double_t met_px = met.px();
  Double_t met_py = met.py();
  //if (mu_event == 1) met_py -= 6.2;
  //if (el_event == 1) met_py -= 7.7;
  met4vector.SetPxPyPzE(met_px, met_py, 0.0, 0.0);
  if (angle_test) met4vector.SetPxPyPzE(neutrino_gen.Px(), neutrino_gen.Py(), 0.0, neutrino_gen.Et());//ONLY A TEST
 
  met_small = 0;
  undefined_met = 0;
	
  METzCalculator calc;
  calc.SetMET(met4vector);
  calc.SetLepton(lepton);
  calc.SetTruthInfo(neutrino_gen);
  if (angle_test) calc.SetLepton(lepton_gen);//ONLY A TEST

  if (abs(reco_lep_generation) == 1) calc.SetLeptonType("electron");// std::cout << "DEBUGGING: lepton is an electron.\n";}
  if (abs(reco_lep_generation) == 2) calc.SetLeptonType("muon");// std::cout << "DEBUGGING: lepton is a muon.\n";}

  Double_t the_mass = lepton.M(); //ONLY FOR A TEST; COMMENTED OUT FOR A FURTHER TEST

  if (abs(reco_lep_generation) == 1) the_mass = el_mass;
  if (abs(reco_lep_generation) == 2) the_mass = mu_mass;
 
  calc.SetWmass(Wmass);
  TLorentzVector wlep_gen = lepton_gen + neutrino_gen;//JUST TO TEST
  //Double_t wgenmass = wlep_gen.M();//JUST TO TEST
  if (angle_test) calc.SetWmass(wlep_gen.M());

  if (!angle_test) calc.SetJets(subjet0, subjet1);
  if (angle_test)  calc.SetJets(quark_gen, antiquark_gen);

  calc.SetnMass(0.0);
  Double_t neu_gen_mass = neutrino_gen.M();
  if (angle_test && (neu_gen_mass < 0.0 || the_mass < 0.0)) bad_event = true;
  if (angle_test) calc.SetnMass(neu_gen_mass);
  
  //if (gen_lep_generation < 0) wgenmass = Wplus_gen.M();//JUST TO TEST
  // if (gen_lep_generation > 0) wgenmass = Wminus_gen.M();//JUST TO TEST

  Double_t met_pz = calc.Calculate(7);
  Double_t other_met_pz = calc.getOther();

  Double_t fill_pz0 = calc.Calculate(0);
  Double_t fill_pz1 = calc.Calculate(1);
  Double_t fill_pz2 = calc.Calculate(2);
  Double_t fill_pz3 = calc.Calculate(3);
  Double_t fill_pz4 = calc.Calculate(4);
  Double_t fill_pz5 = calc.Calculate(5);
  Double_t fill_pz6 = calc.Calculate(6);
  Double_t fill_pz7 = calc.Calculate(7);
 
      
  if (calc.IsComplex()) imaginary_neutrino = 1;
  else imaginary_neutrino = 0;
  
  if (met_pz != met_pz) {bad_event = true;}// std::cout << "BAD EVENT: met pz not well defined.\n"; undefined_met = 1;}
  //else {}
 
  reco_konstant = Wmass * Wmass - the_mass * the_mass + 2.0 * (lepton.Px() * met4vector.Px() + lepton.Py() * met4vector.Py());
  gen_konstant = wlep_gen.M() * wlep_gen.M() - lepton_gen.M() * lepton_gen.M() + 2.0*(lepton_gen.Px()*met4vector.Px() + lepton_gen.Py()*met4vector.Py()) - neu_gen_mass*neu_gen_mass;//ONLY FOR A TEST

  Double_t met_e = (met_pz * lepton.Pz() + 0.5 * reco_konstant) / lepton.E();
  Double_t other_met_e = (other_met_pz * lepton.Pz() + 0.5 * reco_konstant) / lepton.E();

  if (angle_test) met_e = (met_pz * lepton_gen.Pz() + 0.5 * gen_konstant) / lepton_gen.E();//ONLY TO TEST
  if (angle_test) other_met_e = (other_met_pz * lepton_gen.Pz() + 0.5 * gen_konstant) / lepton_gen.E();//ONLY TO TEST
  
  if (!bad_event) met4vector.SetPxPyPzE(met_px, met_py, met_pz, met_e);
  if (!bad_event) other_met4vector.SetPxPyPzE(met_px, met_py, other_met_pz, other_met_e);

  if (angle_test && !bad_event) met4vector.SetPxPyPzE(neutrino_gen.Px(), neutrino_gen.Py(), met_pz, met_e); //JUST TO TEST
  if (angle_test && !bad_event) other_met4vector.SetPxPyPzE(neutrino_gen.Px(), neutrino_gen.Py(), other_met_pz, other_met_e); //JUST TO TEST
  
  if (met4vector.Et() < missing_et_cutoff) {bad_event = true; missing_et_cutoff_count++;}
  //if (met4vector.M() < -0.0) {bad_event = true;}

  Double_t calculated_mass = met4vector.M();

  if (!bad_event && talking && calculated_mass != 0.0) 
    {
      outfile << "input is nu x: " << met4vector.Px() << " nu y: " << met4vector.Py() << " nu z: " << met4vector.Pz() << " nu e: " << met4vector.E() << std::endl;
      outfile << "input is lep x: " << lepton.Px() << " lep y: " << lepton.Py() << " lep z: " << lepton.Pz() << " lep e: " << lepton.E() << std::endl;
      outfile << "lepton mass = " << the_mass << std::endl;
      outfile << "calculated pz is " << met_pz << " and calculated energy is " << met_e << std::endl;
      outfile << "or perhaps calculated pz is " << other_met_pz << " and calculated energy is " << other_met_e << std::endl;
      outfile << "calculated mass is " << calculated_mass << std::endl;
      
      //std::cout << "BAD EVENT: reco lep generation is: " << reco_lep_generation << std::endl;
      //std::cout << "BAD EVENT: electron mark is: " << el_mark << std::endl;
      //std::cout << "BAD EVENT: muon mark is: " << mu_mark << std::endl;
      //std::cout << "BAD EVENT: THAT'S ALL FOLKS!\n\n";
    }

  /* if (no_leptons == 1 || extra_muon == 1 || extra_electron == 1)
    {
      loose_electron = 0;
      loose_muon = 0;
      met_small = 0;
      }*/
  if (no_jets == 1) bad_event = true;
  //  {
  //  lone some_jet = 0;
  //  loose_jet = 0;
  //  met_small = 0;
  //  }
  if (lonesome_jet == 1)
    {
      loose_jet = 0;
      met_small = 0;
    }


  ////////SANITY CHECK
  // if (no_leptons == 1)
  //  {
  //    //counter = -1;
  //    Double_t pt_mu_max = 0.0;
  //    Double_t pt_el_max = 0.0;
  //    Double_t eta_of_mu_max = 0.0;
  //    Double_t eta_of_el_max = 0.0;
//	
  //    for (const pat::Muon &mu_index : *muons)
//	{
//	  if (mu_index.pt() > pt_mu_max)// && fabs(mu_index.eta()) < muon_rapidity)
//	    {
//	      pt_mu_max = mu_index.pt();
//	      eta_of_mu_max = mu_index.eta(); eta_of_mu_max = eta_of_mu_max;
//	    }
//	}
//
  //    for (const pat::Electron &el_index : *electrons)
//	{
//	  if (el_index.et() > pt_el_max)// && fabs(el_index.eta()) < electron_rapidity)
//	    {
//	      pt_el_max = el_index.pt();
//	      eta_of_el_max = el_index.eta(); eta_of_el_max = eta_of_el_max;
//	    }
//	}
//
//	}


  //MARKER GAUGE BOSONS

  //   ____                          ____                            
  //  / ___| __ _ _   _  __ _  ___  | __ )  ___  ___  ___  _ __  ___ 
  // | |  _ / _` | | | |/ _` |/ _ \ |  _ \ / _ \/ __|/ _ \| '_ \/ __|
  // | |_| | (_| | |_| | (_| |  __/ | |_) | (_) \__ \ (_) | | | \__ \ |
  //  \____|\__,_|\__,_|\__, |\___| |____/ \___/|___/\___/|_| |_|___/
  //                    |___/                                        

  if (bad_event) return; //{std::cout << "bad event reached gauge bosons.\n"; return;}
   if (bad_event) std::cout << "you shouldn't be seeing this.";

  TLorentzVector Whad;
  Whad.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);
  TLorentzVector Wlep;
  Wlep.SetPxPyPzE(-99.9,-99.9,-99.9,-99.9);

  Int_t Wlep_set = 0;
  
  if (mu_event == 1 && el_event == 0)
    {
      //if (muon_value_set == 0) {std::cout << "you screwed up again.\n"; exit(23);}
      Wlep_set = 2;
      Wlep = muon + met4vector;
    }
  //else if (muon_value_set == 1) {std::cout << "you screwed up again.\n"; exit(24);}
  
  if (mu_event == 0 && el_event == 1)
    {
      //if (electron_value_set == 0) {std::cout << "you screwed up again.\n"; exit(25);}
      Wlep_set = 1;
      Wlep = electron + met4vector;
    }
  //else if (electron_value_set == 1) {std::cout << "you screwed up again.\n"; exit(26);}
  
  //if (merged_event == 0) Whad = jet0 + jet1;
  //if (merged_event == 1)
  Whad = subjet0 + subjet1;

  
  //MARKER: TRUTH MATCHING STUFF IS GETTING IN THE WAY
  // _____           _   _       __  __       _       _     _             
  //|_   _| __ _   _| |_| |__   |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
  //  | || '__| | | | __| '_ \  | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
  //  | || |  | |_| | |_| | | | | |  | | (_| | || (__| | | | | | | | (_| |
  //  |_||_|   \__,_|\__|_| |_| |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
  //                                                                |___/ 
#if 0
  //////////////////////////////////////TIME TO TEST WMASS FIT
  //
  min_mass_diff = 1000.0;
  Int_t counter0 = -1;
  Int_t counter1 = -1;
  for (const pat::Jet &jet_index0 : *ak8_jets)
    {
      counter0++;
      for (const pat::Jet &jet_index1 : *ak8_jets)
	{
	  counter1++;
	  if (counter0 == counter1) continue;
	  //TLorentzVector dummy0, dummy1;
	  TLorentzVector dummy;
	  dummy.SetPxPyPzE(jet_index0.px() + jet_index1.px(), jet_index0.py() + jet_index1.py(), jet_index0.pz() + jet_index1.pz(), jet_index0.energy() + jet_index1.energy());
	  //dummy0.SetPxPyPzE(jet_index0.px(), jet_index0.py(), jet_index0.pz(), jet_index0.energy());
	  //dummy1.SetPxPyPzE(jet_index1.px(), jet_index1.py(), jet_index1.pz(), jet_index1.energy());
	  //TLorentzVector sum = dummy0 + dummy1;
	  //Double_t two_body_mass = std::sqrt(sum.Dot(sum));
	  Double_t two_body_mass = std::sqrt(dummy.Dot(dummy));
	  Double_t mass_diff = fabs(two_body_mass - Wmass);
	  if (mass_diff < min_mass_diff) min_mass_diff = mass_diff;
	}
      counter1 = -1;
    }
  //if (min_mass_diff > 200.0) bad_event = true;

  ///////////////////////////////////////////////////////actual truth matching

  min_quark_delta = 1.0;

  for (const pat::Jet &jet_index : *ak8_jets)
    {
      TLorentzVector dummy;
      dummy.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      Double_t delta_r = dummy.DeltaR(quark_gen);
      if (delta_r < min_quark_delta) min_quark_delta = delta_r;
    }

  min_antiquark_delta = 1.0;

  for (const pat::Jet &jet_index : *ak8_jets)
    {
      TLorentzVector dummy;
      dummy.SetPxPyPzE(jet_index.px(), jet_index.py(), jet_index.pz(), jet_index.energy());
      Double_t delta_r = dummy.DeltaR(antiquark_gen);
      if (delta_r < min_antiquark_delta) min_antiquark_delta = delta_r;
    }
  //
  ///////////////////////////////////////////////////
  
  //////////////////////////////////////TIME TO TEST WMASS FIT for subjets
  //
  min_submass_diff = 1000.0;
  Int_t subcounter0 = -1;
  Int_t subcounter1 = -1;
  for (const pat::Jet &subjet_index0 : *sub_jets)
    {
      subcounter0++;
      for (const pat::Jet &subjet_index1 : *sub_jets)
	{
	  subcounter1++;
	  if (subcounter0 == subcounter1) continue;
	  //TLorentzVector dummy0, dummy1;
	  TLorentzVector dummy;
	  dummy.SetPxPyPzE(subjet_index0.px() + subjet_index1.px(), subjet_index0.py() + subjet_index1.py(), subjet_index0.pz() + subjet_index1.pz(), subjet_index0.energy() + subjet_index1.energy());
	  //dummy0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
	  //dummy1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
	  //TLorentzVector sum = dummy0 + dummy1;
	  //Double_t two_body_mass = std::sqrt(sum.Dot(sum));
	  Double_t two_body_mass = std::sqrt(dummy.Dot(dummy));
	  Double_t mass_diff = fabs(two_body_mass - Wmass);
	  if (mass_diff < min_submass_diff) min_submass_diff = mass_diff;
	}
      subcounter1 = -1;
    }
  //
  /////////////////////////////////////////////////
  
  //re-do merged stuff, just to test.
  min_submass_diff = 1000.0;
  subcounter0 = -1;
  subcounter1 = -1;
  for (const pat::Jet &subjet_index0 : *sub_jets)
    {
      subcounter0++;
      for (const pat::Jet &subjet_index1 : *sub_jets)
	{
	  subcounter1++;
	  if (subcounter1 == subcounter0) continue;
	  TLorentzVector dummy0;
	  TLorentzVector dummy1;
	  dummy0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
	  dummy1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
	  TLorentzVector dummy_sum = dummy0 + dummy1;
	  Double_t two_body_mass = std::sqrt(dummy_sum.Dot(dummy_sum));
	  Double_t mass_diff = fabs(two_body_mass - Wmass);
	  if (mass_diff < min_submass_diff) min_submass_diff = mass_diff;
	  subjet0.SetPxPyPzE(subjet_index0.px(), subjet_index0.py(), subjet_index0.pz(), subjet_index0.energy());
	  subjet1.SetPxPyPzE(subjet_index1.px(), subjet_index1.py(), subjet_index1.pz(), subjet_index1.energy());
	}
      subcounter1 = -1;
      }
  /////////////////////////////////////////
#endif
  
  //MARKER RHO
  //
  //
  // ######  #     # #######
  // #     # #     # #     #
  // #     # #     # #     #
  // ######  ####### #     #
  // #   #   #     # #     #
  // #    #  #     # #     #
  // #     # #     # #######
  //
  //
#if 0  
  Handle<double> rho_Handle;
  iEvent.getByToken(rhoToken_, rho_Handle);
  double rho_var = *rho_Handle;

  //MARKER JEC PAYLOAD
  //                                                                         
  //   mmm  mmmmmm   mmm         mmmmm    mm  m     m m       mmmm    mm   mmmm  
  //     #  #      m"   "        #   "#   ##   "m m"  #      m"  "m   ##   #   "m
  //     #  #mmmmm #             #mmm#"  #  #   "#"   #      #    #  #  #  #    #
  //     #  #      #             #       #mm#    #    #      #    #  #mm#  #    #
  // "mmm"  #mmmmm  "mmm"        #      #    #   #    #mmmmm  #mm#  #    # #mmm" 
  //
  
  // Run 2016 F into two JEC payloads. 
  //   IOV EF: [276831,278801] corresponds to Summer16_23Sep2016EFV3_DATA (For Runs E/early F)
  //   IOV FG: [278802,280385] corresponds to Summer16_23Sep2016GV3_DATA (For Runs lateF/G)
  std::vector<std::string>  jecPayloadsAK4chsFinal;
  //std::vector<std::string>  jecPayloadsAK8chsFinal;
  //std::vector<std::string>  jecPayloadsAK4pupFinal;
  //std::vector<std::string>  jecPayloadsAK8pupFinal;

  int run_Number = iEvent.id().run();

  if (isRun2016F_ && run_Number < 278801)
    {
      std::cout << "Using Run2016F early JEC" << std::endl;
      jecPayloadsAK4chsFinal = jecPayloadsAK4chs_;
      //jecPayloadsAK8chsFinal = jecPayloadsAK8chs_;
      //jecPayloadsAK4pupFinal = jecPayloadsAK4pup_;
      //jecPayloadsAK8pupFinal = jecPayloadsAK8pup_;
    } 
  else if (isRun2016F_ && run_Number > 278801)
    {
      std::cout << "Using Run2016F late JEC" << std::endl;
      jecPayloadsAK4chsFinal = jecPayloadsAK4chsSecondary_;
      //jecPayloadsAK8chsFinal = jecPayloadsAK8chsSecondary_;
      //jecPayloadsAK4pupFinal = jecPayloadsAK4pupSecondary_;
      //jecPayloadsAK8pupFinal = jecPayloadsAK8pupSecondary_;
    } 
  else
    {
      std::cout << "Using Primary JEC from cfg" << std::endl;
      jecPayloadsAK4chsFinal = jecPayloadsAK4chs_;
      //jecPayloadsAK8chsFinal = jecPayloadsAK8chs_;
      //jecPayloadsAK4pupFinal = jecPayloadsAK4pup_;
      //jecPayloadsAK8pupFinal = jecPayloadsAK8pup_;
    }
  
  if (true)
    {
      std::cout << "jecPayloadsAK4chs_.size()              " << jecPayloadsAK4chs_.size() << std::endl;
      std::cout << "jecPayloadsAK4chsSecondary_.size()     " << jecPayloadsAK4chsSecondary_.size() << std::endl;
      std::cout << "jecPayloadsAK4chsFinal.size()          " << jecPayloadsAK4chsFinal.size() << std::endl;
    }//ADDED:..AND ENDING HERE, BRAND NEW

  // AK4chs JEC 
  std::vector<JetCorrectorParameters> vect_param_ak4;
  for (std::vector<std::string>::const_iterator iterator = jecPayloadsAK4chsFinal.begin(),
	 iterator_end = jecPayloadsAK4chsFinal.end(); iterator != iterator_end - 1; ++iterator)
    {
      JetCorrectorParameters params(*iterator);
      vect_param_ak4.push_back(params);
    }
  
  JetCorrectorAK4chs   = boost::shared_ptr<FactorizedJetCorrector>  (new FactorizedJetCorrector(vect_param_ak4));
  JetCorrUncertAK4chs  = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecPayloadsAK4chsFinal.back()));

 // jet resolution from text files //ADDED:STARTING HERE...
  JME::JetResolution jet_resolution_AK4CHS;
  //JME::JetResolution jet_resolution_AK8CHS;
  jet_resolution_AK4CHS = JME::JetResolution(ak4_textfile_);
  //jet_resolution_AK8CHS = JME::JetResolution(ak8_textfile_);//ADDED:...AND ENDING HERE, NEW
  
  // jet resolution scale factor from text files
  JME::JetResolutionScaleFactor JER_scale_factor;
  JER_scale_factor = JME::JetResolutionScaleFactor(jet_scale_fac_textfile_);

#endif  
  // MARKER AK4 CHS JETS
  //
  //     _    _  ___  _      ____ _   _ ____        _ _____ _____ ____  
  //    / \  | |/ / || |    / ___| | | / ___|      | | ____|_   _/ ___| 
  //   / _ \ | ' /| || |_  | |   | |_| \___ \   _  | |  _|   | | \___ \ |
  //  / ___ \| . \|__   _| | |___|  _  |___) | | |_| | |___  | |  ___) |
  // /_/   \_\_|\_\  |_|    \____|_| |_|____/   \___/|_____| |_| |____/
  //

  /*
  int count_AK4 = 0;
  TLorentzVector Min_DeltaR_ak4_lep_4vector;
  TLorentzVector Min_DeltaR_btag_ak4_lep_4vector;
  
  double Min_DeltaR_ak4_lep_b_disc = -99.0;
  double Min_DeltaR_btag_ak4_lep_b_disc = -99.0;
  
  double Min_DeltaR_ak4_lep  = 99.0;
  double Min_DeltaR_btag_ak4_lep = 99.0;

  bool loose_tag  = false;
  bool middle_tag = false;
  bool tight_tag  = false;

  double ak4_baseline           = 0.0;
  double ak4_baseline_up    = 0.0;
  double ak4_baseline_down    = 0.0;
  double ak4_base_smear  = 0.0;
  double ak4_base_smear_up   = 0.0;
  double ak4_base_smear_down   = 0.0;

  edm::Handle<pat::JetCollection> AK4_MINI_JET;
  iEvent.getByToken(ak4jetToken_, AK4_MINI_JET);

  ///edm::Handle<reco::GenJetCollection> AK4GENJET;  
  ///iEvent.getByToken(ak4genjetToken_, AK4GENJET);


  for (const pat::Jet &ak4_iterator : *AK4_MINI_JET)
    { //NOTE: ak4_iterator is an iterator over ak4_mini_jet leading to AK4 dRMinLep 4 vector
       double iterator_eta  = ak4_iterator.eta();
    
      if (ak4_iterator.pt() < 30.0 || fabs(iterator_eta) > 10.0) {std::cout << "failed pt less than 30 or eta in range of 2.4\n"; continue;} //changing rapidity max from 2.4 to 10

      //------------------------------------
      // Noise jet ID
      //------------------------------------    
      double neutr_had_frac     = ak4_iterator.neutralHadronEnergyFraction();
      double neutr_em_frac      = ak4_iterator.neutralEmEnergyFraction();
      double chrg_had_frac      = ak4_iterator.chargedHadronEnergyFraction();
      double chrg_em_frac       = ak4_iterator.chargedEmEnergyFraction();
      double total_multiplicity = ak4_iterator.chargedMultiplicity() + ak4_iterator.neutralMultiplicity();
      double neutr_multiplicity = ak4_iterator.neutralMultiplicity();
      double chrg_multiplicity  = ak4_iterator.chargedMultiplicity(); 

      bool basic_cut =  
	(fabs(iterator_eta) <= 2.4 && neutr_had_frac < 0.99 && neutr_em_frac < 0.99 && total_multiplicity > 1.0 && chrg_had_frac > 0.0  && chrg_multiplicity > 0.0 && chrg_em_frac < 0.99) 
	|| (fabs(iterator_eta) <= 2.7 && fabs(iterator_eta) > 2.4 && neutr_had_frac < 0.99 && neutr_em_frac < 0.99 && total_multiplicity > 1.0) 
	|| (fabs(iterator_eta) <= 3.0 && fabs(iterator_eta) > 2.7 && neutr_em_frac < 0.9 && neutr_multiplicity > 2.0) 
	|| (fabs(iterator_eta)  > 3.0 && neutr_em_frac < 0.9 && neutr_multiplicity > 10.0);
      
      std::cout << "  goodJet = " << basic_cut << std::endl;

      if (!basic_cut) {std::cout << "  bad AK4 jet. skip. (pt " << ak4_iterator.pt() << " eta " << iterator_eta << " total multiplicity " << total_multiplicity << ")" << std::endl; continue;}
      
     
      
      //------------------------------------
      // AK4CHS JEC correction 
      //------------------------------------
      reco::Candidate::LorentzVector ak4_jet_uncorr = ak4_iterator.correctedP4(0);
      JetCorrectorAK4chs->setJetEta(ak4_jet_uncorr.eta());
      JetCorrectorAK4chs->setJetPt (ak4_jet_uncorr.pt());
      JetCorrectorAK4chs->setJetE  (ak4_jet_uncorr.energy());
      JetCorrectorAK4chs->setJetA  (ak4_iterator.jetArea());
      JetCorrectorAK4chs->setRho   (rho_var);
      JetCorrectorAK4chs->setNPV   (num_verticies);
      
      double correction_factor = JetCorrectorAK4chs->getCorrection();

      reco::Candidate::LorentzVector ak4_jet_corr = correction_factor * ak4_jet_uncorr;

      double reco_ak4_pt   = ak4_jet_corr.pt();
      double reco_ak4_phi  = ak4_jet_corr.phi();
      double reco_ak4_eta  = ak4_jet_corr.eta();
      double reco_ak4_mass = ak4_jet_corr.mass();
      
      std::cout << "uncorrected AK4 jet pt " << ak4_jet_uncorr.pt() << " corrected jet pt " << reco_ak4_pt << std::endl;
    
      //------------------------------------
      // AK4CHS JEC uncertainty
      //------------------------------------
      double downward_correction = 1.0;
      double upward_correction = 1.0;
      
      JetCorrUncertAK4chs->setJetPhi(reco_ak4_phi);
      JetCorrUncertAK4chs->setJetEta(reco_ak4_eta);     
      JetCorrUncertAK4chs->setJetPt (reco_ak4_pt);

      downward_correction = correction_factor - JetCorrUncertAK4chs->getUncertainty(0);
      upward_correction   = correction_factor + JetCorrUncertAK4chs->getUncertainty(1);

      std::cout << "  downward correction " << downward_correction << " upward correction " << upward_correction << std::endl;

      //------------------------------------
      // GenJet  matched
      //------------------------------------ 
      double gen_ak4_pt = 0.0;
      
      const reco::GenJet* gen_level_jetsource = ak4_iterator.genJet();
      if (gen_level_jetsource) {gen_ak4_pt = gen_level_jetsource->pt(); std::cout << "  ak4 gen level jet source pt " << gen_ak4_pt << " mass " << gen_level_jetsource->mass() << std::endl;}
       
      //------------------------------------
      // JER SF
      //------------------------------------
      double pt_smear      = 1.0;
      double pt_smear_up   = 1.0;
      double pt_smear_down = 1.0;
     
	
      double JER_scale_fact      = 2.0;//JER_scale_factor.getScaleFactor({{JME::Binning::JetEta, ak4_jet_corr.eta()}});
      double JER_up_scale_fact   = 2.0;//JER_scale_factor.getScaleFactor({{JME::Binning::JetEta, ak4_jet_corr.eta()}}, Variation::UP);
      double JER_down_scale_fact = 2.0;//JER_scale_factor.getScaleFactor({{JME::Binning::JetEta, ak4_jet_corr.eta()}}, Variation::DOWN);
      std::cout << "  JER Scale factors (Nominal / Up / Down) : " << JER_scale_fact << " / " << JER_up_scale_fact << " / " << JER_down_scale_fact << std::endl;
      
      // double gen_ak4_pt     = GenJetMatched.Perp();
      double Delta_pt      = (reco_ak4_pt - gen_ak4_pt) * (JER_scale_fact      - 1.0);
      double Delta_pt_up   = (reco_ak4_pt - gen_ak4_pt) * (JER_up_scale_fact   - 1.0);
      double Delta_pt_down = (reco_ak4_pt - gen_ak4_pt) * (JER_down_scale_fact - 1.0);
      pt_smear      = std::max((double)0.0, (reco_ak4_pt + Delta_pt)      / reco_ak4_pt);
      pt_smear_up   = std::max((double)0.0, (reco_ak4_pt + Delta_pt_up)   / reco_ak4_pt);
      pt_smear_down = std::max((double)0.0, (reco_ak4_pt + Delta_pt_down) / reco_ak4_pt);
	

      std::cout << "  pt smear " << pt_smear << " pt smear up " << pt_smear_up << " pt smear down " << pt_smear_down << std::endl;


      //------------------------------------
      // AK4 variables 
      //------------------------------------
      std::cout << "pt for jet4 is " << reco_ak4_pt << std::endl;
      std::cout << "mass for jet4 is " << reco_ak4_mass << std::endl;
      std::cout << "eta for jet4 is " << reco_ak4_eta << std::endl;
      std::cout << "phi for jet4 is " << reco_ak4_phi << std::endl;
      //double rapidity     = ak4_iterator.rapidity();
      //double ndau         = ak4_iterator.numberOfDaughters();
      double b_discriminator        = ak4_iterator.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
      //double puid         = ak4_iterator.userFloat("pileupJetId:fullDiscriminant");
 
      //------------------------------------
      // HT calculation
      //------------------------------------

      if (reco_ak4_pt > 30.0)                       ak4_baseline        +=   reco_ak4_pt;
      if (upward_correction * reco_ak4_pt > 30.0)   ak4_baseline_up     +=   upward_correction * reco_ak4_pt;
      if (downward_correction * reco_ak4_pt > 30.0) ak4_baseline_down   +=   downward_correction * reco_ak4_pt;
      if (pt_smear * reco_ak4_pt > 30.0)            ak4_base_smear      +=   pt_smear * reco_ak4_pt;
      if (pt_smear_up * reco_ak4_pt > 30.0)         ak4_base_smear_up   +=   pt_smear_up * reco_ak4_pt;
      if (pt_smear_down * reco_ak4_pt > 30.0)       ak4_base_smear_down +=   pt_smear_down * reco_ak4_pt;

      //------------------------------------
      // Find AK4 jet closest to lepton
      //------------------------------------ 
      double DeltaR_ak4_lep = deltaR(reco_ak4_eta, reco_ak4_phi, lepton.Eta(), lepton.Phi());

      if (reco_ak4_pt > 40.0 && fabs(reco_ak4_eta) < 2.4 && basic_cut && DeltaR_ak4_lep < Min_DeltaR_ak4_lep)
	{
	      Min_DeltaR_ak4_lep = DeltaR_ak4_lep;
	      Min_DeltaR_ak4_lep_4vector.SetPtEtaPhiM(reco_ak4_pt, reco_ak4_eta, reco_ak4_phi, reco_ak4_mass);
	      Min_DeltaR_ak4_lep_b_disc = b_discriminator;
	}

      //------------------------------------
      // Find Loose b-tagged AK4 jet closest to lepton
      //------------------------------------ 
      if (reco_ak4_pt > 40.0 && fabs(reco_ak4_eta) < 2.4 && basic_cut && b_discriminator > 0.460 && DeltaR_ak4_lep < Min_DeltaR_btag_ak4_lep)
	{
	      Min_DeltaR_btag_ak4_lep = DeltaR_ak4_lep;
	      Min_DeltaR_btag_ak4_lep_4vector.SetPtEtaPhiM(reco_ak4_pt, reco_ak4_eta, reco_ak4_phi, reco_ak4_mass); 
	      Min_DeltaR_btag_ak4_lep_b_disc = b_discriminator;
	}

      //------------------------------------
      // Check if there is a b-tagged AK4 jet in the lepton hemisphere
      //------------------------------------   
      if (reco_ak4_pt > 40.0 && fabs(reco_ak4_eta) < 2.4 && basic_cut && fabs(deltaPhi(reco_ak4_phi, lepton.Phi())) < 3.14 / 2.0)
	{             
	      if (b_discriminator > 0.460) loose_tag  = true;
	      if (b_discriminator > 0.800) middle_tag = true;
	      if (b_discriminator > 0.935) tight_tag  = true;
	}
    } //end AK4 loop

      std::cout << "AK4 summary:" << std::endl;
      std::cout << "  closest ak4 jet to lepton:" << std::endl;
      std::cout << "    pt =  " << Min_DeltaR_ak4_lep_4vector.Perp() << std::endl;
      std::cout << "    b discriminator =  " << Min_DeltaR_ak4_lep_b_disc << std::endl;
      std::cout << "    dR  = " << Min_DeltaR_ak4_lep << std::endl;
      std::cout << "  closest loose b-tagged ak4 jet to lepton:" << std::endl;
      std::cout << "    pt =  " << Min_DeltaR_btag_ak4_lep_4vector.Perp() << std::endl;
      std::cout << "    b discriminator =  " << Min_DeltaR_btag_ak4_lep_b_disc << std::endl;
      std::cout << "    dR  = " << Min_DeltaR_btag_ak4_lep << std::endl;
      std::cout << "  b-tagged jet in hemisphere around lepton?" << std::endl;
      std::cout << "    loose_tag  " << loose_tag  << std::endl; 
      std::cout << "    middle_tag " << middle_tag << std::endl; 
      std::cout << "    tight_tag  " << tight_tag  << std::endl; 

  */

  // MARKER AK4 CHS JETS
  //
  //     _    _  ___  _      ____ _   _ ____        _ _____ _____ ____  
  //    / \  | |/ / || |    / ___| | | / ___|      | | ____|_   _/ ___| 
  //   / _ \ | ' /| || |_  | |   | |_| \___ \   _  | |  _|   | | \___ \ |
  //  / ___ \| . \|__   _| | |___|  _  |___) | | |_| | |___  | |  ___) |
  // /_/   \_\_|\_\  |_|    \____|_| |_|____/   \___/|_____| |_| |____/
  //
#if 0  
  int count_AK4 = 0;
  TLorentzVector AK4_dRMinLep_p4;
  TLorentzVector AK4_btagged_dRMinLep_p4;
  
  double AK4_dRMinLep_bdisc = -99; AK4_dRMinLep_bdisc = AK4_dRMinLep_bdisc;
  double AK4_btagged_dRMinLep_bdisc = -99;
  
  double AK4_dRMinLep_deltaR  = 99; AK4_dRMinLep_deltaR = AK4_dRMinLep_deltaR;
  double AK4_btagged_dRMinLep = 99;

  double AK4_dRMinLep_ptsmear    = 1; AK4_dRMinLep_ptsmear = AK4_dRMinLep_ptsmear;
  double AK4_dRMinLep_ptsmearUp  = 1; AK4_dRMinLep_ptsmearUp = AK4_dRMinLep_ptsmearUp;
  double AK4_dRMinLep_ptsmearDn  = 1; AK4_dRMinLep_ptsmearDn = AK4_dRMinLep_ptsmearDn;
  double AK4_dRMinLep_ptuncorr   = 0; AK4_dRMinLep_ptuncorr = AK4_dRMinLep_ptuncorr;
  double AK4_dRMinLep_corr       = 1; AK4_dRMinLep_corr = AK4_dRMinLep_corr;
  double AK4_dRMinLep_corrUp     = 1; AK4_dRMinLep_corrUp = AK4_dRMinLep_corrUp;
  double AK4_dRMinLep_corrDn     = 1; AK4_dRMinLep_corrDn = AK4_dRMinLep_corrDn;

  bool ak4_btag_loose  = false;
  bool ak4_btag_medium = false;
  bool ak4_btag_tight  = false;

  double HT_AK4_pt30           = 0;
  double HT_AK4_pt30_corrUp    = 0;
  double HT_AK4_pt30_corrDn    = 0;
  double HT_AK4_pt30_smearNom  = 0;
  double HT_AK4_pt30_smearUp   = 0;
  double HT_AK4_pt30_smearDn   = 0;

  edm::Handle<pat::JetCollection> AK4MINI;
  iEvent.getByToken(ak4jetToken_, AK4MINI);

  edm::Handle<reco::GenJetCollection> AK4GENJET;  
  iEvent.getByToken(ak4genjetToken_, AK4GENJET);


  //if (verbose_) std::cout << "AK4 jet loop" << std::endl;

  for (const pat::Jet &ijet : *AK4MINI)
    {   
      if (ijet.pt() < 15.0 || fabs(ijet.eta()) > 3.0) continue; 
      //if (verbose_) std::cout << " jet " << count_AK4 << std::endl;
      count_AK4++;

      //------------------------------------
      // Remove leptons from AK4 jets
      //------------------------------------    
      auto uncorrJetObj = ijet.correctedJet(0);
      reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
      // now loop on pf candidates
      //// Jet constituent indices for lepton matching
      std::vector<int> constituentIndices;
      auto jetConstituents = uncorrJetObj.daughterPtrVector();
      //if (verbose_) std::cout << "   -> before lepton cleaning uncorr pt,eta,phi,m = " << uncorrJet.pt() << ", " << uncorrJet.eta() << ", " << uncorrJet.phi() << ", " << uncorrJet.mass() << std::endl;
      /* for ( auto & constituent : jetConstituents)
	 {
	 // If this constituent is part of a muon, remove the constituent's four vector
	 if (std::find(muFootprint.begin(), muFootprint.end(), constituent ) != muFootprint.end())
	 {
	 uncorrJet -= constituent->p4();
	 //if (verbose_) std::cout << "     -> REMOVED part of muon" << std::endl;
	 }
	 // If this constituent is part of an electron, remove the constituent's four vector
	 if (std::find(elFootprint.begin(), elFootprint.end(), constituent ) != elFootprint.end())
	 {
	 uncorrJet -= constituent->p4();
	 //if (verbose_) std::cout << "     -> REMOVED part of electron" << std::endl;
	 }
	 }*/
      //if (verbose_) std::cout << "   -> after lepton cleaning uncorr pt,eta,phi,m = " << uncorrJet.pt() << ", " << uncorrJet.eta() << ", " << uncorrJet.phi() << ", " << uncorrJet.mass() << std::endl;

      //------------------------------------
      // Noise jet ID
      //------------------------------------    

      double NHF       = ijet.neutralHadronEnergyFraction();
      double NEMF      = ijet.neutralEmEnergyFraction();
      double CHF       = ijet.chargedHadronEnergyFraction();
      // double MUF       = ijet.muonEnergyFraction();
      double CEMF      = ijet.chargedEmEnergyFraction();
      double NumConst  = ijet.chargedMultiplicity() + ijet.neutralMultiplicity();
      double NM        = ijet.neutralMultiplicity();
      double CM        = ijet.chargedMultiplicity(); 
      double eta       = ijet.eta(); 

      bool goodJet_looseJetID =  
	(fabs(eta) <= 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst >1 && CHF > 0.0 && CM > 0 && CEMF < 0.99) 
	|| (fabs(eta) <= 2.7 && fabs(eta) > 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst > 1) 
	|| (fabs(eta) <= 3.0 && fabs(eta) > 2.7 && NHF < 0.98 && NEMF > 0.01 && NM > 2) 
	|| (fabs(eta)  > 3.0 && NEMF < 0.9 && NM > 10);
      //if (verbose_ && goodJet_looseJetID) std::cout << "   -> goodJet " << std::endl;

      if (!goodJet_looseJetID)
	{
	  //if(verbose_) std::cout << "   -> bad AK4 jet. skip.  ( pt " << ijet.pt() << " eta " << ijet.eta() << " NumConst " << NumConst << " )" << std::endl;
	  continue;
	}

      //------------------------------------
      // AK4CHS JEC correction 
      //------------------------------------
      JetCorrectorAK4chs->setJetEta(uncorrJet.eta());
      JetCorrectorAK4chs->setJetPt (uncorrJet.pt());
      JetCorrectorAK4chs->setJetE  (uncorrJet.energy());
      JetCorrectorAK4chs->setJetA  (ijet.jetArea());
      JetCorrectorAK4chs->setRho   (rho_var);
      JetCorrectorAK4chs->setNPV   (num_verticies);
      double corr = JetCorrectorAK4chs->getCorrection();

      reco::Candidate::LorentzVector corrJet = corr * uncorrJet;
      //if (verbose_) std::cout << "   -> after JEC pt,eta,phi,m = " << corrJet.pt() << ", " << corrJet.eta() << ", " << corrJet.phi() << ", " << corrJet.mass() << std::endl;
    
      if (corrJet.pt() < 15) continue;  

      //------------------------------------
      // AK4CHS JEC uncertainty
      //------------------------------------
      double corrDn = 1.0;
      JetCorrUncertAK4chs->setJetPhi(corrJet.phi());
      JetCorrUncertAK4chs->setJetEta(corrJet.eta());
      JetCorrUncertAK4chs->setJetPt(corrJet.pt());
      corrDn = corr - JetCorrUncertAK4chs->getUncertainty(0);
      double corrUp = 1.0;
      JetCorrUncertAK4chs->setJetPhi(corrJet.phi());
      JetCorrUncertAK4chs->setJetEta(corrJet.eta());
      JetCorrUncertAK4chs->setJetPt(corrJet.pt());
      corrUp = corr + JetCorrUncertAK4chs->getUncertainty(1);

      //if (verbose_) std::cout << "   -> corr " << corr << " corrDn " << corrDn << " corrUp " << corrUp << std::endl;

      //------------------------------------
      // AK4 JER SF
      //------------------------------------
   
      double ptsmear   = 1;
      double ptsmearUp = 1;
      double ptsmearDn = 1;

      if (!iEvent.isRealData())
	{
	  //if (verbose_) std::cout << "   Get JER SF" << std::endl;

	  // get genjet
	  double genpt = 0;
	  TLorentzVector GenJetMatched;
	  const reco::GenJet* genJet = ijet.genJet();
	  //bool foundgenjet = false;
	  if (genJet)
	    {
	      //foundgenjet = true;
	      genpt = genJet->pt();
	      GenJetMatched.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass());
	      //if (verbose_) std::cout << "      -> Found ak4 genJet pt " << genJet->pt() << " mass " << genJet->mass() << std::endl;
	    }
	  //else if (verbose_) std::cout << "      -> Did not find genJet" << std::endl;
    
	  // Set parameters needed for jet resolution and scale factors
	  JME::JetParameters jer_parameters;
	  jer_parameters.setJetPt (corrJet.pt());
	  jer_parameters.setJetEta(corrJet.eta());
	  jer_parameters.setRho   (rho_var);

	  // Get resolution
	  double res = jet_resolution_AK4CHS.getResolution(jer_parameters); 

	  // Get scale factors
	  double jer_sf    = JER_scale_factor.getScaleFactor(jer_parameters);
	  double jer_sf_up = JER_scale_factor.getScaleFactor(jer_parameters, Variation::UP);
	  double jer_sf_dn = JER_scale_factor.getScaleFactor(jer_parameters, Variation::DOWN);
	  //if (verbose_) std::cout << "      -> JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << "    & Resolution :" << res << std::endl;
     
	  // Get Smearings  
	  // --- If well matched, smear based on GenJet, If not well matched,  gaussian smear based on resolution
	  TLorentzVector AK4JetP4;
	  AK4JetP4.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());
	  double DeltaR_gen_reco  = AK4JetP4.DeltaR( GenJetMatched );
	  double DeltaPt_gen_reco = AK4JetP4.Pt() - GenJetMatched.Pt();
	  double jet_distance_param = 0.4; 
	  //if (verbose_) std::cout << "      -> gen pt " << GenJetMatched.Pt() << " reco pt " << AK4JetP4.Pt() << "  delta " << DeltaPt_gen_reco << std::endl;

	  if (genJet && (DeltaR_gen_reco<jet_distance_param/2.0) && (std::abs(DeltaPt_gen_reco)<(3*res*AK4JetP4.Pt())))
	    {
	      //if (verbose_) std::cout << "      -> Well matched (recojet,genjet)" << std::endl;
	      double recopt    = corrJet.pt();
	      // double genpt     = GenJetMatched.Perp();
	      double deltapt   = (recopt-genpt)*(jer_sf - 1.0);
	      double deltaptUp = (recopt-genpt)*(jer_sf_up - 1.0);
	      double deltaptDn = (recopt-genpt)*(jer_sf_dn - 1.0);

	      ptsmear   = std::max((double)0.0, (recopt + deltapt)/recopt     );
	      ptsmearUp = std::max((double)0.0, (recopt + deltaptUp)/recopt   );
	      ptsmearDn = std::max((double)0.0, (recopt + deltaptDn)/recopt   );
	    }
	  else
	    {
	      /*if (verbose_)
		{
		std::cout << "      -> Not well matched. DeltaR_gen_reco " << DeltaR_gen_reco << " DeltaPt_gen_reco " << DeltaPt_gen_reco << " 3*res*AK4JetP4.Pt()) " << 3*res*AK4JetP4.Pt();
		if (!foundgenjet) std::cout << ". Did not find genjet" << std::endl;
		else std::cout << std::endl;
		}*/
	      double sigma   = std::sqrt(jer_sf * jer_sf - 1)       * res;  
	      double sigmaUp = std::sqrt(jer_sf_up * jer_sf_up - 1) * res;
	      double sigmaDn = std::sqrt(jer_sf_dn * jer_sf_dn - 1) * res;

	      TRandom3 rand1(0);
	      ptsmear   = std::max( (double)0.0, 1 + rand1.Gaus(0, sigma  ) );
	      ptsmearUp = std::max( (double)0.0, 1 + rand1.Gaus(0, sigmaUp) );
	      ptsmearDn = std::max( (double)0.0, 1 + rand1.Gaus(0, sigmaDn) );
	    }
	}

      //if (verbose_) std::cout << "   -> ptsmear " << ptsmear << " ptsmearUp " <<  ptsmearDn << " ptsmearDn " <<  ptsmearUp << std::endl;

      //------------------------------------
      // AK4 variables 
      //------------------------------------
      double pt           = corrJet.pt();
      double mass         = corrJet.mass();
      eta                 = corrJet.eta();
      double phi_ak4      = corrJet.phi();
      //double rapidity     = ijet.rapidity();
      //double ndau         = ijet.numberOfDaughters();
      double bdisc        = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
      //double puid         = ijet.userFloat("pileupJetId:fullDiscriminant");
 
      //  if (corrJet.pt()>20 ){
      //   vAK4pt    ->push_back( corrJet.pt()    ); 
      //   vAK4eta   ->push_back( corrJet.eta()   ); 
      //   vAK4phi   ->push_back( corrJet.phi()   ); 
      //   vAK4m     ->push_back( corrJet.mass()  ); 
      //   vAK4bdisc ->push_back( bdisc           );

      //   SemiLeptAK4pt    ->push_back( corrJet.pt()    ); 
      //   SemiLeptAK4eta   ->push_back( corrJet.eta()   ); 
      //   SemiLeptAK4phi   ->push_back( corrJet.phi()   ); 
      //   SemiLeptAK4m     ->push_back( corrJet.mass()  ); 
      //   SemiLeptAK4bdisc ->push_back( bdisc           );
      // }  

      //------------------------------------
      // HT calculation
      //------------------------------------

      if (corrJet.pt() > 30)            HT_AK4_pt30           +=   pt;
      if (corrUp * corrJet.pt() > 30)   HT_AK4_pt30_corrUp    +=   corrUp * uncorrJet.pt();
      if (corrDn * corrJet.pt() > 30)   HT_AK4_pt30_corrDn    +=   corrDn * uncorrJet.pt();
      if (ptsmear * corrJet.pt() > 30)  HT_AK4_pt30_smearNom  +=   ptsmear * corrJet.pt();
      if (ptsmearUp * corrJet.pt() > 30) HT_AK4_pt30_smearUp   +=   ptsmearUp * corrJet.pt();
      if (ptsmearDn * corrJet.pt() > 30) HT_AK4_pt30_smearDn   +=   ptsmearDn * corrJet.pt();

      //------------------------------------
      // Find AK4 jet closest to lepton
      //------------------------------------ 
      double deltaRlep = deltaR(corrJet.eta(), corrJet.phi(), lepton.Eta(), lepton.Phi() );

      if (pt > 25.0 && fabs(eta) < 3.0 && goodJet_looseJetID)
	{
	  if (deltaRlep<AK4_dRMinLep_deltaR )
	    {
	      AK4_dRMinLep_deltaR = deltaRlep;
	      AK4_dRMinLep_p4.SetPtEtaPhiM( pt, eta, phi_ak4, mass );
	      AK4_dRMinLep_bdisc = bdisc;
	      AK4_dRMinLep_ptsmear   = ptsmear;
	      AK4_dRMinLep_ptsmearUp = ptsmearUp;
	      AK4_dRMinLep_ptsmearDn = ptsmearDn;
	      AK4_dRMinLep_ptuncorr  = uncorrJet.pt();
	      AK4_dRMinLep_corr    = corr;
	      AK4_dRMinLep_corrUp  = corrUp;
	      AK4_dRMinLep_corrDn  = corrDn;
	    }
	}

      //------------------------------------
      // Find Loose b-tagged AK4 jet closest to lepton
      //------------------------------------ 
      if (pt > 25.0 && fabs(eta) < 3.0 && goodJet_looseJetID && bdisc > 0.5426)
      {
	if (deltaRlep<AK4_btagged_dRMinLep)
	  {
	    AK4_btagged_dRMinLep = deltaRlep;
	    AK4_btagged_dRMinLep_p4.SetPtEtaPhiM(pt, eta, phi_ak4, mass);
	    AK4_btagged_dRMinLep_bdisc = bdisc;
	  }
      }

      //------------------------------------
      // Check if there is a b-tagged AK4 jet in the lepton hemisphere
      //------------------------------------ 
      double deltaPhiLep = fabs( deltaPhi(phi,  lepton.Phi()));  
      if (pt > 25.0 && fabs(eta) < 3.0 && goodJet_looseJetID)
	{              
	  if (deltaPhiLep<  3.14/2.0)
	    {
	      if (bdisc > 0.5426) ak4_btag_loose  = true;
	      if (bdisc > 0.8484) ak4_btag_medium = true;
	      if (bdisc > 0.9535) ak4_btag_tight  = true;
	    }
	}
    } //end AK4 loop

  if (true)
    {
      std::cout << "AK4 summary:" << std::endl;
      std::cout << "  closest ak4 jet to lepton:" << std::endl;
      std::cout << "    pt =  " << AK4_dRMinLep_p4.Perp() << std::endl;
      std::cout << "    bdisc =  " << AK4_dRMinLep_bdisc << std::endl;
      std::cout << "    dR  = " << AK4_dRMinLep_deltaR << std::endl;
      std::cout << "  closest loose b-tagged ak4 jet to lepton:" << std::endl;
      std::cout << "    pt =  " << AK4_btagged_dRMinLep_p4.Perp() << std::endl;
      std::cout << "    bdisc =  " << AK4_btagged_dRMinLep_bdisc << std::endl;
      std::cout << "    dR  = " << AK4_btagged_dRMinLep << std::endl;
      std::cout << "  b-tagged jet in hemisphere around lepton?" << std::endl;
      std::cout << "    ak4_btag_loose  " << ak4_btag_loose   << std::endl; 
      std::cout << "    ak4_btag_medium " << ak4_btag_medium  << std::endl; 
      std::cout << "    ak4_btag_tight  " << ak4_btag_tight   << std::endl; 

      std::cout << "HT_AK4_pt30          " << HT_AK4_pt30          << std::endl;
      std::cout << "HT_AK4_pt30_corrUp   " << HT_AK4_pt30_corrUp   << std::endl;
      std::cout << "HT_AK4_pt30_corrDn   " << HT_AK4_pt30_corrDn   << std::endl;
      std::cout << "HT_AK4_pt30_smearNom " << HT_AK4_pt30_smearNom << std::endl;
      std::cout << "HT_AK4_pt30_smearUp  " << HT_AK4_pt30_smearUp  << std::endl;
      std::cout << "HT_AK4_pt30_smearDn  " << HT_AK4_pt30_smearDn  << std::endl;
    }
#endif


  
  // MARKER ANGLES
  //    ##    #    #   ####   #       ######   ####
  //   #  #   ##   #  #    #  #       #       #
  //  #    #  # #  #  #       #       #####    ####
  //  ######  #  # #  #  ###  #       #            #
  //  #    #  #   ##  #    #  #       #       #    #
  //  #    #  #    #   ####   ######  ######   ####
  //if (mu_event == 1 && el_event == 0)
  //{

  if (bad_event) return; //{std::cout << "bad event reached angles\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
  
  if (reco_lep_generation > 0) calculateAngles(lepton, met4vector, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);//changing temporarily to gen
  if (reco_lep_generation < 0) calculateAngles(met4vector, lepton, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);   

  if (reco_lep_generation > 0 && angle_test) calculateAngles(lepton_gen, met4vector, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);//changing temporarily to gen
  if (reco_lep_generation < 0 && angle_test) calculateAngles(met4vector, lepton_gen, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);   

  costheta1 = d_costheta1;
  costheta2 = d_costheta2;
  costhetastar = d_costhetastar;
  phi = d_phi;
  phi1 = d_phi1;
  phi2 = d_phi2;
      
  if (reco_lep_generation > 0) calculateAngles(lepton, other_met4vector, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
  if (reco_lep_generation < 0) calculateAngles(other_met4vector, lepton, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);

  if (reco_lep_generation > 0 && angle_test) calculateAngles(lepton_gen, other_met4vector, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
  if (reco_lep_generation < 0 && angle_test) calculateAngles(other_met4vector, lepton_gen, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);

  other_costheta1 = d_costheta1;
  other_costheta2 = d_costheta2;
  other_costhetastar = d_costhetastar;
  other_phi = d_phi;
  other_phi1 = d_phi1;
  other_phi2 = d_phi2;
  //}
  /*
    else if (mu_event == 0 && el_event == 1)
    {
    if (reco_lep_generation > 0) calculateAngles(electron, met4vector, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
    if (reco_lep_generation < 0) calculateAngles(met4vector, electron, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
    
    costheta1 = d_costheta1;
    costheta2 = d_costheta2;
    costhetastar = d_costhetastar;
    phi = d_phi;
    phi1 = d_phi1;
    phi2 = d_phi2;
      
    if (reco_lep_generation > 0) calculateAngles(electron, other_met4vector, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
    if (reco_lep_generation < 0) calculateAngles(other_met4vector, electron, subjet0, subjet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);  

    other_costheta1 = d_costheta1;
    other_costheta2 = d_costheta2;
    other_costhetastar = d_costhetastar;
    other_phi = d_phi;
    other_phi1 = d_phi1;
    other_phi2 = d_phi2;
    }
      
    else
    {
    std::cout << "neither muon nor electron event.\n";
    bad_event = true;
    }
  */
  if (reco_lep_generation == 0 || reco_lep_generation < -2 || reco_lep_generation > 2)
    {
      costheta1 = 90.0;
      costheta2 = 90.0;
      costhetastar = 90.0;
      phi = 90.0;
      phi1 = 90.0;
      phi2 = 90.0;
    }

  if (costheta1 != costheta1) costheta1 = 2.0;
  if (costheta2 != costheta2) costheta2 = 2.0;
  if (costhetastar != costhetastar) costhetastar = 2.0;
  if (phi != phi) phi = 4.0;
  if (phi1 != phi1) phi1 = 4.0;
  if (phi2 != phi2) phi2 = 4.0;

  if (other_costheta1 != other_costheta1) other_costheta1 = 2.0;
  if (other_costheta2 != other_costheta2) other_costheta2 = 2.0;
  if (other_costhetastar != other_costhetastar) other_costhetastar = 2.0;
  if (other_phi != other_phi) other_phi = 4.0;
  if (other_phi1 != other_phi1) other_phi1 = 4.0;
  if (other_phi2 != other_phi2) other_phi2 = 4.0;

  if (gen_lep_generation  > 0) calculateAngles(lepton_gen, neutrino_gen, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
  if (gen_lep_generation  < 0) calculateAngles(neutrino_gen, lepton_gen, quark_gen, antiquark_gen, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);

  costheta1_gen = d_costheta1;
  costheta2_gen = d_costheta2;
  phi_gen = d_phi;
  costhetastar_gen = d_costhetastar;
  phi1_gen = d_phi1;
  phi2_gen = d_phi2;

   if (gen_lep_generation == 0)
    {
      costheta1_gen = 90.0;
      costheta2_gen = 90.0;
      costhetastar_gen = 90.0;
      phi_gen = 90.0;
      phi1_gen = 90.0;
      phi2_gen = 90.0;
    }

  /*
  if (phi != phi) {std::cout << "undefined phi.\n"; bad_event = true;}
  if (phi2 != phi2) {std::cout << "undefined phi1.\n"; bad_event = true;}
  if (phi1 != phi1) {std::cout << "undefined phi2.\n"; bad_event = true;}
  if (costheta1 != costheta1) {std::cout << "undefined costheta1.\n"; bad_event = true;}
  if (costheta2 != costheta2) {std::cout << "undefined costheta2.\n"; bad_event = true;}
  if (costhetastar != costhetastar) {std::cout << "undefined costhetastar.\n"; bad_event = true;}
  */
  //if (phi == phi && phi2 == phi2 && phi1 == phi1 && costheta1 == costheta1 && costheta2 == costheta2 && costhetastar == costhetastar) bad_event = true;
    
   if (costheta1 == 2.0 || costheta2 == 2.0 || costhetastar == 2.0 || phi == 4.0 || phi1 == 4.0 || phi2 == 4.0) bad_event = true;

  // 
  // MARKER INTERMEDIATE
  //   _       _                               _ _       _       
  //  (_)_ __ | |_ ___ _ __ _ __ ___   ___  __| (_) __ _| |_ ___ 
  //  | | '_ \| __/ _ \ '__| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ |
  //  | | | | | ||  __/ |  | | | | | |  __/ (_| | | (_| | ||  __/
  //  |_|_| |_|\__\___|_|  |_| |_| |_|\___|\__,_|_|\__,_|\__\___|
  //
    
  
  //if (mu_event == 1 && el_event == 0)
  //{
  if (reco_lep_generation > 0) intermediate_steps(lepton, met4vector, subjet0, subjet1, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z);
      
  if (reco_lep_generation < 0) intermediate_steps(met4vector, lepton, subjet0, subjet1, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z);
   
  leptons_in_lep_px = d_leptons_in_lep_px;
  leptons_in_lep_py = d_leptons_in_lep_py;
  leptons_in_lep_pz = d_leptons_in_lep_pz;
  partons_in_lep_px = d_partons_in_lep_px;
  partons_in_lep_py = d_partons_in_lep_py;
  partons_in_lep_pz = d_partons_in_lep_pz;
  parton1_in_lep_px = d_parton1_in_lep_px;
  parton2_in_lep_px = d_parton2_in_lep_px;
  parton1_in_lep_py = d_parton1_in_lep_py;
  parton2_in_lep_py = d_parton2_in_lep_py;
  parton1_in_lep_pz = d_parton1_in_lep_pz;
  parton2_in_lep_pz = d_parton2_in_lep_pz;
  lepton1_in_lep_px = d_lepton1_in_lep_px;
  lepton1_in_lep_py = d_lepton1_in_lep_py;
  lepton1_in_lep_pz = d_lepton1_in_lep_pz;
  lepton1_dotted_x = d_lepton1_dotted_x;
  lepton1_dotted_y = d_lepton1_dotted_y;
  lepton1_dotted_z = d_lepton1_dotted_z;
  leptons_in_had_px = d_leptons_in_had_px;
  leptons_in_had_py = d_leptons_in_had_py;
  leptons_in_had_pz = d_leptons_in_had_pz;
  lepton1_in_had_px = d_lepton1_in_had_px;
  lepton1_in_had_py = d_lepton1_in_had_py;
  lepton1_in_had_pz = d_lepton1_in_had_pz;
  lepton2_in_had_px = d_lepton2_in_had_px;
  lepton2_in_had_py = d_lepton2_in_had_py;
  lepton2_in_had_pz = d_lepton2_in_had_pz;
  parton1_in_had_px = d_parton1_in_had_px;
  parton1_in_had_py = d_parton1_in_had_py;
  parton1_in_had_pz = d_parton1_in_had_pz;
  parton1_dotted_x = d_parton1_dotted_x;
  parton1_dotted_y = d_parton1_dotted_y;
  parton1_dotted_z = d_parton1_dotted_z;
  complicated1_px = d_complicated1_px;
  complicated1_py = d_complicated1_py;
  complicated1_pz = d_complicated1_pz;
  complicated2_px = d_complicated2_px;
  complicated2_py = d_complicated2_py;
  complicated2_pz = d_complicated2_pz;
  lepton_sumWWframe_X = d_lepton_sumWWframe_X;
  lepton_sumWWframe_Y = d_lepton_sumWWframe_Y;
  lepton_sumWWframe_Z = d_lepton_sumWWframe_Z;
  lepton1WWframe_X = d_lepton1WWframe_X;
  lepton1WWframe_Y = d_lepton1WWframe_Y;
  lepton1WWframe_Z = d_lepton1WWframe_Z;
  parton_sumWWframe_X = d_parton_sumWWframe_X;
  parton_sumWWframe_Y = d_parton_sumWWframe_Y;
  parton_sumWWframe_Z = d_parton_sumWWframe_Z;
  parton1WWframe_X = d_parton1WWframe_X;
  parton1WWframe_Y = d_parton1WWframe_Y;
  parton1WWframe_Z = d_parton1WWframe_Z;
  /*
    costheta1 = d_costheta1;
    costheta2 = d_costheta2;
    costhetastar = d_costhetastar;
    phi = d_phi;
    phi1 = d_phi1;
    phi2 = d_phi2;
  */
  boostWWframe_X = d_boostWWframe_X;
  boostWWframe_Y = d_boostWWframe_Y;
  boostWWframe_Z = d_boostWWframe_Z;
  boostWlep_X = d_boostWlep_X;
  boostWlep_Y = d_boostWlep_Y;
  boostWlep_Z = d_boostWlep_Z;
  boostWhad_X = d_boostWhad_X;
  boostWhad_Y = d_boostWhad_Y;
  boostWhad_Z = d_boostWhad_Z;

    
  //}
  /*
    else if (mu_event == 0 && el_event == 1)
    {
    if (reco_lep_generation > 0) intermediate_steps(electron, met4vector, subjet0, subjet1, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z);
      
    if (reco_lep_generation < 0) intermediate_steps(met4vector, electron, subjet0, subjet1, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z);
   
    leptons_in_lep_px = d_leptons_in_lep_px;
    leptons_in_lep_py = d_leptons_in_lep_py;
    leptons_in_lep_pz = d_leptons_in_lep_pz;
    partons_in_lep_px = d_partons_in_lep_px;
    partons_in_lep_py = d_partons_in_lep_py;
    partons_in_lep_pz = d_partons_in_lep_pz;
    parton1_in_lep_px = d_parton1_in_lep_px;
    parton2_in_lep_px = d_parton2_in_lep_px;
    parton1_in_lep_py = d_parton1_in_lep_py;
    parton2_in_lep_py = d_parton2_in_lep_py;
    parton1_in_lep_pz = d_parton1_in_lep_pz;
    parton2_in_lep_pz = d_parton2_in_lep_pz;
    lepton1_in_lep_px = d_lepton1_in_lep_px;
    lepton1_in_lep_py = d_lepton1_in_lep_py;
    lepton1_in_lep_pz = d_lepton1_in_lep_pz;
    lepton1_dotted_x = d_lepton1_dotted_x;
    lepton1_dotted_y = d_lepton1_dotted_y;
    lepton1_dotted_z = d_lepton1_dotted_z;
    leptons_in_had_px = d_leptons_in_had_px;
    leptons_in_had_py = d_leptons_in_had_py;
    leptons_in_had_pz = d_leptons_in_had_pz;
    lepton1_in_had_px = d_lepton1_in_had_px;
    lepton1_in_had_py = d_lepton1_in_had_py;
    lepton1_in_had_pz = d_lepton1_in_had_pz;
    lepton2_in_had_px = d_lepton2_in_had_px;
    lepton2_in_had_py = d_lepton2_in_had_py;
    lepton2_in_had_pz = d_lepton2_in_had_pz;
    parton1_in_had_px = d_parton1_in_had_px;
    parton1_in_had_py = d_parton1_in_had_py;
    parton1_in_had_pz = d_parton1_in_had_pz;
    parton1_dotted_x = d_parton1_dotted_x;
    parton1_dotted_y = d_parton1_dotted_y;
    parton1_dotted_z = d_parton1_dotted_z;
    complicated1_px = d_complicated1_px; 
    complicated1_py = d_complicated1_py;
    complicated1_pz = d_complicated1_pz;
    complicated2_px = d_complicated2_px;
    complicated2_py = d_complicated2_py; 
    complicated2_pz = d_complicated2_pz;
    lepton_sumWWframe_X = d_lepton_sumWWframe_X;
    lepton_sumWWframe_Y = d_lepton_sumWWframe_Y;
    lepton_sumWWframe_Z = d_lepton_sumWWframe_Z;
    lepton1WWframe_X = d_lepton1WWframe_X;
    lepton1WWframe_Y = d_lepton1WWframe_Y;
    lepton1WWframe_Z = d_lepton1WWframe_Z;
    parton_sumWWframe_X = d_parton_sumWWframe_X;
    parton_sumWWframe_Y = d_parton_sumWWframe_Y;
    parton_sumWWframe_Z = d_parton_sumWWframe_Z;
    parton1WWframe_X = d_parton1WWframe_X;
    parton1WWframe_Y = d_parton1WWframe_Y;
    parton1WWframe_Z = d_parton1WWframe_Z;
    / *
    costheta1 = d_costheta1;
    costheta2 = d_costheta2;
    costhetastar = d_costhetastar;
    phi = d_phi;
    phi1 = d_phi1;
    phi2 = d_phi2;
    * /
    boostWWframe_X = d_boostWWframe_X;
    boostWWframe_Y = d_boostWWframe_Y;
    boostWWframe_Z = d_boostWWframe_Z;
    boostWlep_X = d_boostWlep_X;
    boostWlep_Y = d_boostWlep_Y;
    boostWlep_Z = d_boostWlep_Z;
    boostWhad_X = d_boostWhad_X;
    boostWhad_Y = d_boostWhad_Y;
    boostWhad_Z = d_boostWhad_Z;
      
    //if (reco_lep_generation > 0) intermediate_steps(electron, met4vector, subjet0, subjet1, );
    //if (reco_lep_generation < 0) intermediate_steps(met4vector, electron, subjet0, subjet1, );

  
    }
  */

  // MARKER TRUTH MATCHING
  //   _____           _   _       __  __       _       _     _             
  //  |_   _| __ _   _| |_| |__   |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
  //    | || '__| | | | __| '_ \  | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
  //    | || |  | |_| | |_| | | | | |  | | (_| | || (__| | | | | | | | (_| |
  //    |_||_|   \__,_|\__|_| |_| |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
  //                                                                  |___/ 

  if (bad_event) return; //{std::cout << "bad event reached truth matching.\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
   
  lepton_deltaR = lepton.DeltaR(lepton_gen);
  lepton_deltaPhi = lepton.DeltaPhi(lepton_gen);
  
  neutrino_deltaR = met4vector.DeltaR(neutrino_gen);
  neutrino_deltaPhi = met4vector.DeltaPhi(neutrino_gen);
  neutrino_deltaPt = met4vector.Pt() - neutrino_gen.Pt();

  other_neutrino_deltaR = other_met4vector.DeltaR(neutrino_gen);
  other_neutrino_deltaPhi = other_met4vector.DeltaPhi(neutrino_gen);

  TLorentzVector truth_met4vector, false_met4vector;

  if (neutrino_deltaR > other_neutrino_deltaR)
    {
      truth_met4vector = other_met4vector;
      truth_neutrino_deltaR = other_neutrino_deltaR;
      truth_costheta1 = other_costheta1;
      truth_costheta2 = other_costheta2;
      truth_costhetastar = other_costhetastar;
      truth_phi = other_phi;
      truth_phi1 = other_phi1;
      truth_phi2 = other_phi2;

      false_met4vector = met4vector;
      false_neutrino_deltaR = neutrino_deltaR;
      false_costheta1 = costheta1;
      false_costheta2 = costheta2;
      false_costhetastar = costhetastar;
      false_phi = phi;
      false_phi1 = phi1;
      false_phi2 = phi2;
    }
  else
    {
      truth_met4vector = met4vector;
      truth_neutrino_deltaR = neutrino_deltaR;
      truth_costheta1 = costheta1;
      truth_costheta2 = costheta2;
      truth_costhetastar = costhetastar;
      truth_phi = phi;
      truth_phi1 = phi1;
      truth_phi2 = phi2;

      false_met4vector = other_met4vector;
      false_neutrino_deltaR = other_neutrino_deltaR;
      false_costheta1 = other_costheta1;
      false_costheta2 = other_costheta2;
      false_costhetastar = other_costhetastar;
      false_phi = other_phi;
      false_phi1 = other_phi1;
      false_phi2 = other_phi2;
    }

 

  Double_t min = 10000.00;
  Double_t jetDeltaR[4];
  Int_t tag = -1;

  jetDeltaR[0] = jet0.DeltaR(quark_gen);
  jetDeltaR[1] = jet1.DeltaR(quark_gen);
  jetDeltaR[2] = jet0.DeltaR(antiquark_gen);
  jetDeltaR[3] = jet1.DeltaR(antiquark_gen);

  for (int i = 0; i < 4; i++)
    {
      if (jetDeltaR[i] < min)
	{
	  min = jetDeltaR[i];
	  tag = i;
	}
    }

  if (tag == 0)      {jet0_deltaR = min; jet1_deltaR = jetDeltaR[3];}
  else if (tag == 1) {jet1_deltaR = min; jet0_deltaR = jetDeltaR[2];}
  else if (tag == 2) {jet0_deltaR = min; jet1_deltaR = jetDeltaR[1];}
  else if (tag == 3) {jet1_deltaR = min; jet0_deltaR = jetDeltaR[0];}
  else {std::cout << "not even close. The nearest jet is more than " << min << " away.\n"; exit(1);}


  min = 10000.00;   
  Double_t subjetDeltaR[4];
  tag = -1;

  subjetDeltaR[0] = subjet0.DeltaR(quark_gen);
  subjetDeltaR[1] = subjet1.DeltaR(quark_gen);
  subjetDeltaR[2] = subjet0.DeltaR(antiquark_gen);
  subjetDeltaR[3] = subjet1.DeltaR(antiquark_gen);

  for (int i = 0; i < 4; i++)
    {
      if (subjetDeltaR[i] < min)
	{
	  min = subjetDeltaR[i];
	  tag = i;
	}
    }

  if (tag == 0)      {subjet0_deltaR = min; subjet1_deltaR = subjetDeltaR[3];}
  else if (tag == 1) {subjet1_deltaR = min; subjet0_deltaR = subjetDeltaR[2];}
  else if (tag == 2) {subjet0_deltaR = min; subjet1_deltaR = subjetDeltaR[1];}
  else if (tag == 3) {subjet1_deltaR = min; subjet0_deltaR = subjetDeltaR[0];}
  else {std::cout << "not even close. The nearest subjet is more than " << min << " away.\n"; exit(1);}
  
  // MARKER ASSIGN
  //    mm    mmmm   mmmm  mmmmm    mmm  mm   m
  //    ##   #"   " #"   "   #    m"   " #"m  #
  //   #  #  "#mmm  "#mmm    #    #   mm # #m #
  //   #mm#      "#     "#   #    #    # #  # #
  //  #    # "mmm#" "mmm#" mm#mm   "mmm" #   ##


  /////////////////////////////////////LET'S TEST THE WWMASS
  //
  TLorentzVector quark_sum = quark_gen + antiquark_gen;
  TLorentzVector lepton_gen_sum = lepton_gen + neutrino_gen;
  TLorentzVector lepton_sum = lepton + met4vector;

  hadronic_gen_2body_mass = std::sqrt(quark_sum.Dot(quark_sum));
  leptonic_gen_2body_mass = std::sqrt(lepton_gen_sum.Dot(lepton_gen_sum));
  leptonic_2body_mass = std::sqrt(lepton_sum.Dot(lepton_sum));
  //
  ///////////////////////////////////////////////////

  /*//////
    TEMPLATE
    _px = .Px();
    _py = .Py();
    _pz = .Pz();
    _e = .E();
    _eta = .Eta();
    _pt = .Pt();
    _theta = .Theta();
    _phi = .Phi();
    _y = .Rapidity();
    _mass = .M();
    END TEMPLATE */

  electrn_reco.Px = electron.Px();
  electrn_reco.Py = electron.Py();
  electrn_reco.Pz = electron.Pz();
  electrn_reco.E = electron.E();
  electrn_reco.Eta = electron.Eta();
  electrn_reco.Pt = electron.Pt();
  electrn_reco.Theta = electron.Theta();
  electrn_reco.Phi = electron.Phi();
  electrn_reco.Y = electron.Rapidity();
  electrn_reco.M = electron.M();

  muon_reco.Px = muon.Px();
  muon_reco.Py = muon.Py();
  muon_reco.Pz = muon.Pz();
  muon_reco.E = muon.E();
  muon_reco.Eta = muon.Eta();
  muon_reco.Pt = muon.Pt();
  muon_reco.Theta = muon.Theta();
  muon_reco.Phi = muon.Phi();
  muon_reco.Y = muon.Rapidity();
  muon_reco.M = muon.M();

  lepton_reco.Px = lepton.Px();
  lepton_reco.Py = lepton.Py();
  lepton_reco.Pz = lepton.Pz();
  lepton_reco.E = lepton.E();
  lepton_reco.Eta = lepton.Eta();
  lepton_reco.Pt = lepton.Pt();
  lepton_reco.Theta = lepton.Theta();
  lepton_reco.Phi = lepton.Phi();
  lepton_reco.Y = lepton.Rapidity();
  lepton_reco.M = lepton.M();

  jt0.Px = jet0.Px();
  jt0.Py = jet0.Py();
  jt0.Pz = jet0.Pz();
  jt0.E = jet0.E();
  jt0.Eta = jet0.Eta();
  jt0.Pt = jet0.Pt();
  jt0.Theta = jet0.Theta();
  jt0.Phi = jet0.Phi();
  jt0.Y = jet0.Rapidity();
  jt0.M = jet0.M();

  jt1.Px = jet1.Px();
  jt1.Py = jet1.Py();
  jt1.Pz = jet1.Pz();
  jt1.E = jet1.E();
  jt1.Eta = jet1.Eta();
  jt1.Pt = jet1.Pt();
  jt1.Theta = jet1.Theta();
  jt1.Phi = jet1.Phi();
  jt1.Y = jet1.Rapidity();
  jt1.M = jet1.M();

  subjt0.Px = subjet0.Px();
  subjt0.Py = subjet0.Py();
  subjt0.Pz = subjet0.Pz();
  subjt0.E = subjet0.E();
  subjt0.Eta = subjet0.Eta();
  subjt0.Pt = subjet0.Pt();
  subjt0.Theta = subjet0.Theta();
  subjt0.Phi = subjet0.Phi();
  subjt0.Y = subjet0.Rapidity();
  subjt0.M = subjet0.M();

  subjt1.Px = subjet1.Px();
  subjt1.Py = subjet1.Py();
  subjt1.Pz = subjet1.Pz();
  subjt1.E = subjet1.E();
  subjt1.Eta = subjet1.Eta();
  subjt1.Pt = subjet1.Pt();
  subjt1.Theta = subjet1.Theta();
  subjt1.Phi = subjet1.Phi();
  subjt1.Y = subjet1.Rapidity();
  subjt1.M = subjet1.M();

  ntrno_reco.Px = met4vector.Px();
  ntrno_reco.Py = met4vector.Py();
  ntrno_reco.Pz = met4vector.Pz();
  ntrno_reco.E = met4vector.E();
  ntrno_reco.Eta = met4vector.Eta();
  ntrno_reco.Pt = met4vector.Pt();
  ntrno_reco.Theta = met4vector.Theta();
  ntrno_reco.Phi = met4vector.Phi();
  ntrno_reco.Y = met4vector.Rapidity();
  ntrno_reco.M = met4vector.M();

  other_ntrno_reco.Px = other_met4vector.Px();
  other_ntrno_reco.Py = other_met4vector.Py();
  other_ntrno_reco.Pz = other_met4vector.Pz();
  other_ntrno_reco.E = other_met4vector.E();
  other_ntrno_reco.Eta = other_met4vector.Eta();
  other_ntrno_reco.Pt = other_met4vector.Pt();
  other_ntrno_reco.Theta = other_met4vector.Theta();
  other_ntrno_reco.Phi = other_met4vector.Phi();
  other_ntrno_reco.Y = other_met4vector.Rapidity();
  other_ntrno_reco.M = other_met4vector.M();

  truth_ntrno_reco.Px = truth_met4vector.Px();
  truth_ntrno_reco.Py = truth_met4vector.Py();
  truth_ntrno_reco.Pz = truth_met4vector.Pz();
  truth_ntrno_reco.E = truth_met4vector.E();
  truth_ntrno_reco.Eta = truth_met4vector.Eta();
  truth_ntrno_reco.Pt = truth_met4vector.Pt();
  truth_ntrno_reco.Theta = truth_met4vector.Theta();
  truth_ntrno_reco.Phi = truth_met4vector.Phi();
  truth_ntrno_reco.Y = truth_met4vector.Rapidity();
  truth_ntrno_reco.M = truth_met4vector.M();

  false_ntrno_reco.Px = false_met4vector.Px();
  false_ntrno_reco.Py = false_met4vector.Py();
  false_ntrno_reco.Pz = false_met4vector.Pz();
  false_ntrno_reco.E = false_met4vector.E();
  false_ntrno_reco.Eta = false_met4vector.Eta();
  false_ntrno_reco.Pt = false_met4vector.Pt();
  false_ntrno_reco.Theta = false_met4vector.Theta();
  false_ntrno_reco.Phi = false_met4vector.Phi();
  false_ntrno_reco.Y = false_met4vector.Rapidity();
  false_ntrno_reco.M = false_met4vector.M();

  Whad_reco.Px = Whad.Px();
  Whad_reco.Py = Whad.Py();
  Whad_reco.Pz = Whad.Pz();
  Whad_reco.E = Whad.E();
  Whad_reco.Eta = Whad.Eta();
  Whad_reco.Pt = Whad.Pt();
  Whad_reco.Theta = Whad.Theta();
  Whad_reco.Phi = Whad.Phi();
  Whad_reco.Y = Whad.Rapidity();
  Whad_reco.M = Whad.M();

  Wlep_reco.Px = Wlep.Px();
  Wlep_reco.Py = Wlep.Py();
  Wlep_reco.Pz = Wlep.Pz();
  Wlep_reco.E = Wlep.E();
  Wlep_reco.Eta = Wlep.Eta();
  Wlep_reco.Pt = Wlep.Pt();
  Wlep_reco.Theta = Wlep.Theta();
  Wlep_reco.Phi = Wlep.Phi();
  Wlep_reco.Y = Wlep.Rapidity();
  Wlep_reco.M = Wlep.M();

  //assign gen particles
    
  Wpls_gen.Pt = Wplus_gen.Pt();
  Wmns_gen.Pt = Wminus_gen.Pt();
  Zntrl_gen.Pt = Zneutral_gen.Pt();
  qrk_gen.Pt = quark_gen.Pt();
  antiqrk_gen.Pt = antiquark_gen.Pt();
  lptn_gen.Pt = lepton_gen.Pt();
  ntrno_gen.Pt = neutrino_gen.Pt();

  Wpls_gen.Px = Wplus_gen.Px();
  Wmns_gen.Px = Wminus_gen.Px();
  Zntrl_gen.Px = Zneutral_gen.Px();
  qrk_gen.Px = quark_gen.Px();
  antiqrk_gen.Px = antiquark_gen.Px();
  lptn_gen.Px = lepton_gen.Px();
  ntrno_gen.Px = neutrino_gen.Px();

  Wpls_gen.Py = Wplus_gen.Py();
  Wmns_gen.Py = Wminus_gen.Py();
  Zntrl_gen.Py = Zneutral_gen.Py();
  qrk_gen.Py = quark_gen.Py();
  antiqrk_gen.Py = antiquark_gen.Py();
  lptn_gen.Py = lepton_gen.Py();
  ntrno_gen.Py = neutrino_gen.Py();

  Wpls_gen.Pz = Wplus_gen.Pz();
  Wmns_gen.Pz = Wminus_gen.Pz();
  Zntrl_gen.Pz = Zneutral_gen.Pz();
  qrk_gen.Pz = quark_gen.Pz();
  antiqrk_gen.Pz = antiquark_gen.Pz();
  lptn_gen.Pz = lepton_gen.Pz();
  ntrno_gen.Pz = neutrino_gen.Pz();

  Wpls_gen.E = Wplus_gen.E();
  Wmns_gen.E = Wminus_gen.E();
  Zntrl_gen.E = Zneutral_gen.E();
  qrk_gen.E = quark_gen.E();
  antiqrk_gen.E = antiquark_gen.E();
  lptn_gen.E = lepton_gen.E();
  ntrno_gen.E = neutrino_gen.E();

  Wpls_gen.Eta = Wplus_gen.Eta();
  Wmns_gen.Eta = Wminus_gen.Eta();
  Zntrl_gen.Eta = Zneutral_gen.Eta();
  qrk_gen.Eta = quark_gen.Eta();
  antiqrk_gen.Eta = antiquark_gen.Eta();
  lptn_gen.Eta = lepton_gen.Eta();
  ntrno_gen.Eta = neutrino_gen.Eta();

  Wpls_gen.Theta = Wplus_gen.Theta();
  Wmns_gen.Theta = Wminus_gen.Theta();
  Zntrl_gen.Theta = Zneutral_gen.Theta();
  qrk_gen.Theta = quark_gen.Theta();
  antiqrk_gen.Theta = antiquark_gen.Theta();
  lptn_gen.Theta = lepton_gen.Theta();
  ntrno_gen.Theta = neutrino_gen.Theta();

  Wpls_gen.Phi = Wplus_gen.Phi();
  Wmns_gen.Phi = Wminus_gen.Phi();
  Zntrl_gen.Phi = Zneutral_gen.Phi();
  qrk_gen.Phi = quark_gen.Phi();
  antiqrk_gen.Phi = antiquark_gen.Phi();
  lptn_gen.Phi = lepton_gen.Phi();
  ntrno_gen.Phi = neutrino_gen.Phi();

  Wpls_gen.Y = Wplus_gen.Rapidity();
  Wmns_gen.Y = Wminus_gen.Rapidity();
  Zntrl_gen.Y = Zneutral_gen.Rapidity();
  qrk_gen.Y = quark_gen.Rapidity();
  antiqrk_gen.Y = antiquark_gen.Rapidity();
  lptn_gen.Y = lepton_gen.Rapidity();
  ntrno_gen.Y = neutrino_gen.Rapidity();

  Wpls_gen.M = Wplus_gen.M();
  Wmns_gen.M = Wminus_gen.M();
  Zntrl_gen.M = Zneutral_gen.M();
  qrk_gen.M = quark_gen.M();
  antiqrk_gen.M = antiquark_gen.M();
  lptn_gen.M = lepton_gen.M();
  ntrno_gen.M = neutrino_gen.M();

  Tau21 = tau21; Tau32 = tau32; PTau21 = puppi_tau21; PTau32 = puppi_tau32;
  
#ifdef MULTI_WEIGHTS
  //MARKER ANOMALOUS WEIGHTS
  edm::Handle<Double_t> weight1st;
  iEvent.getByToken(mgreweight1_, weight1st);
  anom_weight1 = *weight1st;

  edm::Handle<Double_t> weight2nd;
  iEvent.getByToken(mgreweight2_, weight2nd);
  anom_weight2 = *weight2nd;

  edm::Handle<Double_t> weight3rd;
  iEvent.getByToken(mgreweight3_, weight3rd);
  anom_weight3 = *weight3rd;

  edm::Handle<Double_t> weight4th;
  iEvent.getByToken(mgreweight4_, weight4th);
  anom_weight4 = *weight4th;

  edm::Handle<Double_t> weight5th;
  iEvent.getByToken(mgreweight5_, weight5th);
  anom_weight5 = *weight5th;

  edm::Handle<Double_t> weight6th;
  iEvent.getByToken(mgreweight6_, weight6th);
  anom_weight6 = *weight6th;

  edm::Handle<Double_t> weight7th;
  iEvent.getByToken(mgreweight7_, weight7th);
  anom_weight7 = *weight7th;

  /*edm::Handle<Double_t> weight8th;
  iEvent.getByToken(mgreweight8_, weight8th);
  anom_weight8 = *weight8th;

  edm::Handle<Double_t> weight9th;
  iEvent.getByToken(mgreweight9_, weight9th);
  anom_weight9 = *weight9th;

  edm::Handle<Double_t> weight10th;
  iEvent.getByToken(mgreweight10_, weight10th);
  anom_weight10 = *weight10th;

  edm::Handle<Double_t> weight11th;
  iEvent.getByToken(mgreweight11_, weight11th);
  anom_weight11 = *weight11th;

  edm::Handle<Double_t> weight12th;
  iEvent.getByToken(mgreweight12_, weight12th);
  anom_weight12 = *weight12th;

  edm::Handle<Double_t> weight13th;
  iEvent.getByToken(mgreweight13_, weight13th);
  anom_weight13 = *weight13th;

  edm::Handle<Double_t> weight14th;
  iEvent.getByToken(mgreweight14_, weight14th);
  anom_weight14 = *weight14th;

  edm::Handle<Double_t> weight15th;
  iEvent.getByToken(mgreweight15_, weight15th);
  anom_weight15 = *weight15th;

  edm::Handle<Double_t> weight16th;
  iEvent.getByToken(mgreweight16_, weight16th);
  anom_weight16 = *weight16th;

  edm::Handle<Double_t> weight17th;
  iEvent.getByToken(mgreweight17_, weight17th);
  anom_weight17 = *weight17th;

  edm::Handle<Double_t> weight18th;
  iEvent.getByToken(mgreweight18_, weight18th);
  anom_weight18 = *weight18th;

  edm::Handle<Double_t> weight19th;
  iEvent.getByToken(mgreweight19_, weight19th);
  anom_weight19 = *weight19th;

  edm::Handle<Double_t> weight20th;
  iEvent.getByToken(mgreweight20_, weight20th);
  anom_weight20 = *weight20th;
/ *
  edm::Handle<Double_t> weight21st;
  iEvent.getByToken(mgreweight21_, weight21st);
  anom_weight21 = *weight21st;

  edm::Handle<Double_t> weight22nd;
  iEvent.getByToken(mgreweight22_, weight22nd);
  anom_weight22 = *weight22nd;

  edm::Handle<Double_t> weight23rd;
  iEvent.getByToken(mgreweight23_, weight23rd);
  anom_weight23 = *weight23rd;

  edm::Handle<Double_t> weight24th;
  iEvent.getByToken(mgreweight24_, weight24th);
  anom_weight24 = *weight24th;

  edm::Handle<Double_t> weight25th;
  iEvent.getByToken(mgreweight25_, weight25th);
  anom_weight25 = *weight25th;

  edm::Handle<Double_t> weight26th;
  iEvent.getByToken(mgreweight26_, weight26th);
  anom_weight26 = *weight26th;

  edm::Handle<Double_t> weight27th;
  iEvent.getByToken(mgreweight27_, weight27th);
  anom_weight27 = *weight27th;

  edm::Handle<Double_t> weight28th;
  iEvent.getByToken(mgreweight28_, weight28th);
  anom_weight28 = *weight28th;

  edm::Handle<Double_t> weight29th;
  iEvent.getByToken(mgreweight29_, weight29th);
  anom_weight29 = *weight29th;

  edm::Handle<Double_t> weight30th;
  iEvent.getByToken(mgreweight30_, weight30th);
  anom_weight30 = *weight30th;

  edm::Handle<Double_t> weight31st;
  iEvent.getByToken(mgreweight31_, weight31st);
  anom_weight31 = *weight31st;

  edm::Handle<Double_t> weight32nd;
  iEvent.getByToken(mgreweight32_, weight32nd);
  anom_weight32 = *weight32nd;

  edm::Handle<Double_t> weight33rd;
  iEvent.getByToken(mgreweight33_, weight33rd);
  anom_weight33 = *weight33rd;

  edm::Handle<Double_t> weight34th;
  iEvent.getByToken(mgreweight34_, weight34th);
  anom_weight34 = *weight34th;

  edm::Handle<Double_t> weight35th;
  iEvent.getByToken(mgreweight35_, weight35th);
  anom_weight35 = *weight35th;

  edm::Handle<Double_t> weight36th;
  iEvent.getByToken(mgreweight36_, weight36th);
  anom_weight36 = *weight36th;

  edm::Handle<Double_t> weight37th;
  iEvent.getByToken(mgreweight37_, weight37th);
  anom_weight37 = *weight37th;

  edm::Handle<Double_t> weight38th;
  iEvent.getByToken(mgreweight38_, weight38th);
  anom_weight38 = *weight38th;

  edm::Handle<Double_t> weight39th;
  iEvent.getByToken(mgreweight39_, weight39th);
  anom_weight39 = *weight39th;

  edm::Handle<Double_t> weight40th;
  iEvent.getByToken(mgreweight40_, weight40th);
  anom_weight40 = *weight40th;

  edm::Handle<Double_t> weight41st;
  iEvent.getByToken(mgreweight41_, weight41st);
  anom_weight41 = *weight41st;

  edm::Handle<Double_t> weight42nd;
  iEvent.getByToken(mgreweight42_, weight42nd);
  anom_weight42 = *weight42nd;

  edm::Handle<Double_t> weight43rd;
  iEvent.getByToken(mgreweight43_, weight43rd);
  anom_weight43 = *weight43rd;

  edm::Handle<Double_t> weight44th;
  iEvent.getByToken(mgreweight44_, weight44th);
  anom_weight44 = *weight44th;

  edm::Handle<Double_t> weight45th;
  iEvent.getByToken(mgreweight45_, weight45th);
  anom_weight45 = *weight45th;

  edm::Handle<Double_t> weight46th;
  iEvent.getByToken(mgreweight46_, weight46th);
  anom_weight46 = *weight46th;

  edm::Handle<Double_t> weight47th;
  iEvent.getByToken(mgreweight47_, weight47th);
  anom_weight47 = *weight47th;

  edm::Handle<Double_t> weight48th;
  iEvent.getByToken(mgreweight48_, weight48th);
  anom_weight48 = *weight48th;

  edm::Handle<Double_t> weight49th;
  iEvent.getByToken(mgreweight49_, weight49th);
  anom_weight49 = *weight49th;

  edm::Handle<Double_t> weight50th;
  iEvent.getByToken(mgreweight50_, weight50th);
  anom_weight50 = *weight50th;
*/
#endif
  //MARKER FILL TREE
  //
  //  mmmmmm   "    ""#    ""#          mmmmmmm                     
  //  #      mmm      #      #             #     m mm   mmm    mmm  
  //  #mmmmm   #      #      #             #     #"  " #"  #  #"  # 
  //  #        #      #      #             #     #     #""""  #"""" 
  //  #      mm#mm    "mm    "mm           #     #     "#mm"  "#mm" 
  //
  //These are histograms primarily for comparing this full simulation to a toy monte carlo that was developed to run simulations on neutrino pz schemes

  if (bad_event) return; //{std::cout << "bad event reached fill tree\n"; return;}
  if (bad_event) std::cout << "you shouldn't be seeing this.";
  
  if (!bad_event) MyTree->Fill();
  if ((abs(reco_lep_generation) == 1 && loose_electron == 0 && extra_electron == 0 && electron_out_of_bounds == 0) || (abs(reco_lep_generation) == 2 && loose_muon == 0 && extra_muon == 0 && muon_out_of_bounds == 0))
    {
      if (Wlep_set == 0)
	{
	  std::cout << "You fucked up, bitch!\n";
	  std::cout << "mu_event is " << mu_event << std::endl;
	  std::cout << "el_event is " << el_event << std::endl;
	  std::cout << "reco_lep_generation is " << reco_lep_generation << std::endl;
	  exit(23);
	}
      
      wp_mass->Fill(Wplus_gen.M());
      wm_mass->Fill(Wminus_gen.M());
      //pt_smear->Fill(neutrino_gen.Pt() - met4vector.Pt());

      Double_t anti0 = antiquark_gen.DeltaR(subjet0);//0
      Double_t anti1 = antiquark_gen.DeltaR(subjet1);//1
      Double_t quark0 = quark_gen.DeltaR(subjet0);//2
      Double_t quark1 = quark_gen.DeltaR(subjet1);//3

      Double_t minimum_value = anti0;
      Int_t min_number = 0;
      
      if (anti1 < minimum_value) {minimum_value = anti1; min_number = 1;}
      if (quark0 < minimum_value) {minimum_value = quark0; min_number = 2;}
      if (quark1 < minimum_value) {minimum_value = quark1; min_number = 3;}

      if (good_jet == 1 && (min_number == 0 || min_number == 3))
	{
	  parton0_pt_smear->Fill(antiquark_gen.Pt() - subjet0.Pt());
	  parton1_pt_smear->Fill(quark_gen.Pt() - subjet1.Pt());

	  parton0_phi_smear->Fill(antiquark_gen.Phi() - subjet0.Phi());
	  parton1_phi_smear->Fill(quark_gen.Phi() - subjet1.Phi());

	  parton0_eta_smear->Fill(antiquark_gen.Eta() - subjet0.Eta());
	  parton1_eta_smear->Fill(quark_gen.Eta() - subjet1.Eta());
	}
      
      if (good_jet == 1 && (min_number == 1 || min_number == 2))
	{
	  parton0_pt_smear->Fill(quark_gen.Pt() - subjet0.Pt());
	  parton1_pt_smear->Fill(antiquark_gen.Pt() - subjet1.Pt());

	  parton0_phi_smear->Fill(quark_gen.Phi() - subjet0.Phi());
	  parton1_phi_smear->Fill(antiquark_gen.Phi() - subjet1.Phi());

	  parton0_eta_smear->Fill(quark_gen.Eta() - subjet0.Eta());
	  parton1_eta_smear->Fill(antiquark_gen.Eta() - subjet1.Eta());
	}
      

      TLorentzVector WW_vector = Wplus_gen + Wminus_gen;
      ww_mass->Fill(WW_vector.M());
      ww_pt->Fill(Wplus_gen.Pt() + Wminus_gen.Pt());
      ww_rap->Fill(WW_vector.Rapidity());

      if (imaginary_neutrino == 0 && good_jet == 1)
	{
	  comp_hist0->Fill(neutrino_gen.Pz() - fill_pz0);
	  comp_hist1->Fill(neutrino_gen.Pz() - fill_pz1);
	  comp_hist2->Fill(neutrino_gen.Pz() - fill_pz2);
	  comp_hist3->Fill(neutrino_gen.Pz() - fill_pz3);
	  comp_hist4->Fill(neutrino_gen.Pz() - fill_pz4);
	  comp_hist5->Fill(neutrino_gen.Pz() - fill_pz5);
	  comp_hist6->Fill(neutrino_gen.Pz() - fill_pz6);
	  comp_hist7->Fill(neutrino_gen.Pz() - fill_pz7);
	}

      wp_pt->Fill(Wplus_gen.Pt());
      wm_pt->Fill(Wminus_gen.Pt());
      if (good_jet == 1) whad_pt->Fill(Whad.Pt());
      wlep_pt->Fill(Wlep.Pt());
      neu_pt->Fill(neutrino_gen.Pt());
      met_pt->Fill(met4vector.Pt());

      wp_rap->Fill(Wplus_gen.Rapidity());
      wm_rap->Fill(Wminus_gen.Rapidity());
      if (good_jet == 1) whad_rap->Fill(Whad.Rapidity());
      wlep_rap->Fill(Wlep.Rapidity());
      neu_rap->Fill(neutrino_gen.Rapidity());
      met_rap->Fill(met4vector.Rapidity());

      wplus_theta->Fill(Wplus_gen.Theta()); wminus_theta->Fill(Wminus_gen.Theta());
      lepton_theta->Fill(lepton_gen.Theta()); neutrino_theta->Fill(neutrino_gen.Theta());
      quark_theta->Fill(quark_gen.Theta()); antiquark_theta->Fill(antiquark_gen.Theta());
      
      wplus_phi->Fill(Wplus_gen.Phi()); wminus_phi->Fill(Wminus_gen.Phi());
      lepton_phi->Fill(lepton_gen.Phi()); neutrino_phi->Fill(neutrino_gen.Phi());
      quark_phi->Fill(quark_gen.Phi()); antiquark_phi->Fill(antiquark_gen.Phi());

      if (reco_lep_generation > 0)
	{
	  if (good_jet == 1) whad_pt_smear->Fill(Wplus_gen.Pt() - Whad.Pt());
	  if (good_jet == 1) whad_phi_smear->Fill(Wplus_gen.Phi() - Whad.Phi());
	  if (good_jet == 1) whad_eta_smear->Fill(Wplus_gen.Eta() - Whad.Eta());
	  wlep_pt_smear->Fill(Wminus_gen.Pt() - Wlep.Pt());
	  wlep_phi_smear->Fill(Wminus_gen.Phi() - Wlep.Phi());
	  wlep_eta_smear->Fill(Wminus_gen.Eta() - Wlep.Eta());
	}
      if (reco_lep_generation < 0)
	{
	  if (good_jet == 1) whad_pt_smear->Fill(Wminus_gen.Pt() - Whad.Pt());
	  if (good_jet == 1) whad_phi_smear->Fill(Wminus_gen.Phi() - Whad.Phi());
	  if (good_jet == 1) whad_eta_smear->Fill(Wminus_gen.Eta() - Whad.Eta());
	  wlep_pt_smear->Fill(Wplus_gen.Pt() - Wlep.Pt());
	  wlep_phi_smear->Fill(Wplus_gen.Phi() - Wlep.Phi());
	  wlep_eta_smear->Fill(Wplus_gen.Eta() - Wlep.Eta());
	}
      
      met_pt_smear->Fill(neutrino_gen.Pt() - met4vector.Pt());
      met_phi_smear->Fill(neutrino_gen.Phi() - met4vector.Phi());
      met_eta_smear->Fill(neutrino_gen.Eta() - met4vector.Eta());
            
      if (abs(reco_lep_generation) == 1)
	{
	  if (Wlep_set != 1) {std::cout << "You fucked up, bitch!\n"; exit(24);}
	  //Double_t diff = lepton.Pt() - electron.Pt();
	  //if (lepton.Pt() != electron.Pt()) {std::cout << "error with assignment: pt difference is: " << diff << std::endl; exit(22);}
	  //if (electron_value_set == 0) {std::cout << "reco_lep_generation = " << reco_lep_generation << std::endl << "el_mark = " << el_mark << " how is the electron not set?\n"; exit(22);}
	  
	  elr_pt->Fill(electron.Pt());
	  elg_pt->Fill(lepton_gen.Pt());

	  elr_rap->Fill(electron.Rapidity());
	  elg_rap->Fill(lepton_gen.Rapidity());

	  el_pt_smear->Fill(lepton_gen.Pt() - electron.Pt());
	  el_phi_smear->Fill(lepton_gen.Phi() - electron.Phi());
	  el_eta_smear->Fill(lepton_gen.Eta() - electron.Eta());
	}
      
      if (abs(reco_lep_generation) == 2)
	{
	   if (Wlep_set != 2) {std::cout << "You fucked up, bitch!\n"; exit(25);}
	  //Double_t diff = lepton.Pt() - muon.Pt();
	  //if (lepton.Pt() != muon.Pt()) {std::cout << "error with assignment: pt difference is: " << diff << std::endl; exit(22);}
	  //if (muon_value_set == 0) {std::cout << "reco_lep_generation = " << reco_lep_generation << std::endl << "mu_mark = " << mu_mark << " how is the muon not set?\n"; exit(22);}
	  
	  mur_pt->Fill(muon.Pt());    
	  mug_pt->Fill(lepton_gen.Pt());

	  mur_rap->Fill(muon.Rapidity());    
	  mug_rap->Fill(lepton_gen.Rapidity());

	  mu_pt_smear->Fill(lepton_gen.Pt() - muon.Pt());
	  mu_phi_smear->Fill(lepton_gen.Phi() - muon.Phi());
	  mu_eta_smear->Fill(lepton_gen.Eta() - muon.Eta());
	}
     
      //good_events++;
      //outfile <<  "got good event number " << good_events << std::endl;
    }
  
  if (talking) outfile << "this was event number " << good_events << std::endl;
 
  
//#endif
  /////////////////////////////////////////////////////////////////

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Double_t& costheta1, Double_t& costheta2, Double_t& phi, Double_t& costhetastar, Double_t& phi1, Double_t& phi2)
{


  TLorentzVector thep4H = thep4M11 + thep4M12 + thep4M21 + thep4M22;
  TLorentzVector thep4Z1 = thep4M11 + thep4M12;
  TLorentzVector thep4Z2 = thep4M21 + thep4M22;

  Double_t norm;

  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(thep4Z1);
  TLorentzVector thep4Z2inXFrame(thep4Z2);      
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());

  // calculate phi1, phi2, costhetastar
  ///phi1 = theZ1X_p3.Phi();
  ///phi2 = theZ2X_p3.Phi();

  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  /////////////////////////////////////////////// 
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();

  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  //find the decay axis
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
  TVector3 unitx_1(-p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z());
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  //boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  //create z and y axes
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
  TVector3 p4M21Z1_p3(p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z());
  TVector3 p4M22Z1_p3(p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z());
  TVector3 unitz_1 = p4M21Z1_p3.Cross(p4M22Z1_p3);
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);

  //caculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11(p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z());
  TVector3 unitM11 = p3M11.Unit();
  Double_t x_m11 = unitM11.Dot(unitx_1); Double_t y_m11 = unitM11.Dot(unity_1); Double_t z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();
  //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();

  //set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2(-p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z());
  norm = 1/(unitx_2.Mag());
  unitx_2 *= norm;
  //boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3(p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z());
  TVector3 p4M12Z2_p3(p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z());
  TVector3 unitz_2 = p4M11Z2_p3.Cross(p4M12Z2_p3);
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  //calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21(p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z());
  TVector3 unitM21 = p3M21.Unit();
  Double_t x_m21 = unitM21.Dot(unitx_2); Double_t y_m21 = unitM21.Dot(unity_2); Double_t z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();

  // calculate phi
  //calculating phi_n
  TLorentzVector n_p4Z1inXFrame(p4Z1);
  TLorentzVector n_p4M11inXFrame(p4M11);
  n_p4Z1inXFrame.Boost(boostX);
  n_p4M11inXFrame.Boost(boostX);        
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
  TVector3 n_unitz_1(n_p4Z1inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross(n_unitz_1);
  TVector3 n_unity_1 = n_unitz_1.Cross(n_p4M11inXFrame_unit);
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1);

  TLorentzVector n_p4M21inXFrame(p4M21);
  n_p4M21inXFrame.Boost(boostX);
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime(n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1));

  ///////-----------------new way of calculating phi-----------------///////
  //Double_t phi_n =  n_p4M21inXFrame_unitprime.Phi();
  //
    //std::cout << "---------------------------" << std::endl;
    //std::cout << "phi: " << phi << std::endl;
    //std::cout << "phi_n: " << phi_n << std::endl;
    //std::cout << "phi + phi_n: " << (phi+phi_n) << std::endl;
  //
  /// and then calculate phi1
  TVector3 n_p4PartoninXFrame_unit(0.0, 0.0, 1.0);
  TVector3 n_p4PartoninXFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1));
  // negative sign is for arrow convention in paper
  phi1 = (n_p4PartoninXFrame_unitprime.Phi());

  // and the calculate phi2
  TLorentzVector n_p4Z2inXFrame(p4Z2);
  n_p4Z2inXFrame.Boost(boostX);
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  ///////TLorentzVector n_p4M21inXFrame(p4M21);
  //////n_p4M21inXFrame.Boost(boostX);        
  ////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
  TVector3 n_unitz_2(n_p4Z2inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross(n_unitz_2);
  TVector3 n_unity_2 = n_unitz_2.Cross(n_p4M21inXFrame_unit);
  TVector3 n_unitx_2 = n_unity_2.Cross(n_unitz_2);
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2));
  phi2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());

  //
    //Double_t phi12_0 = phi1 + phi2;
    //if (phi12_0 > TMath::Pi()) phi12 = phi12_0 - 2*TMath::Pi();
    //else if (phi12_0 < (-1.)*TMath::Pi()) phi12 = phi12_0 + 2*TMath::Pi();
    //else phi12 = phi12_0;
  //

}


void intermediate_steps(TLorentzVector lepton1, TLorentzVector lepton2, TLorentzVector parton1, TLorentzVector parton2, Double_t& leptons_in_lep_px, Double_t& leptons_in_lep_py, Double_t& leptons_in_lep_pz, Double_t& partons_in_lep_px, Double_t& partons_in_lep_py, Double_t& partons_in_lep_pz, Double_t& parton1_in_lep_px, Double_t& parton2_in_lep_px, Double_t& parton1_in_lep_py, Double_t& parton2_in_lep_py, Double_t& parton1_in_lep_pz, Double_t& parton2_in_lep_pz, Double_t& lepton1_in_lep_px, Double_t& lepton1_in_lep_py, Double_t& lepton1_in_lep_pz, Double_t& lepton1_dotted_x, Double_t& lepton1_dotted_y, Double_t& lepton1_dotted_z, Double_t& leptons_in_had_px, Double_t& leptons_in_had_py, Double_t& leptons_in_had_pz, Double_t& lepton1_in_had_px, Double_t& lepton1_in_had_py, Double_t& lepton1_in_had_pz, Double_t& lepton2_in_had_px, Double_t& lepton2_in_had_py, Double_t& lepton2_in_had_pz, Double_t& parton1_in_had_px, Double_t& parton1_in_had_py, Double_t& parton1_in_had_pz, Double_t& parton1_dotted_x, Double_t& parton1_dotted_y, Double_t& parton1_dotted_z, Double_t& complicated1_px, Double_t& complicated1_py, Double_t& complicated1_pz, Double_t& complicated2_px, Double_t& complicated2_py, Double_t& complicated2_pz, Double_t& lepton_sumWWframe_X, Double_t& lepton_sumWWframe_Y, Double_t& lepton_sumWWframe_Z, Double_t& lepton1WWframe_X, Double_t& lepton1WWframe_Y, Double_t& lepton1WWframe_Z, Double_t& parton_sumWWframe_X, Double_t& parton_sumWWframe_Y, Double_t& parton_sumWWframe_Z, Double_t& parton1WWframe_X, Double_t& parton1WWframe_Y, Double_t& parton1WWframe_Z, Double_t& costhetastar, Double_t& costheta1, Double_t& phi, Double_t& costheta2, Double_t& phi1, Double_t& phi2, Double_t& boostWWframe_X, Double_t& boostWWframe_Y, Double_t& boostWWframe_Z, Double_t& boostWlep_X, Double_t& boostWlep_Y, Double_t& boostWlep_Z, Double_t& boostWhad_X, Double_t& boostWhad_Y, Double_t& boostWhad_Z)
{


  TLorentzVector fermion_sum = lepton1 + lepton2 + parton1 + parton2; //4-vector sum of all 4 particles XXX 
  TLorentzVector lepton_sum = lepton1 + lepton2; //4-vector sum of lepton and neutrino XXX
  TLorentzVector parton_sum = parton1 + parton2; //4-vector sum of partons XXX

  Double_t norm;

  TVector3 boostWWframe = -(fermion_sum.BoostVector()); //boost to WW rest frame XXX
  boostWWframe_X = boostWWframe.X(); //newer
  boostWWframe_Y = boostWWframe.Y();
  boostWWframe_Z = boostWWframe.Z();
  
  TLorentzVector lepton_sum_WWframe(lepton_sum); //4-vector sum of leptons boosted to WW rest frame (below) XXX  
  TLorentzVector parton_sum_WWframe(parton_sum); //4-vector sum of partons boosted to WW rest frame (below) XXX 
  lepton_sum_WWframe.Boost(boostWWframe);
  parton_sum_WWframe.Boost(boostWWframe);
  

  leptons_in_lep_px = lepton_sum_WWframe.X(); //new
  leptons_in_lep_py = lepton_sum_WWframe.Y();
  leptons_in_lep_pz = lepton_sum_WWframe.Z();
  
  if (leptons_in_lep_px != leptons_in_lep_px) leptons_in_lep_px = 800.0;
  if (leptons_in_lep_py != leptons_in_lep_py) leptons_in_lep_py = 700.0;
  if (leptons_in_lep_pz != leptons_in_lep_pz) leptons_in_lep_pz = 6000.0;
  
 TVector3 lepton_sum_WWframe3vec = TVector3(leptons_in_lep_px, leptons_in_lep_py, leptons_in_lep_pz); //3-vector sum of leptons in WW frame XXX
  
  //////////////////////////////////////////////////////////////////
  
  // calculate phi1, phi2, costhetastar

  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////
  
  //TLorentzVector fermion_sum, lepton_sum, lepton1, lepton2, parton_sum, parton1, parton2;
  
  costhetastar = lepton_sum_WWframe3vec.CosTheta(); //ANGLE THETA* 

  // now helicity angles................................
  // ...................................................
  TVector3 boostWlep = -(lepton_sum.BoostVector()); //boost to leptonic W rest frame XXX
  
  boostWlep_X = boostWlep.X(); //newer
  boostWlep_Y = boostWlep.Y();
  boostWlep_Z = boostWlep.Z();
  
  TLorentzVector parton_sumWlep_frame(parton_sum); //4-vector sum of partons boosted to leptonic W rest frame (below) XXX
  parton_sumWlep_frame.Boost(boostWlep);
  //find the decay axis
 

  partons_in_lep_px = parton_sumWlep_frame.X(); //new
  partons_in_lep_py = parton_sumWlep_frame.Y();
  partons_in_lep_pz = parton_sumWlep_frame.Z();
  
  if (partons_in_lep_px != partons_in_lep_px) partons_in_lep_px = 40000.0;
  if (partons_in_lep_py != partons_in_lep_py) partons_in_lep_py = 50000.0;
  if (partons_in_lep_pz != partons_in_lep_pz) partons_in_lep_pz = -20000.0;
  
  TVector3 unitx_1(-partons_in_lep_px, -partons_in_lep_py, -partons_in_lep_pz); //unit 3-vector (below) of partons boosted to leptonic W rest frame XXX
  norm = 1.0 / (unitx_1.Mag());
  unitx_1 *= norm;
  
  //boost daughters of z2
  TLorentzVector parton1Wlep_frame(parton1); //first parton boosted (below) to leptonic W rest frame XXX
  TLorentzVector parton2Wlep_frame(parton2); //second parton boosted (below) to leptonic W rest frame XXX
  parton1Wlep_frame.Boost(boostWlep);
  parton2Wlep_frame.Boost(boostWlep);
  
  //create z and y axes

  parton1_in_lep_px = parton1Wlep_frame.X(); //new
  parton2_in_lep_px = parton2Wlep_frame.X();
  
  parton1_in_lep_py = parton1Wlep_frame.Y();
  parton2_in_lep_py = parton2Wlep_frame.Y();
  
  parton1_in_lep_pz = parton1Wlep_frame.Z();
  parton2_in_lep_pz = parton2Wlep_frame.Z();
  
  if (parton1_in_lep_px != parton1_in_lep_px) parton1_in_lep_px = 30000.0;
  if (parton1_in_lep_py != parton1_in_lep_py) parton1_in_lep_py = 40000.0;
  if (parton1_in_lep_pz != parton1_in_lep_pz) parton1_in_lep_pz = -13000.0;

  if (parton2_in_lep_px != parton2_in_lep_px) parton2_in_lep_px = 30000.0;
  if (parton2_in_lep_py != parton2_in_lep_py) parton2_in_lep_py = 30000.0;
  if (parton2_in_lep_pz != parton2_in_lep_pz) parton2_in_lep_pz = -8000.0;
  
  TVector3 parton1Wlep_frame3vec(parton1_in_lep_px, parton1_in_lep_py, parton1_in_lep_pz); //3-vector of first parton in leptonic frame XXX
  TVector3 parton2Wlep_frame3vec(parton2_in_lep_px, parton2_in_lep_py, parton2_in_lep_pz); //3-vector of second parton in leptonic frame XXX

  
  TVector3 unitz_1 = parton1Wlep_frame3vec.Cross(parton2Wlep_frame3vec); //unit 3-vector (below) of cross product of partons in leptonic frame XXX
  norm = 1.0 / (unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1); //cross product of sum of partons with cross of partons, all in Wlep frame: (parton1 + parton2) X (parton1 X parton2) XXX

  ///unit1: x = partons in lep ; z = cross product of partons and leptons ; y = z cross x

  //caculate theta1
  TLorentzVector lepton1Wlep_frame(lepton1); //first lepton boosted (below) to leptonic W rest frame XXX
  lepton1Wlep_frame.Boost(boostWlep);

  lepton1_in_lep_px = lepton1Wlep_frame.X(); // new
  lepton1_in_lep_py = lepton1Wlep_frame.Y();
  lepton1_in_lep_pz = lepton1Wlep_frame.Z();
  
  if (lepton1_in_lep_px != lepton1_in_lep_px) lepton1_in_lep_px = -500.0;
  if (lepton1_in_lep_py != lepton1_in_lep_py) lepton1_in_lep_py = -500.0;
  if (lepton1_in_lep_pz != lepton1_in_lep_pz) lepton1_in_lep_pz = -300.0;
  
  TVector3 lepton1Wlep_frame3vec(lepton1_in_lep_px, lepton1_in_lep_py, lepton1_in_lep_pz); //3-vector of 1st lepton in leptonic frame XXX
  TVector3 lepton1Wlep_frame_unit3vec = lepton1Wlep_frame3vec.Unit(); //unit 3-vector of 1st lepton in leptonic frame XXX

  Double_t x_m11 = lepton1Wlep_frame_unit3vec.Dot(unitx_1);
  Double_t y_m11 = lepton1Wlep_frame_unit3vec.Dot(unity_1);
  Double_t z_m11 = lepton1Wlep_frame_unit3vec.Dot(unitz_1); //dot products of unit 3-vector of 1st lepton in lepton frame and parton units in lepton frame XXX

  lepton1_dotted_x = x_m11; // new
  lepton1_dotted_y = y_m11;
  lepton1_dotted_z = z_m11;
  
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11); //3-vector made of dot products of 1st lepton XXX


  ///////////////////////////////////////////////
  
  costheta1 = M11_Z1frame.CosTheta(); //ANGLE THETA1
 
  phi = M11_Z1frame.Phi(); ///ANGLE PHI

  //////////////////////////////////////////////////////////////////////// now on to last part

  //set axes for other system
  TVector3 boostWhad = -(parton_sum.BoostVector()); //boost to hadronic W rest frame XXX
  boostWhad_X = boostWhad.X(); //newer
  boostWhad_Y = boostWhad.Y();
  boostWhad_Z = boostWhad.Z();

  TLorentzVector lepton_sumWhad_frame(lepton_sum); //sum of leptonic vectors boosted to hadronic frame (below) XXX
  lepton_sumWhad_frame.Boost(boostWhad);

  leptons_in_had_px = lepton_sumWhad_frame.X(); //new
  leptons_in_had_py = lepton_sumWhad_frame.Y();
  leptons_in_had_pz = lepton_sumWhad_frame.Z();
  
  if (leptons_in_had_px != leptons_in_had_px) leptons_in_had_px = -15000.0;
  if (leptons_in_had_py != leptons_in_had_py) leptons_in_had_py = -15000.0;
  if (leptons_in_had_pz != leptons_in_had_pz) leptons_in_had_pz = -20000.0;
  
  TVector3 unitx_2(-leptons_in_had_px, -leptons_in_had_py, -leptons_in_had_pz); //unit 3-vector (below) of leptonic/hadronic boost XXX
  
  norm = 1.0 / (unitx_2.Mag());
  unitx_2 *= norm;
  
  //boost daughters of z2 //////////////////
  TLorentzVector lepton1Whad_frame(lepton1); //first lepton boosted to hadronic frame XXX
  TLorentzVector lepton2Whad_frame(lepton2); //second lepton boosted to hadronic frame XXX
  lepton1Whad_frame.Boost(boostWhad);
  lepton2Whad_frame.Boost(boostWhad);

  lepton1_in_had_px = lepton1Whad_frame.X(); // new
  lepton2_in_had_px = lepton2Whad_frame.X();
  
  lepton1_in_had_py = lepton1Whad_frame.Y();
  lepton2_in_had_py = lepton2Whad_frame.Y();
  
  lepton1_in_had_pz = lepton1Whad_frame.Z();
  lepton2_in_had_pz = lepton2Whad_frame.Z();

  ////////////////
  if (lepton1_in_had_px != lepton1_in_had_px) lepton1_in_had_px = -8000.0;
  if (lepton1_in_had_py != lepton1_in_had_py) lepton1_in_had_py = -10000.0;
  if (lepton1_in_had_pz != lepton1_in_had_pz) lepton1_in_had_pz = -14000.0;

  if (lepton2_in_had_px != lepton2_in_had_px) lepton2_in_had_px = -8000.0;
  if (lepton2_in_had_py != lepton2_in_had_py) lepton2_in_had_py = +7000.0;
  if (lepton2_in_had_pz != lepton2_in_had_pz) lepton2_in_had_pz = -14000.0;
  
  TVector3 lepton1Whad_frame3vec(lepton1_in_had_px, lepton1_in_had_py, lepton1_in_had_pz); //3-vector of 1st lepton in hadronic frame XXX
  TVector3 lepton2Whad_frame3vec(lepton2_in_had_px, lepton2_in_had_py, lepton2_in_had_pz); //3-vector of 2nd lepton in hadronic frame XXX
  
  TVector3 unitz_2 = lepton1Whad_frame3vec.Cross(lepton2Whad_frame3vec); //unit of cross of leptons in hadronic frame XXX
  norm = 1.0 / (unitz_2.Mag());
  unitz_2 *= norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2); //cross of other units XXX

  ///unit2: x = lepton sum in had; z = lepton1 cross lepton 2, all in had frame; y = z cross z
  
  //calcuate theta2
  TLorentzVector parton1Whad_frame(parton1); //first parton boosted to hadronic frame XXX
  parton1Whad_frame.Boost(boostWhad);

  parton1_in_had_px = parton1Whad_frame.X(); // new
  parton1_in_had_py = parton1Whad_frame.Y();
  parton1_in_had_pz = parton1Whad_frame.Z();
  
  if (parton1_in_had_px != parton1_in_had_px) parton1_in_had_px = -50.0;
  if (parton1_in_had_py != parton1_in_had_py) parton1_in_had_py = -50.0;
  if (parton1_in_had_pz != parton1_in_had_pz) parton1_in_had_pz = -60.0;
  
  TVector3 parton1Whad_frame3vec(parton1_in_had_px, parton1_in_had_py, parton1_in_had_pz); //3-vector of 1st parton in hadronic frame XXX
  TVector3 parton1Whad_frame_unit3vec = parton1Whad_frame3vec.Unit(); //unit 3-vector of 1st parton in hadronic frame XXX

  
  Double_t x_m21 = parton1Whad_frame_unit3vec.Dot(unitx_2);
  Double_t y_m21 = parton1Whad_frame_unit3vec.Dot(unity_2);
  Double_t z_m21 = parton1Whad_frame_unit3vec.Dot(unitz_2); //dot products of unit 3-vector of 1st parton in hadronic frame and lepton units in parton frame XXX

  parton1_dotted_x = x_m21; // new
  parton1_dotted_y = y_m21;
  parton1_dotted_z = z_m21;
  
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21); //3-vector made of dot products of 1st parton
  costheta2 = M21_Z2frame.CosTheta(); //ANGLE THETA2

  ////////////////////////// now a little bit more 
  
  // calculate phi
  //calculating phi_n
  TLorentzVector lepton_sumWWframe(lepton_sum); //sum of leptons, boosted to WW frame XXX
  TLorentzVector lepton1WWframe(lepton1); //1st lepton, boosted to WW frame XXX
  lepton_sumWWframe.Boost(boostWWframe);
  lepton1WWframe.Boost(boostWWframe);

  
  lepton_sumWWframe_X = lepton_sumWWframe.X(); //new
  lepton_sumWWframe_Y = lepton_sumWWframe.Y();
  lepton_sumWWframe_Z = lepton_sumWWframe.Z();
  
  lepton1WWframe_X = lepton1WWframe.X(); //new
  lepton1WWframe_Y = lepton1WWframe.Y();
  lepton1WWframe_Z = lepton1WWframe.Z();
  
  if (lepton_sumWWframe_X != lepton_sumWWframe_X) lepton_sumWWframe_X = -5000.0;
  if (lepton_sumWWframe_Y != lepton_sumWWframe_Y) lepton_sumWWframe_Y = -5000.0;
  if (lepton_sumWWframe_Z != lepton_sumWWframe_Z) lepton_sumWWframe_Z = -5000.0;

  if (lepton1WWframe_X != lepton1WWframe_X) lepton1WWframe_X = -5000.0;
  if (lepton1WWframe_Y != lepton1WWframe_Y) lepton1WWframe_Y = -5000.0;
  if (lepton1WWframe_Z != lepton1WWframe_Z) lepton1WWframe_Z = -5000.0;
  
  TVector3 lepton_sumWWframe_3vec(lepton_sumWWframe_X, lepton_sumWWframe_Y, lepton_sumWWframe_Z); //new
  TVector3 lepton1WWframe_3vec(lepton1WWframe_X, lepton1WWframe_Y, lepton1WWframe_Z);  //new
  
  //TVector3 lepton_sumWWframe_unit3vec = lepton_sumWWframe.Vect().Unit(); //unit 3-vector of sum of leptons, boosted to WW frame XXX
  TVector3 lepton_sumWWframe_unit3vec = lepton_sumWWframe_3vec.Unit(); //unit 3-vector of sum of leptons, boosted to WW frame XXX
  
  //TVector3 lepton1WWframe_unit3vec = lepton1WWframe.Vect().Unit(); //unit 3-vector of 1st lepton, boosted to WW frame XXX
  TVector3 lepton1WWframe_unit3vec = lepton1WWframe_3vec.Unit(); //unit 3-vector of 1st lepton, boosted to WW frame XXX
  
  TVector3 n_unitz_1(lepton_sumWWframe_unit3vec); //unit 3-vector of sum of leptons, boosted to WW frame XXX (yes, it's there twice)
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
 
  TVector3 n_unity_1 = n_unitz_1.Cross(lepton1WWframe_unit3vec); //cross product of 1st lepton and sum of leptons, in WW frame XXX 
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1); //cross product lepton sum and above cross product, all in WW frame XXX

  ////unit1n: z = lepton sum; y = z cross lepton1; x = y cross z

  //HERE?

  TLorentzVector parton1WWframe(parton1); //1st parton, boosted to WW frame XXX
  parton1WWframe.Boost(boostWWframe);

  parton1WWframe_X = parton1WWframe.X(); //new
  parton1WWframe_Y = parton1WWframe.Y();
  parton1WWframe_Z = parton1WWframe.Z();
  
  if (parton1WWframe_X != parton1WWframe_X) parton1WWframe_X = -5000.0;
  if (parton1WWframe_Y != parton1WWframe_Y) parton1WWframe_Y = -5000.0;
  if (parton1WWframe_Z != parton1WWframe_Z) parton1WWframe_Z = -5000.0;
  
  TVector3 parton1WWframe_3vec(parton1WWframe_X, parton1WWframe_Y, parton1WWframe_Z);  //new
  
  //TVector3 parton1WWframe_unit3vec = parton1WWframe.Vect().Unit(); //unit 3-vector of 1st parton, boosted to WW frame XXX
  TVector3 parton1WWframe_unit3vec = parton1WWframe_3vec.Unit(); //unit 3-vector of 1st parton, boosted to WW frame XXX
  //rotate into other plane
 

  ///////-----------------new way of calculating phi-----------------///////
  
  /// and then calculate phi1
  TVector3 z_hat(0.0, 0.0, 1.0);
  TVector3 z_component_n1_units(z_hat.Dot(n_unitx_1), z_hat.Dot(n_unity_1), z_hat.Dot(n_unitz_1)); //n_units dotted into z-axis??? XXX
  // negative sign is for arrow convention in paper
  phi1 = (z_component_n1_units.Phi()); //ANGLE PHI1

  complicated1_px = z_hat.Dot(n_unitx_1); //new
  complicated1_py = z_hat.Dot(n_unity_1);
  complicated1_pz = z_hat.Dot(n_unitz_1);

  // and the calculate phi2 
  TLorentzVector parton_sumWWframe(parton_sum); //sum of partons in WW frame XXX
  parton_sumWWframe.Boost(boostWWframe);

  parton_sumWWframe_X = parton_sumWWframe.X(); //new HERE NEXT
  parton_sumWWframe_Y = parton_sumWWframe.Y();
  parton_sumWWframe_Z = parton_sumWWframe.Z();
  
  if (parton_sumWWframe_X != parton_sumWWframe_X) parton_sumWWframe_X = -5000.0;
  if (parton_sumWWframe_Y != parton_sumWWframe_Y) parton_sumWWframe_Y = -5000.0;
  if (parton_sumWWframe_Z != parton_sumWWframe_Z) parton_sumWWframe_Z = -5000.0;
  
  TVector3 parton_sumWWframe_3vec(parton_sumWWframe_X, parton_sumWWframe_Y, parton_sumWWframe_Z); //new
  
  //TVector3 parton_sumWWframe_unit3vec = parton_sumWWframe.Vect().Unit(); //unit 3-vector of sum of partons in WW frame XXX
  TVector3 parton_sumWWframe_unit3vec = parton_sumWWframe_3vec.Unit(); //unit 3-vector of sum of partons in WW frame XXX 
  
  TVector3 n_unitz_2(parton_sumWWframe_unit3vec); //unit 3-vector of sum of partons in WW frame XXX (yes, it's there twice)
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
 
  TVector3 n_unity_2 = n_unitz_2.Cross(parton1WWframe_unit3vec); //sum of partons crossed into 2nd lepton, in WW frame XXX
  TVector3 n_unitx_2 = n_unity_2.Cross(n_unitz_2); //above cross product crossed into 2nd lepton, all in WW frame XXX

  //unit2n: z = parton sum; y = z cross parton1; x = y cross z

  TVector3 z_component_n2_units(z_hat.Dot(n_unitx_2), z_hat.Dot(n_unity_2), z_hat.Dot(n_unitz_2)); //above wierd thing in units of whatever XXX
  phi2 = (z_component_n2_units.Phi()); //ANGLE PHI2

  complicated2_px = z_hat.Dot(n_unitx_2); //new
  complicated2_py = z_hat.Dot(n_unity_2);
  complicated2_pz = z_hat.Dot(n_unitz_2);

}



//define this as a plug-in
DEFINE_FWK_MODULE(ntuplizer);
