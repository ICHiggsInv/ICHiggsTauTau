#ifndef ICHiggsTauTau_Module_HhhSync_h
#define ICHiggsTauTau_Module_HhhSync_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "TTree.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SVFitService.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"

#include <string>

namespace ic {

class HhhSync : public ModuleBase {
 private:
  std::string output_name_;


  SVFitService sv_service_;
  CLASS_MEMBER(HhhSync, bool, is_embedded)
  CLASS_MEMBER(HhhSync, std::string, select_category)
  CLASS_MEMBER(HhhSync, std::string, met_label)
  CLASS_MEMBER(HhhSync, int, select_sel_mode)
  CLASS_MEMBER(HhhSync, double, jet_eta)

  TFile *lOFile;
  TTree *lOTree;

  int   lRun;
  int   lLumi;
  int   lEvt;

  //Event Variables
  int   lNPV;
  int   lNPU;
  float lRho;
  
  //Event Weights
  float lMCWeight;
  float lPUWeight;
  float lFakeWeight;

  float lTrigweight_1;
  float lTrigweight_2;
  float lIdweight_1;
  float lIdweight_2;
  float lIsoweight_1;
  float lIsoweight_2;

  float lEffWeight;
  float lWeight;
  float lEmbeddedWeight;
  float lSignalWeight;
  //SV Fit variables
  float lMSV;
  float lPtSV;
  float lEtaSV;
  float lPhiSV;
  float lMVis;
  float lMSVUp;
  float lMSVDown;
  
  ///First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
  float lPt1;
  float lPhi1;
  float lEta1;
  float lM1;
  int   lq1;
  float lIso1;
  float lMVA1;
  float lD01;
  float lDZ1;
  bool  lPassId1;
  bool  lPassIso1;
  float lMt1;

  float lPtTT;

  ///Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
  float lPt2;
  float lPhi2;
  float lEta2;
  float lM2;
  int   lq2;
  float lIso2;
  float lD02;
  float lDZ2;
  float l3Hits_2;
  float lagainstElectronMVA3raw_2;
  float lbyIsolationMVA2raw_2;
  float lagainstMuonLoose2_2;
  float lagainstMuonMedium2_2;
  float lagainstMuonTight2_2;
  float lMVA2;
  bool  lPassId2;
  bool  lPassIso2;
  float lMt2;
  
  //Met related variables
  float lMet;
  float lMetPhi;

  float lL1Met;
  float lL1MetPhi;
  float lL1MetCorr;
  
  float lCaloMet;
  float lCaloMetPhi;
  float lCaloMetCorr;
  float lCaloMetPhiCorr;
  
  float lMVAMet;
  float lMVAMetPhi;
  float lPZetaVis;
  float lPZetaMiss;
  //MET covariance matrices
  float lMetCov00;
  float lMetCov01;
  float lMetCov10;
  float lMetCov11;
  //MVAMet covariance matrices
  float lMVACov00;
  float lMVACov01;
  float lMVACov10;
  float lMVACov11;

  //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
  float lJPt1;
  float lJEta1;
  float lJPhi1;
  float lJPtRaw1;
  float lJPtUnc1;
  float lJMVA1;
  float lLRM1;
  int lCTM1;
  bool  lJPass1;

  //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
  float lJPt2;
  float lJEta2;
  float lJPhi2;
  float lJPtRaw2;
  float lJPtUnc2;
  float lJMVA2;
  float lLRM2;
  int lCTM2;
  bool  lJPass2;
  
  //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
  float lBTagPt1;
  float lBTagEta1;
  float lBTagPhi1;
  float lBTagCSV1;
  
  float lBTagPt2;
  float lBTagEta2;
  float lBTagPhi2;
  float lBTagCSV2;
  
  float lBTagPt3;
  float lBTagEta3;
  float lBTagPhi3;
  float lBTagCSV3;
  
  float lMbb;
  float lMttbb;
  
  //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets 
  float lMJJ;
  float lJDEta;
  int   lNJetInGap;
  float lMVA;
  
  //Variables that go into the VBF MVA
  float lJDPhi;
  float lDiJetPt;
  float lDiJetPhi;
  float lHDJetPhi;
  float lVisJetEta;
  float lPtVis;
  
  //number of btags passing btag id ( pt > 20 )
  int   lNBTag;

  //number of jets passing jet id ( pt > 30 )
  int   lNJets;
  int   lNJetsPt20;
      
  
  float em_gf_mva_;
  float em_vbf_mva_;

  ic::channel channel_;
  TMVA::Reader *reader;
  std::vector<float> vbfvars;


 public:
  HhhSync(std::string const& name, std::string const& output_name, ic::channel channel);
  virtual ~HhhSync();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();



};

}

#endif
