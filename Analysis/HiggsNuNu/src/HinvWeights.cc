#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWeights.h"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerObject.hh"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"


namespace ic {

  HinvWeights::HinvWeights(std::string const& name) : ModuleBase(name),
    mc_(mc::summer12_53X),
    era_(era::data_2012_moriond) {
    do_trg_weights_     = false;
    trg_applied_in_mc_  = false;
    do_idiso_weights_   = false;
    do_w_soup_          = false;
    do_dy_soup_          = false;
    hist_trigSF_METL1 = 0;
    hist_trigSF_METHLT = 0;
    hist_trigSF_MjjHLT = 0;
    hist_trigSF_JetHLT = 0;

  }

  HinvWeights::~HinvWeights() {
    ;
  }

  int HinvWeights::PreAnalysis() {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "PreAnalysis Info for HinvWeights" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Era: " << Era2String(era_) << std::endl;
    std::cout << "MC: " << MC2String(mc_) << std::endl;
    std::cout << "Do Trg Weights?: \t\t" << do_trg_weights_ << std::endl;
    std::cout << "Trg Sel Applied?: \t\t" << trg_applied_in_mc_ << std::endl;
    std::cout << "Do ID & iso weights?: \t\t" << do_idiso_weights_ << std::endl;

    if (do_w_soup_) {
      std::cout << "Making W Soup:" << std::endl;
      std::cout << "nInc = " << n_inc_ << std::endl;
      w1_ = (n_inc_*f1_) / ( (n_inc_*f1_) + n1_ );
      w2_ = (n_inc_*f2_) / ( (n_inc_*f2_) + n2_ );
      w3_ = (n_inc_*f3_) / ( (n_inc_*f3_) + n3_ );
      w4_ = (n_inc_*f4_) / ( (n_inc_*f4_) + n4_ );
      std::cout << "f1 = " << f1_ << "\t" << "n1 = " << n1_ << "\t" << "w1 = " << w1_ << std::endl;
      std::cout << "f2 = " << f2_ << "\t" << "n2 = " << n2_ << "\t" << "w2 = " << w2_ << std::endl;
      std::cout << "f3 = " << f3_ << "\t" << "n3 = " << n3_ << "\t" << "w3 = " << w3_ << std::endl;
      std::cout << "f4 = " << f4_ << "\t" << "n4 = " << n4_ << "\t" << "w4 = " << w4_ << std::endl;
    }
    if (do_dy_soup_) {
      std::cout << "Making DY Soup:" << std::endl;
      std::cout << "nInc = " << zn_inc_ << std::endl;
      zw1_ = (zn_inc_*zf1_) / ( (zn_inc_*zf1_) + zn1_ );
      zw2_ = (zn_inc_*zf2_) / ( (zn_inc_*zf2_) + zn2_ );
      zw3_ = (zn_inc_*zf3_) / ( (zn_inc_*zf3_) + zn3_ );
      zw4_ = (zn_inc_*zf4_) / ( (zn_inc_*zf4_) + zn4_ );
      std::cout << "f1 = " << zf1_ << "\t" << "n1 = " << zn1_ << "\t" << "w1 = " << zw1_ << std::endl;
      std::cout << "f2 = " << zf2_ << "\t" << "n2 = " << zn2_ << "\t" << "w2 = " << zw2_ << std::endl;
      std::cout << "f3 = " << zf3_ << "\t" << "n3 = " << zn3_ << "\t" << "w3 = " << zw3_ << std::endl;
      std::cout << "f4 = " << zf4_ << "\t" << "n4 = " << zn4_ << "\t" << "w4 = " << zw4_ << std::endl;
    }

    //get trigger scale factor histograms from file
    triggerSF_ = new TFile("data/scale_factors/DataMCWeight_53X_v1.root");
    hist_trigSF_METL1 = (TH1F*)gDirectory->Get("METL1");
    hist_trigSF_METHLT = (TH1F*)gDirectory->Get("METHLT");
    hist_trigSF_MjjHLT = (TH1F*)gDirectory->Get("MjjHLT");
    hist_trigSF_JetHLT = (TH1F*)gDirectory->Get("JetHLT");

    std::cout << " -- Content of histogram METL1 : " << std::endl;
    for (int i(0); i< hist_trigSF_METL1->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_METL1->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_METL1->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_METL1->GetBinContent(i)
		<< std::endl;
    }

    std::cout << " -- Content of histogram METHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_METHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_METHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_METHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_METHLT->GetBinContent(i);
		
      //change first bins to "1": no data/MC scale factor applied...
      if (i>0 && hist_trigSF_METHLT->GetBinContent(i) == 0) {
	hist_trigSF_METHLT->SetBinContent(i,1);
	std::cout << " -> changed to " <<  hist_trigSF_METHLT->GetBinContent(i);
      }
      std::cout	<< std::endl;
      
    }

    std::cout << " -- Content of histogram MjjHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_MjjHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_MjjHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_MjjHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_MjjHLT->GetBinContent(i);
      std::cout	<< std::endl;
    }

    std::cout << " -- Content of histogram JetHLT : " << std::endl;
    for (int i(0); i< hist_trigSF_JetHLT->GetNbinsX()+2; ++i){
      std::cout << " -- bin " << i << " [" 
		<< hist_trigSF_JetHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
		<< hist_trigSF_JetHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
		<<  hist_trigSF_JetHLT->GetBinContent(i);

      //change first bins to "1": no data/MC scale factor applied...
      if (i>0 && hist_trigSF_JetHLT->GetBinContent(i) == 0) {
	hist_trigSF_JetHLT->SetBinContent(i,1);
	std::cout << " -> changed to " <<  hist_trigSF_JetHLT->GetBinContent(i);
      }
      std::cout	<< std::endl;
    }


    return 0;
  }

  int HinvWeights::Execute(TreeEvent *event) {

    //std::vector<CompositeCandidate *> const& dilepton = event->GetPtrVec<CompositeCandidate>("emtauCandidates");

    //double weight = 1.0;
    EventInfo * eventInfo = event->GetPtr<EventInfo>("eventInfo");

 //    if (do_btag_weight_) {
//       std::vector<PFJet*> jets = event->GetPtrVec<PFJet>("pfJetsPFlow"); // Make a copy of the jet collection
//       ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
//       //double no_btag_weight = btag_weight.GetWeight(jets, "CSVM", 0, 0, is_2012_);
//       //double inclusive_btag_weight = btag_weight.GetWeight(jets, "CSVM", 1, 99, is_2012_);
//       double no_btag_weight = 1.0;
//       double inclusive_btag_weight = 1.0;
//       btag_weight.ReTag(jets, mc_ == mc::summer12_53X);
//       event->Add("no_btag_weight", no_btag_weight);
//       event->Add("inclusive_btag_weight", inclusive_btag_weight);
//     }

    if (do_trg_weights_) {
      //get METnoMuons:
      Met const* metNoMu = event->GetPtr<Met>("metNoMuons");

      double lValue = metNoMu->pt();
      double lMax = hist_trigSF_METL1->GetXaxis()->GetBinCenter(hist_trigSF_METL1->GetNbinsX());
      double lMin = hist_trigSF_METL1->GetXaxis()->GetBinCenter(1);
      if (lValue > lMax)  lValue = lMax;
      if (lValue < lMin)  lValue = lMin;
      int lBin = hist_trigSF_METL1->GetXaxis()->FindFixBin(lValue);
      double metl1 = hist_trigSF_METL1->GetBinContent(lBin);
      eventInfo->set_weight("trig_metL1",metl1);
      //std::cout << " -- MET L1 " << lValue << " " << metl1 << std::endl;

      lMax = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(hist_trigSF_METHLT->GetNbinsX());
      lMin = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(1);
      if (lValue > lMax)  lValue = lMax;
      if (lValue < lMin)  lValue = lMin;
      lBin = hist_trigSF_METHLT->GetXaxis()->FindFixBin(lValue);
      double methlt = hist_trigSF_METHLT->GetBinContent(lBin);
      eventInfo->set_weight("trig_metHLT",methlt);
      //std::cout << " -- MET HLT " << lValue << " " << methlt << std::endl;

      //get 2 leading jets
      std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjCandidates");
      if (dijet_vec.size() > 0) {
	
	CompositeCandidate const* dijet = dijet_vec.at(0);

	Candidate const* jet1 = dijet->GetCandidate("jet1");
	Candidate const* jet2 = dijet->GetCandidate("jet2");

	lValue = dijet->M();
	lMax = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(hist_trigSF_MjjHLT->GetNbinsX());
	lMin = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(1);
	if (lValue > lMax)  lValue = lMax;
	if (lValue < lMin)  lValue = lMin;
	lBin = hist_trigSF_MjjHLT->GetXaxis()->FindFixBin(lValue);
	double mjjhlt = hist_trigSF_MjjHLT->GetBinContent(lBin);
	eventInfo->set_weight("trig_mjjHLT",mjjhlt);
	//std::cout << " -- Mjj HLT " << lValue << " " << mjjhlt << std::endl;

	lMax = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(hist_trigSF_JetHLT->GetNbinsX());
	lMin = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(1);
	lValue = jet1->pt();
	if (lValue < lMin)  lValue = lMin;
	if (lValue > lMax)  lValue = lMax;
	lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(lValue);
	double jet1hlt = hist_trigSF_JetHLT->GetBinContent(lBin);
	eventInfo->set_weight("trig_jet1HLT",jet1hlt);
	//std::cout << " -- Jet1 HLT " << lValue << " " << jet1hlt << std::endl;

	lValue = jet2->pt();
	if (lValue > lMax)  lValue = lMax;
	if (lValue < lMin)  lValue = lMin;
	lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(lValue);
	double jet2hlt = hist_trigSF_JetHLT->GetBinContent(lBin);
	eventInfo->set_weight("trig_jet2HLT",jet2hlt);
	//std::cout << " -- Jet2 HLT " << lValue << " " << jet2hlt << std::endl;

      }



      //weight *= (ele_trg * tau_trg);
      //event->Add("trigweight_1", ele_trg);
      //event->Add("trigweight_2", tau_trg);
    }

    //eventInfo->set_weight("lepton", weight);

    bool zeroParton = false;

    if (do_w_soup_) {
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      bool count_jets = false;
      unsigned partons = 0;
      for (unsigned i = 0; i < parts.size(); ++i) {
        if (parts[i]->status() != 3) continue;
        unsigned id = abs(parts[i]->pdgid());
        if (count_jets) { 
          if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21) partons++;
        }
        if (id == 24) count_jets = true; 
      }
      if (partons > 4) {
        std::cerr << "Error making soup, event has " << partons << " partons!" << std::endl;
        throw;
      }
      if (partons == 1) eventInfo->set_weight("wsoup", w1_);
      if (partons == 2) eventInfo->set_weight("wsoup", w2_);
      if (partons == 3) eventInfo->set_weight("wsoup", w3_);
      if (partons == 4) eventInfo->set_weight("wsoup", w4_);

      if (partons == 0) zeroParton = true;

    }

    if (do_dy_soup_) {
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      bool count_jets = false;
      unsigned partons = 0;
      for (unsigned i = 0; i < parts.size(); ++i) {
        // std::cout << i << "\t" << parts[i]->status() << "\t" << parts[i]->pdgid() << "\t" << parts[i]->vector() << std::endl;
        if (parts[i]->status() != 3) continue;
        unsigned id = abs(parts[i]->pdgid());
        if (count_jets) { 
          if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21) partons++;
        }
        if (id == 23) count_jets = true; 
      }
      if (partons > 4) {
        std::cerr << "Error making soup, event has " << partons << " partons!" << std::endl;
        throw;
      }
      if (partons == 1) eventInfo->set_weight("dysoup", zw1_);
      if (partons == 2) eventInfo->set_weight("dysoup", zw2_);
      if (partons == 3) eventInfo->set_weight("dysoup", zw3_);
      if (partons == 4) eventInfo->set_weight("dysoup", zw4_);
      if (partons == 0) zeroParton = true;
    }

    event->Add("NoParton",zeroParton);

    return 0;
  }

  int HinvWeights::PostAnalysis() {
 
    return 0;
  }

  void HinvWeights::PrintInfo() {
    ;
  }

  void HinvWeights::SetWTargetFractions(double f0, double f1, double f2, double f3, double f4) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
    f3_ = f3;
    f4_ = f4;

  }
  void HinvWeights::SetWInputYields(double n_inc, double n1, double n2, double n3, double n4) {
    n_inc_ = n_inc;
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
  }

  void HinvWeights::SetDYTargetFractions(double zf0, double zf1, double zf2, double zf3, double zf4) {
    zf0_ = zf0;
    zf1_ = zf1;
    zf2_ = zf2;
    zf3_ = zf3;
    zf4_ = zf4;

  }
  void HinvWeights::SetDYInputYields(double zn_inc, double zn1, double zn2, double zn3, double zn4) {
    zn_inc_ = zn_inc;
    zn1_ = zn1;
    zn2_ = zn2;
    zn3_ = zn3;
    zn4_ = zn4;
  }

  double HinvWeights::Efficiency(double m, double m0, double sigma, double alpha,
    double n, double norm)
  {
    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) sigma);
    double t = (m - m0)/sig;
    if(alpha < 0)
      t = -t;
    double absAlpha = fabs(alpha/sig);
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = absAlpha - n/absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
    double area = leftArea + rightArea;
    if( t <= absAlpha ){
      arg = t / sqrt2;
      if(arg > 5.) ApproxErf = 1;
      else if (arg < -5.) ApproxErf = -1;
      else ApproxErf = TMath::Erf(arg);
      return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
    }
    else{
      return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
        1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
    }
  }
}