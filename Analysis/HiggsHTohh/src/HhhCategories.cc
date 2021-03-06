#include "UserCode/ICHiggsTauTau/Analysis/HiggsHTohh/interface/HhhCategories.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"

#include "TMVA/Reader.h"
#include "TVector3.h"
#include "boost/format.hpp"
#include "TMath.h"

namespace ic {

  HhhCategories::HhhCategories(std::string const& name) : ModuleBase(name), 
      channel_(channel::et), 
      era_(era::data_2012_moriond),
      strategy_(strategy::paper2013) {
      ditau_label_ = "emtauCandidates";
      met_label_ = "pfMVAMet";
      mass_shift_ = 1.0;
      fs_ = NULL;
      write_tree_ = true;
      write_plots_ = true;
      experimental_ = false;
  }

  HhhCategories::~HhhCategories() {
    ;
  }

  void HhhCategories::InitSelection(std::string const& selection) {
    selections_[selection] = false;
  }
  
  void HhhCategories::InitCategory(std::string const& category) {
    categories_[category] = false;
    InitMassPlots(category);
  }

  int HhhCategories::PreAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HhhCategories" << std::endl;
    std::cout << "-------------------------------------" << std::endl;    
    if (fs_) {
      std::cout << boost::format(param_fmt()) % "channel"         % Channel2String(channel_);
      std::cout << boost::format(param_fmt()) % "strategy"        % Strategy2String(strategy_);
      std::cout << boost::format(param_fmt()) % "era"             % Era2String(era_);
      std::cout << boost::format(param_fmt()) % "dilepton_label"  % ditau_label_;
      std::cout << boost::format(param_fmt()) % "met_label"       % met_label_;
      std::cout << boost::format(param_fmt()) % "mass_shift"      % mass_shift_;
      std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;
      std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
      std::cout << boost::format(param_fmt()) % "experimental"    % experimental_;

      if (write_tree_) {
        outtree_ = fs_->make<TTree>("ntuple","ntuple");
        outtree_->Branch("wt",                &wt_);
        outtree_->Branch("wt_ggh_pt_up",      &wt_ggh_pt_up_);
        outtree_->Branch("wt_ggh_pt_down",    &wt_ggh_pt_down_);
        outtree_->Branch("wt_tau_fake_up",    &wt_tau_fake_up_);
        outtree_->Branch("wt_tau_fake_down",  &wt_tau_fake_down_);
        outtree_->Branch("wt_tquark_up",      &wt_tquark_up_);
        outtree_->Branch("wt_tquark_down",    &wt_tquark_down_);
        outtree_->Branch("wt_tau_id_up",      &wt_tau_id_up_);
        outtree_->Branch("wt_tau_id_down",    &wt_tau_id_down_);
        outtree_->Branch("os",                &os_);
        outtree_->Branch("n_vtx",             &n_vtx_);
        outtree_->Branch("m_sv",              &m_sv_);
        outtree_->Branch("m_vis",             &m_vis_);
        outtree_->Branch("pt_h",              &pt_h_);
        outtree_->Branch("pt_tt",             &pt_tt_);
        outtree_->Branch("mt_1",              &mt_1_);
        outtree_->Branch("pzeta",             &pzeta_);
        outtree_->Branch("pt_1",              &pt_1_);
        outtree_->Branch("pt_2",              &pt_2_);
        outtree_->Branch("eta_1",             &eta_1_);
        outtree_->Branch("eta_2",             &eta_2_);
        outtree_->Branch("iso_2",             &iso_2_);
        outtree_->Branch("z_2",               &z_2_);
        outtree_->Branch("m_2",               &m_2_);
        outtree_->Branch("met",               &met_);
        outtree_->Branch("met_phi",           &met_phi_);
        outtree_->Branch("tau_decay_mode",    &tau_decay_mode_);
        outtree_->Branch("n_jets",            &n_jets_);
        outtree_->Branch("n_lowpt_jets",      &n_lowpt_jets_);
        outtree_->Branch("n_bjets",           &n_bjets_);
        outtree_->Branch("n_prebjets",        &n_prebjets_);
        outtree_->Branch("n_loose_bjets",     &n_loose_bjets_);
        outtree_->Branch("n_jetsingap",       &n_jetsingap_);
        outtree_->Branch("jpt_1",             &jpt_1_);
        outtree_->Branch("j1_dm",             &j1_dm_);
        outtree_->Branch("jpt_2",             &jpt_2_);
        outtree_->Branch("jeta_1",            &jeta_1_);
        outtree_->Branch("jeta_2",            &jeta_2_);
        outtree_->Branch("bpt_1",             &bpt_1_);
        outtree_->Branch("beta_1",            &beta_1_);
        outtree_->Branch("bcsv_1",            &bcsv_1_);
        outtree_->Branch("prebjetpt_1",             &prebjetpt_1_);
        outtree_->Branch("prebjetEt_1",             &prebjetEt_1_);
        outtree_->Branch("prebjeteta_1",            &prebjeteta_1_);
        outtree_->Branch("prebjetbcsv_1",            &prebjetbcsv_1_);
        outtree_->Branch("prebjet1_dm",             &prebjet1_dm_);
        outtree_->Branch("prebjetpt_2",             &prebjetpt_2_);
        outtree_->Branch("prebjeteta_2",            &prebjeteta_2_);
        outtree_->Branch("prebjetbcsv_2",            &prebjetbcsv_2_);
        outtree_->Branch("E_1",             &E_1_);
        outtree_->Branch("mjj",               &mjj_);
        outtree_->Branch("mjj_h",               &mjj_h_);
        outtree_->Branch("mjj_tt",               &mjj_tt_);
        outtree_->Branch("jdeta",             &jdeta_);
        outtree_->Branch("prebjet_mjj",               &prebjet_mjj_);
        outtree_->Branch("prebjet_dphi",               &prebjet_dphi_);
        outtree_->Branch("prebjet_deta",             &prebjet_deta_);
        outtree_->Branch("prebjet_dtheta",             &prebjet_dtheta_);
        outtree_->Branch("prebjet_1_met_dphi",             &prebjet_1_met_dphi_);
        outtree_->Branch("prebjet_1_met_dtheta",             &prebjet_1_met_dtheta_);
        outtree_->Branch("prebjet_1_lep1_dphi",             &prebjet_1_lep1_dphi_);
        outtree_->Branch("prebjet_1_lep1_dtheta",             &prebjet_1_lep1_dtheta_);
        outtree_->Branch("prebjet_1_lep1_m",             &prebjet_1_lep1_m_);
        outtree_->Branch("jet_1_met_dphi",             &jet_1_met_dphi_);
        outtree_->Branch("jet_1_met_dtheta",             &jet_1_met_dtheta_);
        outtree_->Branch("mjj_lowpt",         &mjj_lowpt_);
        outtree_->Branch("jdeta_lowpt",       &jdeta_lowpt_);
        outtree_->Branch("n_jetsingap_lowpt", &n_jetsingap_lowpt_);
        outtree_->Branch("l1_met",            &l1_met_);
        outtree_->Branch("calo_nohf_met",     &calo_nohf_met_);
        if (channel_ == channel::em) {
          outtree_->Branch("em_gf_mva",         &em_gf_mva_);
          // outtree_->Branch("em_vbf_mva",        &em_vbf_mva_);
          outtree_->Branch("pzetavis",          &pzetavis_);
          outtree_->Branch("pzetamiss",         &pzetamiss_);
          outtree_->Branch("mt_ll",             &mt_ll_);
          outtree_->Branch("emu_dphi",          &emu_dphi_);
          outtree_->Branch("emu_csv",           &emu_csv_);
          outtree_->Branch("emu_dxy_1",         &emu_dxy_1_);
          outtree_->Branch("emu_dxy_2",         &emu_dxy_2_);
          //outtree_->Branch("emu_dxy_fromref",&emu_dxy_fromref_);
          //outtree_->Branch("emu_dxy", &emu_dxy_);
          //outtree_->Branch("em_gf_mva",&em_gf_mva_);
          outtree_->Branch("em_gf_mva_bdtg",&em_gf_mva_bdtg_);
          outtree_->Branch("em_gf_mva_bdt",&em_gf_mva_bdt_);
          //outtree_->Branch("bdt",&bdt_);
          //outtree_->Branch("emu_dca_1",&emu_dca_1_);
          //outtree_->Branch("emu_dca_2",&emu_dca_2_);
        }
      }
    }
    TH1F::SetDefaultSumw2();
    misc_plots_ = new DynamicHistoSet(fs_->mkdir("misc_plots"));
    misc_2dplots_ = new Dynamic2DHistoSet(fs_->mkdir("misc_2dplots"));
    misc_2dplots_->Create("jpt_vs_hpt", 50, 0, 200, 50, 0, 200);
    misc_2dplots_->Create("hpt_vs_pt_tt", 50, 0, 200, 50, 0, 200);
    misc_2dplots_->Create("gen_vs_hpt", 50, 0, 200, 50, 0, 200);
    misc_2dplots_->Create("gen_vs_pt_tt", 50, 0, 200, 50, 0, 200);
    misc_2dplots_->Create("m_gen_m_sv_nobtag", 75, 0, 1500, 75, 0, 1500);
    misc_2dplots_->Create("m_gen_m_sv_btag", 75, 0, 1500, 75, 0, 1500);
    misc_2dplots_->Create("m_gen_m_vis", 75, 0, 1500, 75, 0, 1500);
    misc_2dplots_->Create("orig_met_m_sv", 75, 0, 750, 75, 0, 1500);
    misc_2dplots_->Create("orig_met_m_vis", 75, 0, 750, 75, 0, 1500);
    misc_plots_->Create("M",150,0,1500);

    if (write_plots_) {
      InitSelection("os");
      InitSelection("os_sel");
      InitSelection("os_con");
      InitSelection("os_con_mt_60-120");
      InitSelection("ss");
      InitSelection("ss_sel");
      InitSelection("ss_con");
      InitSelection("ss_con_mt_60-120");

      InitCategory("inclusive");
      InitCoreControlPlots("inclusive");

      InitCategory("vbf");
      InitCategory("vbf_no_cjv");
      InitCategory("vbf_loose");
      InitCategory("vbf_loose_jets20");

      InitCategory("twojet");
      InitCoreControlPlots("twojet");

      InitCategory("1jet");
      InitCoreControlPlots("1jet");

      InitCategory("1jet_high");
      InitCategory("1jet_low");
      InitCategory("1jet_low_nometcut");
      InitCategory("0jet_high");
      InitCategory("0jet_low");

      InitCategory("prebtag");
      InitCoreControlPlots("prebtag");

      InitCategory("sasha");
      InitCategory("presasha");
      
      InitCategory("btag");
      InitCoreControlPlots("btag");

      InitCategory("btag_low");
      InitCategory("btag_low_loose");

      InitCategory("btag_high");
      InitCategory("btag_high_loose");

      InitCategory("btag_loose");

      InitCategory("nobtag");

      InitCategory("1jet0tag");
      InitCoreControlPlots("1jet0tag");
      InitCategory("1jet1tag");
      InitCoreControlPlots("1jet1tag");
      InitCategory("2jet0tag");
      InitCoreControlPlots("2jet0tag");
      InitCategory("2jet1tag");
      InitCoreControlPlots("2jet1tag");
      InitCategory("2jet2tag");
      InitCoreControlPlots("2jet2tag");
    }
    return 0;
  }

  int HhhCategories::Execute(TreeEvent *event) {

    Reset();

    // Get the objects we need from the event
    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");

    wt_ = eventInfo->total_weight();
    wt_ggh_pt_up_ = 1.0;
    wt_ggh_pt_down_ = 1.0;
    wt_tau_fake_up_ = 1.0;
    wt_tau_fake_down_ = 1.0;
    wt_tquark_up_ = 1.0;
    wt_tquark_down_ = 1.0;
    wt_tau_id_up_ = 1.0;
    wt_tau_id_down_ = 1.0;
    if (event->Exists("wt_ggh_pt_up"))      wt_ggh_pt_up_   = event->Get<double>("wt_ggh_pt_up");
    if (event->Exists("wt_ggh_pt_down"))    wt_ggh_pt_down_ = event->Get<double>("wt_ggh_pt_down");
    if (event->Exists("wt_tau_fake_up"))    wt_tau_fake_up_   = event->Get<double>("wt_tau_fake_up");
    if (event->Exists("wt_tau_fake_down"))  wt_tau_fake_down_ = event->Get<double>("wt_tau_fake_down");
    if (event->Exists("wt_tquark_up"))      wt_tquark_up_   = event->Get<double>("wt_tquark_up");
    if (event->Exists("wt_tquark_down"))    wt_tquark_down_ = event->Get<double>("wt_tquark_down");
    if (event->Exists("wt_tau_id_up"))      wt_tau_id_up_   = event->Get<double>("wt_tau_id_up");
    if (event->Exists("wt_tau_id_down"))    wt_tau_id_down_ = event->Get<double>("wt_tau_id_down");
    
    std::vector<CompositeCandidate *> const& ditau_vec = event->GetPtrVec<CompositeCandidate>(ditau_label_);
    CompositeCandidate const* ditau = ditau_vec.at(0);
    Candidate const* lep1 = ditau->GetCandidate("lepton1");
    Candidate const* lep2 = ditau->GetCandidate("lepton2");
    Met const* met = event->GetPtr<Met>(met_label_);
    std::vector<PFJet*> jets = event->GetPtrVec<PFJet>("pfJetsPFlow");
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::vector<PFJet*> lowpt_jets = jets;
    ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 30.0, 4.7));
    ic::erase_if(lowpt_jets,!boost::bind(MinPtMaxEta, _1, 20.0, 4.7));
    std::vector<PFJet*> prebjets = lowpt_jets;
    ic::erase_if(prebjets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
    std::vector<PFJet*> bjets = prebjets;
    std::vector<PFJet*> loose_bjets = prebjets;
    ic::erase_if(loose_bjets, boost::bind(&PFJet::GetBDiscriminator, _1, "combinedSecondaryVertexBJetTags") < 0.244);
    //Use prebjet collection for candidate jets for h->bb. Sort by CSV discriminator
    std::sort(prebjets.begin(), prebjets.end(), bind(&PFJet::GetBDiscriminator, _1, "combinedSecondaryVertexBJetTags") > bind(&PFJet::GetBDiscriminator, _2, "combinedSecondaryVertexBJetTags"));

    // Instead of changing b-tag value in the promote/demote method we look for a map of bools
    // that say whether a jet should pass the WP or not
    if (event->Exists("retag_result")) {
      auto const& retag_result = event->Get<std::map<std::size_t,bool>>("retag_result"); 
      ic::erase_if(bjets, !boost::bind(IsReBTagged, _1, retag_result));
    } else {
      ic::erase_if(bjets, boost::bind(&PFJet::GetBDiscriminator, _1, "combinedSecondaryVertexBJetTags") < 0.679);
    } 
    
    // Define event properties
    // IMPORTANT: Make sure each property is re-set
    // for each new event
    if (PairOppSign(ditau)) {
      os_ = true;
    } else {
      os_ = false;
    }

    n_vtx_ = eventInfo->good_vertices();

    if (event->Exists("svfitMass")) {
      m_sv_ = event->Get<double>("svfitMass");
    } else {
      m_sv_ = -9999;
    }

    if (event->Exists("svfitHiggs")) {
      pt_h_ = event->Get<Candidate>("svfitHiggs").pt();
    } else {
      pt_h_ = -9999;
    }

    pt_tt_ = (ditau->vector() + met->vector()).pt();
    m_vis_ = ditau->M();
   
    GenParticle const* boson = NULL;
    if (experimental_ && event->Exists("genParticles")) {
      std::vector<GenParticle *> parts = event->GetPtrVec<GenParticle>("genParticles");
      for (unsigned i = 0; i < parts.size(); ++i) {
        unsigned id = abs(parts[i]->pdgid());
        if (parts[i]->status() == 3) {
          if (id == 25 || id == 35 || id == 36) {
            boson = parts[i];
          }
        }
      }
    }

    // This is the HCP hack for the em channel
    // to better align the data with the embedded
    // mass.  
    if (channel_ == channel::em) {
      m_sv_ = m_sv_ * mass_shift_;
      m_vis_ = m_vis_ * mass_shift_;
      em_gf_mva_ = event->Exists("em_gf_mva") ? event->Get<double>("em_gf_mva") : 0.;
			em_gf_mva_bdtg_ = event->Exists("em_gf_mva_bdtg") ? event->Get<double>("em_gf_mva_bdtg") :0.;
			em_gf_mva_bdt_ = event->Exists("em_gf_mva_bdt") ? event->Get<double>("em_gf_mva_bdt") : 0.;
      // em_vbf_mva_ = event->Exists("em_vbf_mva") ? event->Get<double>("em_vbf_mva") : 0.;
    }
    if (event->Exists("mass_scale")) {
      m_sv_ = m_sv_ * event->Get<double>("mass_scale");
      m_vis_ = m_vis_ * event->Get<double>("mass_scale");
    }


    mt_1_ = MT(lep1, met);
    mt_ll_ = MT(ditau, met);
    pzeta_ = PZeta(ditau, met, 0.85);
    pzetavis_ = PZetaVis(ditau);
    pzetamiss_ = PZeta(ditau, met, 0.0);
    emu_dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(lep1->vector(), lep2->vector()));

    pt_1_ = lep1->pt();
    E_1_ = lep1->energy();
    pt_2_ = lep2->pt();
    eta_1_ = lep1->eta();
    eta_2_ = lep2->eta();

    m_2_ = lep2->M();
    met_ = met->pt();
    met_phi_ = met->phi();

    emu_dxy_1_ = 0.0;
    emu_dxy_2_ = 0.0;



    if (strategy_ == strategy::paper2013) {
      if (channel_ == channel::et) {
        Electron const* elec = dynamic_cast<Electron const*>(lep1);
        iso_1_ = PF04IsolationVal(elec, 0.5);
        Tau const* tau = dynamic_cast<Tau const*>(lep2);
        iso_2_ = tau->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      }
      if (channel_ == channel::mt || channel_ == channel::mtmet) {
        Muon const* muon = dynamic_cast<Muon const*>(lep1);
        iso_1_ = PF04IsolationVal(muon, 0.5);
        Tau const* tau = dynamic_cast<Tau const*>(lep2);
        iso_2_ = tau->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      }
      if (channel_ == channel::em) {
        Electron const* elec = dynamic_cast<Electron const*>(lep1);
        iso_1_ = PF04IsolationVal(elec, 0.5);
        Muon const* muon = dynamic_cast<Muon const*>(lep2);
        iso_2_ = PF04IsolationVal(muon, 0.5);
        emu_dxy_1_ = -1. * elec->dxy_vertex();
        emu_dxy_2_ = -1. * muon->dxy_vertex();
      }
    } else {
      iso_1_ = 0.0;
      iso_2_ = 0.0;
    }

    if (strategy_ == strategy::paper2013 && era_ != era::data_2011) {
      auto l1_met = event->GetPtrVec<Candidate>("l1extraMET");
      l1_met_ = l1_met.at(0)->pt();
      auto calo_nohf_met = event->GetPtr<Met>("metNoHF");
      calo_nohf_met_ = calo_nohf_met->pt();
    } else {
      l1_met_ = 0.0;
      calo_nohf_met_ = 0.0;
    }

    Tau const* tau = dynamic_cast<Tau const*>(lep2);
    if (tau) {
      tau_decay_mode_ = tau->decay_mode();
      z_2_ = tau->vz() + (130. / tan(tau->vector().theta()));
    } else {
      tau_decay_mode_ = 0;
      z_2_ = 0.;
    }

    n_jets_ = jets.size();
    n_lowpt_jets_ = lowpt_jets.size();
    n_bjets_ = bjets.size();
    n_prebjets_ = prebjets.size();
    n_loose_bjets_ = loose_bjets.size();
 

    if (n_prebjets_ >= 1) {
      prebjetpt_1_ = prebjets[0]->pt();
      prebjetEt_1_ = std::sqrt(prebjets[0]->pt()*prebjets[0]->pt() + prebjets[0]->M()*prebjets[0]->M());
      prebjeteta_1_ = prebjets[0]->eta();
      prebjet_1_met_dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(prebjets[0]->vector(), met->vector()));
      prebjet_1_met_dtheta_ = std::fabs(prebjets[0]->vector().theta() - met->vector().theta());
      prebjet_1_lep1_dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(prebjets[0]->vector(), lep1->vector()));
      prebjet_1_lep1_dtheta_ = std::fabs(prebjets[0]->vector().theta() -  lep1->vector().theta());
      prebjet_1_lep1_m_ = (prebjets[0]->vector() + lep1->vector()).M();
      std::vector<ic::Tau *> taus = event->GetPtrVec<Tau>("taus");
      std::vector<ic::Jet *> leadjet = { prebjets[0] };
      std::vector<std::pair<ic::Jet *, ic::Tau *>> matches = MatchByDR(leadjet, taus, 0.5, true, true);
      if (matches.size() == 1) {
        prebjet1_dm_ = matches[0].second->decay_mode();
      } else {
        prebjet1_dm_ = -1;
      }
    } else {
      prebjetpt_1_ = -9999;
      prebjetEt_1_ = -9999;
      prebjeteta_1_ = -9999;
      prebjet_1_met_dphi_ = -9999;
      prebjet_1_met_dtheta_ = -9999;
      prebjet_1_lep1_dphi_ = -9999;
      prebjet_1_lep1_dtheta_ = -9999;
      prebjet_1_lep1_m_ = -9999;
    }

    if (n_prebjets_ >= 2) {
      prebjetpt_2_ = prebjets[1]->pt();
      prebjeteta_2_ = prebjets[1]->eta();
      prebjet_mjj_ = (prebjets[0]->vector() + prebjets[1]->vector()).M();
      prebjet_deta_ = fabs(prebjets[0]->eta() - prebjets[1]->eta());
      prebjet_dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(prebjets[0]->vector(), prebjets[1]->vector()));
      prebjet_dtheta_ = std::fabs((prebjets[0]->vector().theta() -  prebjets[1]->vector().theta()));
      mjj_tt_= (prebjets[0]->vector() + prebjets[1]->vector() + ditau->vector() + met->vector()).M();
      if (event->Exists("svfitHiggs")) {
        mjj_h_= (prebjets[0]->vector() + prebjets[1]->vector() + event->Get<Candidate>("svfitHiggs").vector() ).M();
      } else {
        mjj_h_ = -9999;
      }
      //double eta_high = (prebjets[0]->eta() > prebjets[1]->eta()) ? prebjets[0]->eta() : prebjets[1]->eta();
      //double eta_low = (prebjets[0]->eta() > jets[1]->eta()) ? jets[1]->eta() : jets[0]->eta();
      //n_jetsingap_ = 0;
      //if (n_jets_ > 2) {
      //  for (unsigned i = 2; i < jets.size(); ++i) {
      //   if (jets[i]->pt() > 30.0 &&  jets[i]->eta() > eta_low && jets[i]->eta() < eta_high) ++n_jetsingap_;
      //  }
      
    } else {
      prebjetpt_2_ = -9999;
      prebjeteta_2_ = -9999;
      prebjet_mjj_ = -9999;
      prebjet_deta_ = -9999;
      prebjet_dphi_ = -9999;
      prebjet_dtheta_ = -9999;
      mjj_h_ = -9999;
      mjj_tt_ = -9999;
    }
    
    if (n_jets_ >= 1) {
      jpt_1_ = jets[0]->pt();
      jeta_1_ = jets[0]->eta();
      jet_1_met_dphi_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[0]->vector(), met->vector()));
      std::vector<ic::Tau *> taus = event->GetPtrVec<Tau>("taus");
      std::vector<ic::Jet *> leadjet = { jets[0] };
      std::vector<std::pair<ic::Jet *, ic::Tau *>> matches = MatchByDR(leadjet, taus, 0.5, true, true);
      if (matches.size() == 1) {
        j1_dm_ = matches[0].second->decay_mode();
      } else {
        j1_dm_ = -1;
      }
    } else {
      jpt_1_ = -9999;
      jeta_1_ = -9999;
      jet_1_met_dphi_ = -9999;
    }

    if (n_jets_ >= 2) {
      jpt_2_ = jets[1]->pt();
      jeta_2_ = jets[1]->eta();
      mjj_ = (jets[0]->vector() + jets[1]->vector()).M();
      jdeta_ = fabs(jets[0]->eta() - jets[1]->eta());
      double eta_high = (jets[0]->eta() > jets[1]->eta()) ? jets[0]->eta() : jets[1]->eta();
      double eta_low = (jets[0]->eta() > jets[1]->eta()) ? jets[1]->eta() : jets[0]->eta();
      n_jetsingap_ = 0;
      if (n_jets_ > 2) {
        for (unsigned i = 2; i < jets.size(); ++i) {
         if (jets[i]->pt() > 30.0 &&  jets[i]->eta() > eta_low && jets[i]->eta() < eta_high) ++n_jetsingap_;
        }
      }
    } else {
      jpt_2_ = -9999;
      jeta_2_ = -9999;
      mjj_ = -9999;
      jdeta_ = -9999;
      n_jetsingap_ = 9999;
    }


    if (n_lowpt_jets_ >= 2) {
      mjj_lowpt_ = (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).M();
      jdeta_lowpt_ = fabs(lowpt_jets[0]->eta() - lowpt_jets[1]->eta());
      double eta_high = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[0]->eta() : lowpt_jets[1]->eta();
      double eta_low = (lowpt_jets[0]->eta() > lowpt_jets[1]->eta()) ? lowpt_jets[1]->eta() : lowpt_jets[0]->eta();
      n_jetsingap_lowpt_ = 0;
      if (n_lowpt_jets_ > 2) {
        for (unsigned i = 2; i < lowpt_jets.size(); ++i) {
         if (lowpt_jets[i]->pt() > 30.0 &&  lowpt_jets[i]->eta() > eta_low && lowpt_jets[i]->eta() < eta_high) ++n_jetsingap_lowpt_;
        }
      }
    } else {
      mjj_lowpt_ = -9999;
      jdeta_lowpt_ = -9999;
      n_jetsingap_lowpt_ = 9999;
    }

    if (n_bjets_ >= 1) {
      bpt_1_ = bjets[0]->pt();
      beta_1_ = bjets[0]->eta();
    } else {
      bpt_1_ = -9999;
      beta_1_ = -9999;
    }

    if (prebjets.size() >= 1) {
      bcsv_1_ = prebjets[0]->GetBDiscriminator("combinedSecondaryVertexBJetTags");
      prebjetbcsv_1_ = prebjets[0]->GetBDiscriminator("combinedSecondaryVertexBJetTags");

    } else {
      bcsv_1_ = -9999;
      prebjetbcsv_1_ = -9999;
    }
    if (prebjets.size() >= 2) {
      prebjetbcsv_2_ = prebjets[1]->GetBDiscriminator("combinedSecondaryVertexBJetTags");

    } else {
      prebjetbcsv_2_ = -9999;
    }
    emu_csv_ = (bcsv_1_ > 0.244) ? bcsv_1_ : -1.0;

    if ((event->Exists("genParticlesEmbedded") || event->Exists("genParticles")) && os_ && mt_1_ <30.) {
      // recalculate met
      ROOT::Math::PtEtaPhiEVector met_vec = met->vector();
      ROOT::Math::PtEtaPhiEVector ditau_vec = ditau->vector();
      ROOT::Math::PtEtaPhiEVector dimu_vec;
      std::string gen_label = event->Exists("genParticlesEmbedded") ? "genParticlesEmbedded" : "genParticles";
      std::vector<GenParticle *> const& particles = event->GetPtrVec<GenParticle>(gen_label);
      for (unsigned i = 0; i < particles.size(); ++i) {
        if (particles[i]->pdgid() == 23) {
          misc_plots_->Fill("M",particles[i]->vector().M(),wt_);
          if (n_jets_ <= 1 && n_bjets_ > 0) {
            misc_2dplots_->Fill("m_gen_m_sv_btag",particles[i]->vector().M(),m_sv_,wt_);
          }
          if (n_bjets_ == 0) {
            misc_2dplots_->Fill("m_gen_m_sv_nobtag",particles[i]->vector().M(),m_sv_,wt_);
          }
          misc_2dplots_->Fill("m_gen_m_vis",particles[i]->vector().M(),m_vis_,wt_);
          dimu_vec = particles[i]->vector();
          met_vec += ditau_vec;
          met_vec -= dimu_vec;
          misc_2dplots_->Fill("orig_met_m_sv",met_vec.pt(),m_sv_,wt_);
          misc_2dplots_->Fill("orig_met_m_vis",met_vec.pt(),m_vis_,wt_);
        }
      }
    }

    if (write_tree_) outtree_->Fill();
    if (!write_plots_) return 0;

    // Define which selections this event passes
    if (channel_ == channel::et || channel_ == channel::etmet || channel_ == channel::mt || channel_ == channel::mtmet) {
      if (os_ && mt_1_ < 30.0) {
        SetPassSelection("os_sel");
        if (experimental_) {
          misc_2dplots_->Fill("jpt_vs_hpt", (jpt_1_ > 0.) ? jpt_1_ : 0., pt_h_, wt_);
          misc_2dplots_->Fill("hpt_vs_pt_tt", pt_h_, pt_tt_, wt_);
          if (boson) {
            misc_2dplots_->Fill("gen_vs_hpt", boson->pt(), pt_h_, wt_);
            misc_2dplots_->Fill("gen_vs_pt_tt", boson->pt(), pt_tt_, wt_);
          }
        }
      }
      if (os_) SetPassSelection("os");
      if (!os_) SetPassSelection("ss");
      if (os_ && mt_1_ > 70.0) SetPassSelection("os_con");
      if (os_ && mt_1_ > 60.0 && mt_1_ < 120.) SetPassSelection("os_con_mt_60-120");
      if (!os_ && mt_1_ < 30.0) SetPassSelection("ss_sel");
      if (!os_ && mt_1_ > 70.0) SetPassSelection("ss_con");
      if (!os_ && mt_1_ > 60.0 && mt_1_ < 120.) SetPassSelection("ss_con_mt_60-120");
    }

    if (channel_ == channel::em) {
      if (os_ && pzeta_ > -20) SetPassSelection("os_sel");
      if (os_) SetPassSelection("os");
      if (!os_) SetPassSelection("ss");
      if (!os_ && pzeta_ > -20) SetPassSelection("ss_sel");
    }

    // Define the 1- and 0-jet split based on pt_2
    double pt2_split = 40.0; // Tau pT for et,mt and mtmet
    if (channel_ == channel::em) pt2_split = 35.0;  // Mu pT for em


    // Inclusive Category
    SetPassCategory("inclusive");
    FillCoreControlPlots("inclusive");


    // auto lowpt_jets_copy = lowpt_jets;
    // ic::erase_if(lowpt_jets_copy,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
    // if (lowpt_jets_copy.size() >= 2 && n_bjets_ >= 1) SetPassCategory("sasha");
   
   //Some categories for MSSM, leave in for now
   if (n_jets_ <= 1 && prebjets.size() > 0) {
      SetPassCategory("prebtag");
      FillCoreControlPlots("prebtag");
    }

    if (n_jets_ <= 1 && n_bjets_ > 0) {
      SetPassCategory("btag");
      FillCoreControlPlots("btag");
    }
    if(n_prebjets_>1 && n_bjets_>0)
    {
        SetPassCategory("sasha");
    }
    if(n_prebjets_>1)
    {
        SetPassCategory("presasha");
        if(prebjetbcsv_1_ < 0.679 && prebjetbcsv_2_ < 0.679 ) {
            SetPassCategory("2jet0tag");
            FillCoreControlPlots("2jet0tag");
        } else if(prebjetbcsv_1_ > 0.679 && prebjetbcsv_2_ < 0.679 ) {
            SetPassCategory("2jet1tag");
            FillCoreControlPlots("2jet1tag");
        } else if(prebjetbcsv_1_ > 0.679 && prebjetbcsv_2_ > 0.679 ) {
            SetPassCategory("2jet2tag");
            FillCoreControlPlots("2jet2tag");
        }
    }
    if(n_prebjets_==1) {
        if(prebjetbcsv_1_ < 0.898) {
            SetPassCategory("1jet0tag");
            FillCoreControlPlots("1jet0tag");
        } else {
            SetPassCategory("1jet1tag");
            FillCoreControlPlots("1jet1tag");
        }
    }
    
    
    
    if (n_jets_ <= 1 && n_bjets_ > 0 && pt_2_ <= pt2_split) SetPassCategory("btag_low");
    if (n_jets_ <= 1 && n_bjets_ > 0 && pt_2_ > pt2_split)  SetPassCategory("btag_high");

    if (n_jets_ <= 1 && n_loose_bjets_ > 0)                       SetPassCategory("btag_loose");
    if (n_jets_ <= 1 && n_loose_bjets_ > 0 && pt_2_ <= pt2_split) SetPassCategory("btag_low_loose");
    if (n_jets_ <= 1 && n_loose_bjets_ > 0 && pt_2_ > pt2_split)  SetPassCategory("btag_high_loose");


    if (!PassesCategory("vbf") && n_bjets_ == 0) SetPassCategory("nobtag");

    return 0;
  }

  bool HhhCategories::PassesCategory(std::string const& category) const {
    std::map<std::string, bool>::const_iterator it = categories_.find(category);
    if (it != categories_.end()) {
      return it->second;
    } else {
      std::cerr << "Error in HhhCategories::PassesCategory: No category registered with label " << category << std::endl;
      throw;
      return false;
    }
  }


  void HhhCategories::InitMassPlots(std::string const& category) {
    for (std::map<std::string, bool>::const_iterator it = selections_.begin(); it != selections_.end(); ++it) {
      massplots_[category+"_"+it->first] = new MassPlots(fs_->mkdir(category+"_"+it->first));
    }

  }
  void HhhCategories::InitCoreControlPlots(std::string const& category) {
    for (std::map<std::string, bool>::const_iterator it = selections_.begin(); it != selections_.end(); ++it) {
      controlplots_[category+"_"+it->first] = new CoreControlPlots(fs_->mkdir(category+"_"+it->first));
    }
  }

  void HhhCategories::Reset() {
    for (std::map<std::string, bool>::iterator it = selections_.begin(); it != selections_.end(); ++it) {
      it->second = false;
    }
    for (std::map<std::string, bool>::iterator it = categories_.begin(); it != categories_.end(); ++it) {
      it->second = false;
    }
  }

  void HhhCategories::SetPassSelection(std::string const& selection) {
    std::map<std::string, bool>::iterator it = selections_.find(selection);
    if (it != selections_.end()) {
      it->second = true;
    } else {
      std::cerr << "Error in HhhCategories::SetPassSelection: No selection registered with label " << selection << std::endl;
      throw;
    }
  }

  void HhhCategories::SetPassCategory(std::string const& category) {
    std::map<std::string, bool>::iterator it = categories_.find(category);
    if (it != categories_.end()) {
      it->second = true;
      FillMassPlots(category);
      FillYields(category);
    } else {
      std::cerr << "Error in HhhCategories::SetPassCategory: No category registered with label " << category << std::endl;
      throw;
    }
  }

  void HhhCategories::FillMassPlots(std::string const& category) {
    for (std::map<std::string, bool>::iterator it = selections_.begin(); it != selections_.end(); ++it) {
      if (it->second) {
        std::map<std::string, MassPlots*>::iterator p_it = massplots_.find(category+"_"+it->first);
        if (p_it != massplots_.end()) {
          MassPlots *plots = p_it->second;
          plots->m_sv->Fill(m_sv_, wt_);
          plots->m_sv_sm->Fill(m_sv_, wt_);
          plots->m_sv_sm_fine->Fill(m_sv_, wt_);
          plots->m_sv_mssm->Fill(m_sv_, wt_);
          plots->m_sv_mssm_fine->Fill(m_sv_, wt_);
          plots->m_vis->Fill(m_vis_, wt_);
          plots->m_vis_sm->Fill(m_vis_, wt_);
          plots->m_vis_sm_fine->Fill(m_vis_, wt_);
          plots->m_vis_mssm->Fill(m_vis_, wt_);
          plots->m_vis_mssm_fine->Fill(m_vis_, wt_);
        }
      } 
    }
  }

  void HhhCategories::FillYields(std::string const& category) {
    for (std::map<std::string, bool>::iterator it = selections_.begin(); it != selections_.end(); ++it) {
      if (it->second) {
        yields_[category+"_"+it->first] = yields_[category+"_"+it->first] + wt_;
      } 
    }
  }


  void HhhCategories::FillCoreControlPlots(std::string const& category) {
    for (std::map<std::string, bool>::iterator it = selections_.begin(); it != selections_.end(); ++it) {
      if (it->second) {
        std::map<std::string, CoreControlPlots*>::iterator p_it = controlplots_.find(category+"_"+it->first);
        if (p_it != controlplots_.end()) {
          CoreControlPlots *plots = p_it->second;
          plots->n_vtx->Fill(n_vtx_, wt_);
          plots->mt_1->Fill(mt_1_, wt_);
          plots->pzeta->Fill(pzeta_, wt_);
          plots->pt_h->Fill(pt_h_, wt_);
          plots->pt_tt->Fill(pt_tt_, wt_);
          plots->pt_1->Fill(pt_1_, wt_);
          plots->pt_2->Fill(pt_2_, wt_);
          plots->eta_1->Fill(eta_1_, wt_);
          plots->eta_2->Fill(eta_2_, wt_);
          plots->z_2->Fill(z_2_, wt_);
          plots->m_2->Fill(m_2_, wt_);
          plots->met->Fill(met_, wt_);
          plots->l1_met->Fill(l1_met_, wt_);
          plots->calo_nohf_met->Fill(calo_nohf_met_, wt_);
          plots->met_phi->Fill(met_phi_, wt_);
          plots->n_jets->Fill(n_jets_, wt_);
          plots->n_bjets->Fill(n_bjets_, wt_);
          plots->n_jetsingap->Fill(n_jetsingap_, wt_);
          if (tau_decay_mode_ == 0) plots->tau_decay_mode->Fill(0. , wt_);
          if (tau_decay_mode_ == 1 || tau_decay_mode_ == 2) plots->tau_decay_mode->Fill(1., wt_);
          if (tau_decay_mode_ == 10) plots->tau_decay_mode->Fill(2., wt_);
          // Have to be careful here: we fill the jet control plots even
          // though these values may not be defined. This is to be sure
          // when plotting the inclusive category yield is consistent
          // if (n_jets_ >= 1) {
            plots->jpt_1->Fill(jpt_1_, wt_);
            plots->jeta_1->Fill(jeta_1_, wt_);
          // }
          // if (n_jets_ >= 2) {
            plots->jpt_2->Fill(jpt_2_, wt_);
            plots->jeta_2->Fill(jeta_2_, wt_);            
            plots->mjj->Fill(mjj_, wt_);
            plots->jdeta->Fill(jdeta_, wt_);
          // }
          // if (n_bjets_ >= 1) {
            plots->bpt_1->Fill(bpt_1_, wt_);
            plots->beta_1->Fill(beta_1_, wt_);
            plots->bcsv_1->Fill(bcsv_1_, wt_);
          // }
        }
      } 
    }
  }



  int HhhCategories::PostAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HhhCategories" << std::endl;
    std::cout << "-------------------------------------" << std::endl;      
    std::vector<std::string> print_selections;
    print_selections.push_back("os_sel");
    print_selections.push_back("os_con");
    print_selections.push_back("ss_sel");
    print_selections.push_back("ss_con");
    std::vector<std::string> print_cats;
    print_cats.push_back("inclusive");
    print_cats.push_back("vbf_loose");
    print_cats.push_back("vbf");
    // print_cats.push_back("vbf_tight");
    print_cats.push_back("1jet_high");
    print_cats.push_back("1jet_low");
    print_cats.push_back("0jet_high");
    print_cats.push_back("0jet_low");
    print_cats.push_back("btag");
    print_cats.push_back("nobtag");
    print_cats.push_back("sasha");
    print_cats.push_back("presasha");
    print_cats.push_back("1jet0tag");
    print_cats.push_back("1jet1tag");
    print_cats.push_back("2jet0tag");
    print_cats.push_back("2jet1tag");
    print_cats.push_back("2jet2tag");
    std::cout << boost::format("%-20s") % "Selections:";
    for (unsigned i = 0; i < print_selections.size(); ++i) {
      std::cout << boost::format("%-12s") % print_selections[i];
    }
    std::cout << std::endl;
    for (unsigned i = 0; i < print_cats.size(); ++i) {
      std::cout << boost::format("%-20s") % print_cats[i];
      for (unsigned j = 0; j < print_selections.size(); ++j) {
        std::cout << boost::format("%-12s") % (yields_[print_cats[i]+"_"+print_selections[j]]);
      }
      std::cout << std::endl;
    }

    return 0;
  }

  void HhhCategories::PrintInfo() {
    ;
  }
}
