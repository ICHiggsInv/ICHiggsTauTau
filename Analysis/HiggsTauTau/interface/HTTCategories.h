#ifndef ICHiggsTauTau_HiggsTauTau_HTTCategories_h
#define ICHiggsTauTau_HiggsTauTau_HTTCategories_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HTTPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"

#include <string>

namespace ic {

class HTTCategories : public ModuleBase {

 private:
  std::string jets_label_;
  CLASS_MEMBER(HTTCategories, std::string, ditau_label)
  CLASS_MEMBER(HTTCategories, std::string, met_label)
  CLASS_MEMBER(HTTCategories, double, mass_shift)
  CLASS_MEMBER(HTTCategories, ic::channel, channel)
  CLASS_MEMBER(HTTCategories, ic::era, era)
  CLASS_MEMBER(HTTCategories, ic::strategy, strategy)
  CLASS_MEMBER(HTTCategories, bool, write_tree)
  CLASS_MEMBER(HTTCategories, bool, write_plots)
  CLASS_MEMBER(HTTCategories, bool, experimental)
  CLASS_MEMBER(HTTCategories, fwlite::TFileService*, fs)

  std::map<std::string, bool> categories_;
  std::map<std::string, bool> selections_;
  std::map<std::string, MassPlots*> massplots_;
  std::map<std::string, double> yields_;
  std::map<std::string, CoreControlPlots*> controlplots_;
  DynamicHistoSet * misc_plots_;
  Dynamic2DHistoSet * misc_2dplots_;
  
  TTree *outtree_;

  // Event Properties
  double wt_;
  double wt_ggh_pt_up_;
  double wt_ggh_pt_down_;
  double wt_tau_fake_up_;
  double wt_tau_fake_down_;
  double wt_tquark_up_;
  double wt_tquark_down_;
  double wt_tau_id_up_;
  double wt_tau_id_down_;
  bool os_;
  unsigned n_vtx_;
  double m_sv_;
  double m_vis_;
  double pt_h_;
  double pt_tt_;
  double mt_1_;
  double mt_ll_;
  double pzeta_;
  double pzetavis_;
  double pzetamiss_;
  double emu_dphi_;
  double emu_csv_;
  double emu_dxy_1_;
  double emu_dxy_2_;
  double pt_1_;
  double pt_2_;
  double eta_1_;
  double eta_2_;
  double iso_1_;
  double iso_2_;
  double z_2_;
  double m_2_;
  double met_;
  double met_phi_;

  int    tau_decay_mode_;

  unsigned n_jets_;
  unsigned n_lowpt_jets_;
  unsigned n_bjets_;
  unsigned n_loose_bjets_;
  unsigned n_jetsingap_; // Defined if n_jets >= 2
  double jpt_1_;     // Defined if n_jets >= 1
  double jpt_2_;     // Defined if n_jets >= 2
  double jeta_1_;    // Defined if n_jets >= 1
  double jeta_2_;    // Defined if n_jets >= 2
  double bpt_1_;     // Defined if n_bjets >= 1
  double beta_1_;    // Defined if n_bjets >= 1
  double bcsv_1_; 

  int j1_dm_;

  double mjj_;       // Defined if n_jets >= 2
  double jdeta_;     // Defined if n_jets >= 2

  double mjj_lowpt_;       // Defined if n_lowpt_jets >= 2
  double jdeta_lowpt_;     // Defined if n_lowpt_jets >= 2
  unsigned n_jetsingap_lowpt_; // Defined if n_lowpt_jets >= 2

  unsigned n_prebjets_;

  double l1_met_;
  double calo_nohf_met_;

  double em_gf_mva_;
  double em_vbf_mva_;

    // Other VBF MVA variables?

 public:
  HTTCategories(std::string const& name);
  virtual ~HTTCategories();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

  bool PassesCategory(std::string const& category) const;

  void InitSelection(std::string const& selection);
  void InitCategory(std::string const& category);
  void InitMassPlots(std::string const& category);
  void FillMassPlots(std::string const& category);
  void FillYields(std::string const& category);
  void InitCoreControlPlots(std::string const& category);
  void FillCoreControlPlots(std::string const& category);

  void SetPassCategory(std::string const& category);
  void SetPassSelection(std::string const& selection);

  void Reset();


};

}


#endif
