#ifndef ICHiggsTauTau_Zbb_ZbbAnalyser_h
#define ICHiggsTauTau_Zbb_ZbbAnalyser_h
#include <string>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Core/interface/TreeEvent.h"
#include "Core/interface/ModuleBase.h"
#include "Utilities/interface/BTagWeight.h"
namespace ic { class SecondaryVertex; }
namespace ic { class Track; }


namespace ic {

class ZbbAnalyser : public ModuleBase {
 private:
  CLASS_MEMBER(ZbbAnalyser, fwlite::TFileService*, fs)
  CLASS_MEMBER(ZbbAnalyser, bool, set_z_flav)
  TTree *outtree_;
  BTagWeight btag_weight_;

  // Event Properties
  bool is_ee_;
  unsigned zflav_; // 0 = Z+l, 1 = Z+c, 2=Z+b;
  double wt_;
  double wt_1b_inc_;
  double wt_1b_exc_;
  double wt_2b_inc_;

  double wt_1b_inc_d_;
  double wt_1b_exc_d_;
  double wt_2b_inc_d_;

  double wt_1b_inc_u_;
  double wt_1b_exc_u_;
  double wt_2b_inc_u_;

  // global variables
  double met_;
  double m_z_;

  // Leading lepton
  double pt_1_;
  double eta_1_;

  // Subleading lepton
  double pt_2_;
  double eta_2_;

  unsigned n_jets_;
  double jpt_1_;
  double jeta_1_;
  double jnsv_1_;
  double jssv_1_;
  double jpt_2_;
  double jeta_2_;
  double jnsv_2_;
  double jssv_2_;

  unsigned n_b_jets_;

  // Highest pT b-tagged jet
  double bpt_1_;
  double beta_1_;
  double m_sv_1_;

  double bpt_2_;
  double beta_2_;
  double m_sv_2_;

  double dphi_z_bb_;
  double pt_bb_;
  double pt_z_;

  double dphi_z_b_;

  double ElectronWeight(Candidate const* elec1, Candidate const* elec2);
  double MuonWeight(Candidate const* muon1, Candidate const* muon2);
  double SVMass(Jet const* jet,
      std::map<std::size_t, SecondaryVertex *> const& sv_map,
      std::map<std::size_t, Track *> const& trk_map);

 public:
  ZbbAnalyser(std::string const& name);
  virtual ~ZbbAnalyser();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

};
}


#endif
