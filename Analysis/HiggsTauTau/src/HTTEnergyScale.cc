#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTEnergyScale.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
#include "boost/format.hpp"

namespace ic {

  HTTEnergyScale::HTTEnergyScale(std::string const& name) : ModuleBase(name), channel_(channel::et), strategy_(strategy::moriond2013) {
    moriond_corrections_ = false;
  }

  HTTEnergyScale::~HTTEnergyScale() {
    ;
  }

  int HTTEnergyScale::PreAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HTTEnergyScale" << std::endl;
    std::cout << "-------------------------------------" << std::endl;    
    std::cout << boost::format(param_fmt()) % "moriond_corrections" % moriond_corrections_;
    std::cout << boost::format(param_fmt()) % "input_label"         % input_label_;
    std::cout << boost::format(param_fmt()) % "shift"               % shift_;
    return 0;
  }

  int HTTEnergyScale::Execute(TreeEvent *event) {
    std::vector<Tau *> & taus = event->GetPtrVec<Tau>(input_label_);
    std::map<std::size_t, double> tau_scales;
    for (unsigned i = 0; i < taus.size(); ++i) {
      double central_shift = 1.00;
      if (moriond_corrections_) {
        int dm = taus[i]->decay_mode();
        double pt = taus[i]->pt();
        if (dm == 0) {
          central_shift = 1.00;
        } else if (dm == 1 || dm == 2) {
          if (strategy_ == strategy::paper2013) {
            central_shift = 1.012;
          } else {
            central_shift = 1.015 + 0.001 * TMath::Min(TMath::Max(pt-45. ,0.),10.0);
          }
        } else if (dm == 10) {
          if (strategy_ == strategy::paper2013) {
            central_shift = 1.012;
          } else {
            central_shift = 1.012 + 0.001 * TMath::Min(TMath::Max(pt-32. ,0.),18.0);
          }
        }
      } else {
        central_shift = 1.00;
      }
      double total_shift = central_shift * shift_;
      taus[i]->set_pt(taus[i]->pt() * total_shift);
      taus[i]->set_energy(taus[i]->energy() * total_shift);
      tau_scales[taus[i]->id()] = total_shift;
    }
    event->Add("tau_scales", tau_scales);
    return 0;
  }



  int HTTEnergyScale::PostAnalysis() {
    return 0;
  }

  void HTTEnergyScale::PrintInfo() {
    ;
  }

}
