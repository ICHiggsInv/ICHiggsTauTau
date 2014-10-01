#ifndef ICHiggsTauTau_Module_ModifyMet_h
#define ICHiggsTauTau_Module_ModifyMet_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/interface/Candidate.hh"

#include <string>
#include "TMath.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"

namespace ic {
  
  class ModifyMet : public ModuleBase {
  private:
    std::string met_name_;
    std::string lep_name_;
    unsigned lepFlavour_;
    unsigned nLepToAdd_;
    double mindphi_;

    template<class T> void correctMet(TreeEvent *event, Candidate::Vector & aVec){
     
      std::vector<T*> lLeps = event->GetPtrVec<T>(lep_name_);
      std::sort(lLeps.begin(), lLeps.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

      ROOT::Math::PtEtaPhiEVector unmodified = aVec;

      unsigned nLep = (nLepToAdd_ > lLeps.size()) ? lLeps.size() : nLepToAdd_;
      for (unsigned iL(0); iL<nLep;++iL){
	bool addLep = true;
	//filter jets
	if (lepFlavour_==4){
	  double dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(unmodified,lLeps[iL]->vector()));
	  if (dphi>mindphi_ || iL==0 || iL==1) addLep=false;
	}

	if (addLep) aVec += lLeps[iL]->vector();
      }
    };

  public:
    ModifyMet(std::string const& name, std::string met_name, std::string lep_name, 
	      unsigned lepFlavour, unsigned nLepToAdd,
	      double mindphi=1.0);
    virtual ~ModifyMet();



    virtual int PreAnalysis();
    virtual int Execute(TreeEvent *event);
    virtual int PostAnalysis();
    virtual void PrintInfo();
    
  };

}//namespace


#endif
