#ifndef ICHiggsTauTau_Module_LightTree_h
#define ICHiggsTauTau_Module_LightTree_h

#include <string>

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"

#include "TTree.h"

namespace ic {
  class LightTree : public ModuleBase {

  private:

    CLASS_MEMBER(LightTree,fwlite::TFileService*, fs);
        CLASS_MEMBER(LightTree,std::string, met_label);
    CLASS_MEMBER(LightTree,std::string, dijet_label);
    CLASS_MEMBER(LightTree,std::string, sel_label);
    CLASS_MEMBER(LightTree,bool, is_data);
    CLASS_MEMBER(LightTree,bool, ignoreLeptons);
    CLASS_MEMBER(LightTree,bool, dotrigskim);
    CLASS_MEMBER(LightTree,bool, do_promptskim);
    CLASS_MEMBER(LightTree,bool, is_embedded);
    CLASS_MEMBER(LightTree,std::string, trig_obj_label);
    CLASS_MEMBER(LightTree,std::string, trigger_path);

    TTree *outputTree_;
    
    unsigned run_;
    unsigned lumi_;
    unsigned event_;
    double weight_nolep_;
    double total_weight_lepveto_;
    double total_weight_leptight_;
    double puweight_up_scale_;
    double puweight_down_scale_;
    double topweight_up_scale_;
    double topweight_down_scale_;
    double jet1_pt_;
    double jet2_pt_;
    double jet3_pt_;
    double jet1_E_;
    double jet2_E_;
    double jet3_E_;
    double jet1_eta_;
    double jet2_eta_;
    double jet3_eta_;
    double forward_tag_eta_;
    double central_tag_eta_;
    double jet1_phi_;
    double jet2_phi_;
    double jet3_phi_;
    double jet_csv1_;
    double jet_csv2_;
    double jet_csv3_;
    double dijet_M_;
    double dijet_deta_;
    double dijet_sumeta_;
    double dijet_dphi_;
    double met_;
    double met_x_;
    double met_y_;
    double metnomu_x_;
    double metnomu_y_;
    double met_significance_;
    double metnomu_significance_;
    double sumet_;
    double ht_;
    double ht30_;
    double mht_;
    double sqrt_ht_;
    double unclustered_et_;
    double jet1met_dphi_;
    double jet2met_dphi_;
    double jet1metnomu_dphi_;
    double jet2metnomu_dphi_;
    double jetmet_mindphi_;
    double jetmetnomu_mindphi_;
    double alljetsmet_mindphi_;
    double alljetsmetnomu_mindphi_;
    double jetunclet_mindphi_;
    double metunclet_dphi_;
    double metnomuunclet_dphi_;
    double dijetmet_scalarSum_pt_;
    double dijetmet_vectorialSum_pt_;
    double dijetmet_ptfraction_;
    double jet1met_scalarprod_;
    double jet2met_scalarprod_;
    double dijetmetnomu_scalarSum_pt_;
    double dijetmetnomu_vectorialSum_pt_;
    double dijetmetnomu_ptfraction_;
    double jet1metnomu_scalarprod_;
    double jet2metnomu_scalarprod_;
    double cjvjetpt_;
    unsigned n_jets_cjv_30_;
    unsigned n_jets_cjv_20EB_30EE_;
    unsigned n_jets_15_;
    unsigned n_jets_30_;
    double passtrigger_;
    double passparkedtrigger1_;
    double passparkedtrigger2_;
    double l1met_;
    double metnomuons_;

    int nvetomuons_;
    int nselmuons_;
    int nvetoelectrons_;
    int nselelectrons_;
    int ntaus_;
    double m_mumu_;
    double m_ee_;
    double m_mumu_gen_;
    double m_ee_gen_;
    double genlep1_pt_;
    double genlep1_eta_;
    double genlep1_phi_;
    double genlep1_id_;
    double genlep2_pt_;
    double genlep2_eta_;
    double genlep2_phi_;
    double genlep2_id_;
    double mu1_pt_;
    double mu1_eta_;
    double mu1_phi_;
    double mu2_pt_;
    double mu2_eta_;
    double mu2_phi_;
    double ele1_pt_;
    double ele1_eta_;
    double ele1_phi_;
    double tau1_pt_;
    double tau1_eta_;
    double tau1_phi_;
    double lep_mt_;
    int n_vertices_;

  public:
    LightTree(std::string const& name);
    virtual ~LightTree();
    
    virtual int PreAnalysis();
    virtual int Execute(TreeEvent *event);
    virtual int PostAnalysis();
    virtual void PrintInfo();

   };

}


#endif
