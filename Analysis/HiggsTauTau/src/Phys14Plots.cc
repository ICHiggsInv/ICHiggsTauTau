#include "HiggsTauTau/interface/Phys14Plots.h"
#include "Utilities/interface/FnPredicates.h"
#include "Utilities/interface/FnPairs.h"
#include "HiggsTauTau/interface/HTTGenEvent.h"

namespace ic {

Phys14Plots::Phys14Plots(std::string const& name) : ModuleBase(name) {}

Phys14Plots::~Phys14Plots() { ; }

int Phys14Plots::PreAnalysis() {
  using ROOT::Math::Pi;
  if (!fs_) return 0;
  dir_ = new TFileDirectory(fs_->mkdir("phys14"));
  h_n_vtx = dir_->make<TH1F>("n_vtx", "",71, -0.5, 70.5);

  h_gen_h_pt = dir_->make<TH1F>("gen_h_pt", "",     250, 0, 500);
  h_gen_h_eta = dir_->make<TH1F>("gen_h_eta", "",   50, -5, 5);
  h_gen_h_phi = dir_->make<TH1F>("gen_h_phi", "",   50, -Pi(), Pi());
  h_gen_h_mass = dir_->make<TH1F>("gen_h_mass", "", 100, 0, 200);
  // = dir_->make<TH1F>("", "", , , );
  h_gen_th_pt = dir_->make<TH1F>("gen_th_pt", "", 100, 0, 200);
  h_gen_th_eta = dir_->make<TH1F>("gen_th_eta", "", 50, -5, 5);
  h_gen_th_mode = dir_->make<TH1F>("gen_th_mode", "", 20, -0.5, 19.5);


  th_dm_eff_vs_pt = EfficiencyPlot1D(dir_, "th_dm_eff_vs_pt", 100, 0, 200);
  th_dm_eff_vs_eta = EfficiencyPlot1D(dir_, "th_dm_eff_vs_eta", 50, -5, 5);
  th_dm_eff_vs_nvtx = EfficiencyPlot1D(dir_, "th_dm_eff_vs_nvtx", 51, -0.5, 50.5);

  th0_dm_eff_vs_pt = EfficiencyPlot1D(dir_, "th0_dm_eff_vs_pt", 100, 0, 200);
  th0_dm_eff_vs_eta = EfficiencyPlot1D(dir_, "th0_dm_eff_vs_eta", 50, -5, 5);
  th0_dm_eff_vs_nvtx = EfficiencyPlot1D(dir_, "th0_dm_eff_vs_nvtx", 51, -0.5, 50.5);

  th1_dm_eff_vs_pt = EfficiencyPlot1D(dir_, "th1_dm_eff_vs_pt", 100, 0, 200);
  th1_dm_eff_vs_eta = EfficiencyPlot1D(dir_, "th1_dm_eff_vs_eta", 50, -5, 5);
  th1_dm_eff_vs_nvtx = EfficiencyPlot1D(dir_, "th1_dm_eff_vs_nvtx", 51, -0.5, 50.5);

  th10_dm_eff_vs_pt = EfficiencyPlot1D(dir_, "th10_dm_eff_vs_pt", 100, 0, 200);
  th10_dm_eff_vs_eta = EfficiencyPlot1D(dir_, "th10_dm_eff_vs_eta", 50, -5, 5);
  th10_dm_eff_vs_nvtx = EfficiencyPlot1D(dir_, "th10_dm_eff_vs_nvtx", 51, -0.5, 50.5);

  h_th_mode_table =
      dir_->make<TH2F>("th_mode_table", "", 17, -1.5, 15.5, 17, -1.5, 15.5);

  // Want to plot:
  // Probability to pass decay mode reco vs: gen vis pT, gen vis eta, nvtx
  return 0;
}

int Phys14Plots::Execute(TreeEvent* event) {
  auto const* event_info = event->GetPtr<EventInfo>("eventInfo");
  unsigned n_vtx = event_info->good_vertices();
  h_n_vtx->Fill(n_vtx);

  auto const& gevt = event->Get<GenEvent_XToTauTau>("genEvent_XToTauTau");
  if (gevt.boson) {
    h_gen_h_pt->Fill(gevt.boson->pt());
    h_gen_h_eta->Fill(gevt.boson->eta());
    h_gen_h_phi->Fill(gevt.boson->phi());
    h_gen_h_mass->Fill(gevt.boson->M());
  }

  // A vector of the GenEvent infos for any hadronic tau decays
  std::vector<GenEvent_Tau const*> gen_th_vec;
  // A vector of the visible hadronic tau jets (only really needed as input to
  // the matching algorithm)
  std::vector<GenJet const*> gen_th_vis_vec;
  if (gevt.tau_0.hadronic_mode >= 0) gen_th_vec.push_back(&(gevt.tau_0));
  if (gevt.tau_1.hadronic_mode >= 0) gen_th_vec.push_back(&(gevt.tau_1));
  for (auto th : gen_th_vec) gen_th_vis_vec.push_back(th->vis_jet);

  // Get the reco tau collection
  auto reco_th_vec = event->GetPtrVec<Tau>("taus");

  // Do the matching on the unbiased gen. and reco. taus

  auto gen_rec_th_matches =
      MatchByDR(gen_th_vis_vec, reco_th_vec, 0.3, true, true);


  // We have a vector of pairs from MatchByDR, but a more useful
  // format will be a map going from the gen. to the reco. tau
  std::map<GenJet const*, Tau const*> gen_rec_th_map;
  for (auto const& x : gen_rec_th_matches) {
    gen_rec_th_map[x.first] = x.second;
  }

  // Now loop through the hadronic taus
  for (auto gen_th : gen_th_vec) {
    // Get the visible jet
    GenJet const* gen_th_vis = gen_th->vis_jet;
    h_gen_th_pt->Fill(gen_th_vis->pt());
    h_gen_th_eta->Fill(gen_th_vis->eta());
    h_gen_th_mode->Fill(gen_th->hadronic_mode);

    // Use decay mode -1 to mean not reconstructed
    int reco_mode = -1;
    Tau const* rec_th = nullptr;
    if (gen_rec_th_map.count(gen_th_vis)) {
      rec_th = gen_rec_th_map[gen_th_vis];
      if (rec_th->decay_mode() >= 0) reco_mode = rec_th->decay_mode();
    }

    // The mode table definitely seems to require the reco acceptance cuts be
    // applied. Not clear if they should also be applied at the gen level. For
    // now let's just apply it at the reco level
    if (rec_th && MinPtMaxEta(rec_th, th_pt_acc, th_eta_acc)) {
      h_th_mode_table->Fill(gen_th->hadronic_mode, reco_mode);
    }



    // Now we'll do the efficiencies. Here the definition is clearer: gen
    // acceptance is in the denominator
    if (MinPtMaxEta(gen_th_vis, th_pt_acc, th_eta_acc)) {
      // We're in the denominator
      bool pass_dm = false;
      if (rec_th && rec_th->GetTauID("decayModeFinding") > 0.5 &&
          MinPtMaxEta(rec_th, th_pt_acc, th_eta_acc))
        pass_dm = true;
      th_dm_eff_vs_pt.Fill(gen_th_vis->pt(), pass_dm);
      th_dm_eff_vs_eta.Fill(gen_th_vis->eta(), pass_dm);
      th_dm_eff_vs_nvtx.Fill(n_vtx, pass_dm);
      if (gen_th->hadronic_mode == 0) {
        th0_dm_eff_vs_pt.Fill(gen_th_vis->pt(), pass_dm);
        th0_dm_eff_vs_eta.Fill(gen_th_vis->eta(), pass_dm);
        th0_dm_eff_vs_nvtx.Fill(n_vtx, pass_dm);
      }
      if (gen_th->hadronic_mode >= 1 && gen_th->hadronic_mode <= 4) {
        th1_dm_eff_vs_pt.Fill(gen_th_vis->pt(), pass_dm);
        th1_dm_eff_vs_eta.Fill(gen_th_vis->eta(), pass_dm);
        th1_dm_eff_vs_nvtx.Fill(n_vtx, pass_dm);
      }
      if (gen_th->hadronic_mode == 10) {
        th10_dm_eff_vs_pt.Fill(gen_th_vis->pt(), pass_dm);
        th10_dm_eff_vs_eta.Fill(gen_th_vis->eta(), pass_dm);
        th10_dm_eff_vs_nvtx.Fill(n_vtx, pass_dm);
      }
    }
  }
  return 0;
}

int Phys14Plots::PostAnalysis() { return 0; }

void Phys14Plots::PrintInfo() { ; }
}