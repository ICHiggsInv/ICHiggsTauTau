#include <string>
// #include "boost/filesystem.hpp"
// #include "boost/regex.hpp"
// #include "boost/format.hpp"
#include "boost/bind.hpp"
// #include "boost/assign/list_of.hpp"
#include "CombineTools/interface/CombineHarvester.h"
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/HelperFunctions.h"
#include "CombineTools/interface/HttSystematics.h"
#include "CombineTools/interface/RooHttSMPdf.h"
#include "CombineTools/interface/RooHttYield.h"
// #include "TH1F.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

using boost::bind;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;

int main() {
  TH1::AddDirectory(false);
  ch::CombineHarvester cb;

  vector<pair<int, string>> mt_cats_7 = {
      make_pair(1, "muTau_0jet_medium"),
      make_pair(2, "muTau_0jet_high"),
      make_pair(3, "muTau_1jet_medium"),
      make_pair(4, "muTau_1jet_high_lowhiggs"),
      make_pair(5, "muTau_1jet_high_mediumhiggs"),
      make_pair(6, "muTau_vbf")};

  vector<pair<int, string>> mt_cats_8 = {
      make_pair(1, "muTau_0jet_medium"),
      make_pair(2, "muTau_0jet_high"),
      make_pair(3, "muTau_1jet_medium"),
      make_pair(4, "muTau_1jet_high_lowhiggs"),
      make_pair(5, "muTau_1jet_high_mediumhiggs"),
      make_pair(6, "muTau_vbf_loose"),
      make_pair(7, "muTau_vbf_tight")};

  vector<string> masses = {"90",  "95",  "100", "105", "110", "115",
                           "120", "125", "130", "135", "140", "145"};
  // vector<string> masses = {"125"};

  std::cout << "Adding observations...";
  cb.AddObservations({"*"}, {"htt"}, {"7TeV"}, {"mt"}, mt_cats_7);
  cb.AddObservations({"*"}, {"htt"}, {"8TeV"}, {"mt"}, mt_cats_8);
  std::cout << " done\n";

  std::cout << "Adding background processes...";
  cb.AddProcesses({"*"}, {"htt"}, {"7TeV"}, {"mt"},
      {"ZTT", "W", "QCD", "ZL", "ZJ", "TT", "VV"}, mt_cats_7, false);
  cb.AddProcesses({"*"}, {"htt"}, {"8TeV"}, {"mt"},
      {"ZTT", "W", "QCD", "ZL", "ZJ", "TT", "VV"}, mt_cats_8, false);
  std::cout << " done\n";

  std::cout << "Adding signal processes...";
  cb.AddProcesses(masses, {"htt"}, {"7TeV"}, {"mt"},
      {"ggH", "qqH", "WH", "ZH"}, mt_cats_7, true);
  cb.AddProcesses(masses, {"htt"}, {"8TeV"}, {"mt"},
      {"ggH", "qqH", "WH", "ZH"}, mt_cats_8, true);
  std::cout << " done\n";

  /*
  (1) hww(mH) as signal:
  cb.AddProcesses(masses, {"htt"}, {"7TeV"}, {"mt"},
      {"ggH_hww", "qqH_hww"}, mt_cats_7, true);
  cb.AddProcesses(masses, {"htt"}, {"8TeV"}, {"mt"},
      {"ggH_hww", "qqH_hww"}, mt_cats_8, true);

  (2) hww(125) as bkg:
  cb.AddProcesses({"*"}, {"htt"}, {"7TeV"}, {"mt"},
      {"ggH_hww125", "qqH_hww125"}, mt_cats_7, false);
  cb.AddProcesses({"*"}, {"htt"}, {"8TeV"}, {"mt"},
      {"ggH_hww125", "qqH_hww125"}, mt_cats_8, false);

  (2) hww(mH) as bkg:
  for (m : masses) {
    cb.AddProcesses({m}, {"htt"}, {"7TeV"}, {"mt"},
        {"ggH_hww"+m, "qqH_hww"+m}, mt_cats_7, false);
    cb.AddProcesses({m}, {"htt"}, {"8TeV"}, {"mt"},
        {"ggH_hww"+m, "qqH_hww"+m}, mt_cats_8, false);
  }
  */

  std::cout << "Adding systematic uncertainties...";
  ch::AddDefaultSystematics(&cb);
  std::cout << " done\n";

  std::cout << "Extracting histograms from input root files...";
  cb.cp().era({"7TeV"}).backgrounds().ExtractShapes(
      "data/demo/htt_mt.inputs-sm-7TeV-hcg.root", "$CHANNEL/$PROCESS",
      "$CHANNEL/$PROCESS_$SYSTEMATIC");
  cb.cp().era({"8TeV"}).backgrounds().ExtractShapes(
      "data/demo/htt_mt.inputs-sm-8TeV-hcg.root", "$CHANNEL/$PROCESS",
      "$CHANNEL/$PROCESS_$SYSTEMATIC");
  cb.cp().era({"7TeV"}).signals().ExtractShapes(
      "data/demo/htt_mt.inputs-sm-7TeV-hcg.root", "$CHANNEL/$PROCESS$MASS_fine_binning",
      "$CHANNEL/$PROCESS$MASS_$SYSTEMATIC_fine_binning");
  cb.cp().era({"8TeV"}).signals().ExtractShapes(
      "data/demo/htt_mt.inputs-sm-8TeV-hcg.root", "$CHANNEL/$PROCESS$MASS_fine_binning",
      "$CHANNEL/$PROCESS$MASS_$SYSTEMATIC_fine_binning");
  std::cout << " done\n";

  std::cout << "Scaling signal process rates...\n";
  map<string, TGraph> xs;
  for (std::string const& e : {"7TeV", "8TeV"}) {
    for (std::string const& p : {"ggH", "qqH", "WH", "ZH"}) {
      ch::ParseTable(&xs, "data/xsecs_brs/"+p+"_"+e+"_YR3.txt", {p+"_"+e});
    }
  }
  ch::ParseTable(&xs, "data/xsecs_brs/htt_YR3.txt", {"htt"});

  for (std::string const& e : {"7TeV", "8TeV"}) {
    for (std::string const& p : {"ggH", "qqH", "WH", "ZH"}) {
      std::cout << "Scaling for process " << p << " and era " << e << "\n";

      cb.cp().process({p}).era({e}).ForEachProc(
          bind(ch::ScaleProcessRate, _1, &xs, p+"_"+e, "htt"));
    }
  }

  std::cout << "Setting standardised bin names...";
  cb.ForEachObs(ch::SetStandardBinName<ch::Observation>);
  cb.ForEachProc(ch::SetStandardBinName<ch::Process>);
  cb.ForEachNus(ch::SetStandardBinName<ch::Nuisance>);
  std::cout << " done\n";


  // cb.bin({"htt_mt_5_8TeV", "htt_mt_6_8TeV"});
  // cb.nus_type({"shape"}, false);
  set<string> bins =
      cb.GenerateSetFromObs<string>(mem_fn(&ch::Observation::bin));
  RooWorkspace ws("htt", "htt");
  RooRealVar mtt("CMS_th1x", "CMS_th1x", 0, 350);
  RooRealVar mh("MH", "MH", 125, 90, 145);
  for (auto b : bins) {
    ch::CombineHarvester tmp = std::move(cb.cp().bin({b}).signals());
    set<string> sigs =
        tmp.GenerateSetFromProcs<string>(mem_fn(&ch::Process::process));
    for (auto s : sigs) {
      ch::CombineHarvester tmp2 = std::move(tmp.cp().process({s}));
      TH1F data_hist = tmp2.GetObservedShape();
      TH1F shape("tmp", "tmp", data_hist.GetNbinsX(), 0.,
                 static_cast<float>(data_hist.GetNbinsX()));
      RooDataHist dh((b + "_" + s + "_hist").c_str(),
                     (b + "_" + s + "_hist").c_str(), RooArgList(mtt),
                     RooFit::Import(shape, false));
      RooHttSMPdf pdf((b + "_" + s + "_pdf").c_str(),
                      (b + "_" + s + "_pdf").c_str(), RooArgList(mtt), dh, mh);
      RooHttYield yield((b + "_" + s + "_pdf_norm").c_str(),
                        (b + "_" + s + "_pdf_norm").c_str(), mh);
      tmp2.ForEachProc([&](ch::Process* p) {
        double m = boost::lexical_cast<double>(p->mass());
        pdf.AddPoint(m, *((TH1F*)p->shape()));
        yield.AddPoint(m, p->rate());
        p->set_shape(nullptr);
        p->set_rate(1.0);
      });
      unsigned nbins = data_hist.GetXaxis()->GetXbins()->GetSize();
      double const* arr = data_hist.GetXaxis()->GetXbins()->GetArray();
      std::vector<double> binning(nbins);
      binning.assign(arr, arr+nbins);
      pdf.SetBinning(binning);
      ws.import(pdf);
      ws.import(yield);
    }
  }
  cb.AddWorkspace(&ws);
  ws.Print();
  cb.cp().signals().ExtractPdfs("htt", "$CHANNEL_$PROCESS_pdf");

  // cb.cp().mass({"125", "*"}).PrintAll();
  std::cout << "Merging bin errors...\n";
  cb.cp().bin_id({0, 1, 2, 3, 4}).process({"W", "QCD"})
      .MergeBinErrors(0.1, 0.4);
  cb.cp().bin_id({5}).era({"7TeV"}).process({"W"})
      .MergeBinErrors(0.1, 0.4);
  cb.cp().bin_id({5, 6}).era({"8TeV"}).process({"W"})
      .MergeBinErrors(0.1, 0.4);
  cb.cp().bin_id({7}).era({"8TeV"}).process({"W", "ZTT"})
      .MergeBinErrors(0.1, 0.4);
  std::cout << "...done\n";

  std::cout << "Generating bbb uncertainties...\n";
  cb.cp().bin_id({0, 1, 2, 3, 4}).process({"W", "QCD"})
      .AddBinByBin(0.1, true, &cb);
  cb.cp().bin_id({5}).era({"7TeV"}).process({"W"})
      .AddBinByBin(0.1, true, &cb);
  cb.cp().bin_id({5, 6}).era({"8TeV"}).process({"W"})
      .AddBinByBin(0.1, true, &cb);
  cb.cp().bin_id({7}).era({"8TeV"}).process({"W", "ZTT"})
      .AddBinByBin(0.1, true, &cb);
  std::cout << "...done\n";

  TFile output("htt_mt.input.root", "RECREATE");
  cb.cp().mass({"125", "*"}).WriteDatacard("htt_mt_all.txt", output);
  // for (auto b : bins) {
  //   for (auto m : /*masses*/ {"125"}) {
  //     cout << "Writing datacard for bin: " << b << " and mass: " << m
  //               << "\n";
  //     cb.cp().bin({b}).mass({m, "*"}).WriteDatacard(
  //         b + "_" + m + ".txt", output);
  //   }
  // }

  /* Still to do:
     [x] adding systematics
     [ ] adding workspaces, data and pdfs
     [x] scaling signal processes
     [x] merging bin errors
     [x] adding bbb uncertainties
     [ ] interpolated signal points (with fb templates)
     [x] template function to assign common properties between objects
  */
}
