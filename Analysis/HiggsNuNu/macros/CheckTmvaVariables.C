#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>

#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h" 
#include "TGraphErrors.h" 
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"



TH1F * make_histogram(std::string aName, std::string aAxisTitle, std::string aUnit,
		      int nBins, float xMin, float xMax){
  std::string lUnit = aUnit;
  if (lUnit.size() > 0) lUnit = " ("+lUnit+")";
  TH1F *hTmp = new TH1F(aName.c_str(),(";"+aAxisTitle+lUnit).c_str(),nBins,xMin,xMax);
  return hTmp;
};

void MakeAllHistograms(TH1F * hist[], const std::string aLabel){
  hist[0] = make_histogram((aLabel+"_jet1_pt").c_str(),"Jet 1 p_{T}", "GeV", 100,0,1000);
  hist[1] = make_histogram((aLabel+"_jet2_pt").c_str(),"Jet 2 p_{T}", "GeV", 100,0,1000);
  hist[2] = make_histogram((aLabel+"_jet1_eta").c_str(),"Jet 1 #eta", "", 50,-5.,5);
  hist[3] = make_histogram((aLabel+"_jet2_eta").c_str(),"Jet 2 #eta", "", 50,-5,5);
  hist[4] = make_histogram((aLabel+"_dijet_M").c_str(),"M_{jj}", " GeV", 100,0,5000);
  hist[5] = make_histogram((aLabel+"_dijet_deta").c_str(),"#Delta#eta_{jj}", "", 100,0,8);
  hist[6] = make_histogram((aLabel+"_dijet_sumeta").c_str(),"#eta_{j1}+#eta_{j2}", "", 100,-5,5);
  hist[7] = make_histogram((aLabel+"_dijet_dphi").c_str(),"#Delta#phi_{jj}", "", 100,0,3.1416);
  hist[8] = make_histogram((aLabel+"_met").c_str(),"MET", "GeV", 100,0,1000);
  hist[9] = make_histogram((aLabel+"_met_phi").c_str(),"MET #phi", "", 100,-3.1416,3.1416);
  hist[10] = make_histogram((aLabel+"_met_significance").c_str(),"MET significance", "", 100,0,20);
  hist[11] = make_histogram((aLabel+"_sumet").c_str(),"#Sigma E_{T}", "GeV", 100,0,5000);
  hist[12] = make_histogram((aLabel+"_ht").c_str(),"H_{T}", "GeV", 100,0,4000);
  hist[13] = make_histogram((aLabel+"_mht").c_str(),"MH_{T}", "GeV", 100,0,4000);
  hist[14] = make_histogram((aLabel+"_sqrt_ht").c_str(),"#sqrt{H_{T}}", "GeV^{1/2}", 100,0,70);
  hist[15] = make_histogram((aLabel+"_unclustered_et").c_str(),"Unclustered E_{T}", "GeV", 100,0,6000);
  hist[16] = make_histogram((aLabel+"_unclustered_phi").c_str(),"Unclustered #phi", "GeV", 100,-3.1416,3.1416);
  hist[17] = make_histogram((aLabel+"_jet1met_dphi").c_str(),"#Delta#phi(MET,jet1)", "", 100,0,3.1416);
  hist[18] = make_histogram((aLabel+"_jet2met_dphi").c_str(),"#Delta#phi(MET,jet2)", "", 100,0,3.1416);
  hist[19] = make_histogram((aLabel+"_jetmet_mindphi").c_str(),"minimum #Delta#phi(MET,jet)", "", 100,0,3.1416);
  hist[20] = make_histogram((aLabel+"_jetunclet_mindphi").c_str(),"minimum #Delta#phi(unclustered,jet)", "", 100,0,3.1416);
  hist[21] = make_histogram((aLabel+"_metunclet_dphi").c_str(),"#Delta#phi(MET,unclustered)", "", 100,0,3.1416);
  hist[22] = make_histogram((aLabel+"_dijetmet_scalarSum_pt").c_str(), "p_{T}^{jet1}+p_{T}^{jet2}+MET", "GeV" ,100,0,4000);
  hist[23] = make_histogram((aLabel+"_dijetmet_vectorialSum_pt").c_str(),"p_{T}(#vec{j1}+#vec{j2}+#vec{MET})", "GeV", 50,0,1000);
  hist[24] = make_histogram((aLabel+"_dijetmet_ptfraction").c_str(),"p_{T}^{dijet}/(p_{T}^{dijet}+MET)", "", 100,0,1);
  hist[25] = make_histogram((aLabel+"_jet1met_scalarprod").c_str(), "#vec{p_{T}^{jet1}}.#vec{MET}/MET", "GeV" ,100,-1500,1500);
  hist[26] = make_histogram((aLabel+"_jet2met_scalarprod").c_str(), "#vec{p_{T}^{jet2}}.#vec{MET}/MET", "GeV" ,100,-1500,1500);
  hist[27] = make_histogram((aLabel+"_jet1met_scalarprod_frac").c_str(), "#vec{p_{T}^{jet1}}.#vec{MET}/MET^{2}", "" ,50,-8,8);
  hist[28] = make_histogram((aLabel+"_jet2met_scalarprod_frac").c_str(), "#vec{p_{T}^{jet2}}.#vec{MET}/MET^{2}", "" ,50,-8,8);
  hist[29] = make_histogram((aLabel+"_n_jets_cjv_30").c_str(), "CJV jets (30 GeV)", "" ,15,0,15);
  hist[30] = make_histogram((aLabel+"_n_jets_cjv_20EB_30EE").c_str(), "CJV jets (20 EB, 30 EE)", "" ,15,0,15);

}

void FillVariables(TFile *inputFile, const unsigned nVars, TH1F * hist[], double & n_tot){

  double variable[nVars];
  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    variable[iV] = 0;
  }

  double weight = 1;
  double jet1_phi=0,jet2_phi=0,jet1_E=0,jet2_E=0,metx=0,mety=0;
  unsigned nCJV_30=0,nCJV_20_30=0;

  TTree *t = (TTree*)inputFile->Get("TmvaInputTree");
  t->SetBranchAddress("total_weight", &weight);
  t->SetBranchAddress("jet1_pt",                 &variable[0]);
  t->SetBranchAddress("jet2_pt",                 &variable[1]);
  t->SetBranchAddress("jet1_eta",                &variable[2]);
  t->SetBranchAddress("jet2_eta",                &variable[3]);
  t->SetBranchAddress("jet1_phi",                &jet1_phi);
  t->SetBranchAddress("jet2_phi",                &jet2_phi);
  t->SetBranchAddress("jet1_E",                  &jet1_E);
  t->SetBranchAddress("jet2_E",                  &jet2_E);
  t->SetBranchAddress("dijet_M",                 &variable[4]);
  t->SetBranchAddress("dijet_deta",              &variable[5]);
  t->SetBranchAddress("dijet_sumeta",            &variable[6]);
  t->SetBranchAddress("dijet_dphi",              &variable[7]);
  t->SetBranchAddress("met",                     &variable[8]);
  t->SetBranchAddress("met_x",                   &metx);
  t->SetBranchAddress("met_y",                   &mety);
  t->SetBranchAddress("met_phi",                 &variable[9]);
  t->SetBranchAddress("met_significance",        &variable[10]);
  t->SetBranchAddress("sumet",                   &variable[11]);
  t->SetBranchAddress("ht",                      &variable[12]);
  t->SetBranchAddress("mht",                     &variable[13]);
  t->SetBranchAddress("sqrt_ht",                 &variable[14]);
  t->SetBranchAddress("unclustered_et",          &variable[15]);
  t->SetBranchAddress("unclustered_phi",         &variable[16]);
  t->SetBranchAddress("jet1met_dphi",            &variable[17]);
  t->SetBranchAddress("jet2met_dphi",            &variable[18]);
  t->SetBranchAddress("jetmet_mindphi",          &variable[19]);
  t->SetBranchAddress("jetunclet_mindphi",       &variable[20]);
  t->SetBranchAddress("metunclet_dphi",          &variable[21]);
  t->SetBranchAddress("dijetmet_scalarSum_pt",   &variable[22]);
  t->SetBranchAddress("dijetmet_vectorialSum_pt",&variable[23]);
  t->SetBranchAddress("dijetmet_ptfraction",     &variable[24]);
  t->SetBranchAddress("jet1met_scalarprod",      &variable[25]);
  t->SetBranchAddress("jet2met_scalarprod",      &variable[26]);
  t->SetBranchAddress("n_jets_cjv_30",           &nCJV_30);
  t->SetBranchAddress("n_jets_cjv_20EB_30EE",    &nCJV_20_30);

  Int_t nevent = t->GetEntries();
  for (Int_t i=0;i<nevent;i++) {//loop on events
    t->GetEntry(i);

    //ROOT::Math::PtEtaPhiEVector jet1(variable[0],variable[2],jet1_phi,jet1_E);
    //ROOT::Math::PtEtaPhiEVector jet2(variable[1],variable[3],jet2_phi,jet2_E);
    //ROOT::Math::PtEtaPhiEVector metvec(variable[8],0,variable[9],variable[8]);
    //ROOT::Math::PtEtaPhiEVector unclustered(variable[15],0,variable[16],variable[15]);
    //variable[23] = (jet1.px()*metx+jet1.py()*mety)/variable[8];
    //variable[24] = (jet2.px()*metx+jet2.py()*mety)/variable[8];
    //variable[25] = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1,unclustered)),
    //fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2,unclustered)));
  //variable[26] = fabs(ROOT::Math::VectorUtil::DeltaPhi(unclustered,metvec));
    variable[27] = variable[25]/variable[8];
    variable[28] = variable[26]/variable[8];
    variable[29] = static_cast<double>(nCJV_30);
    variable[30] = static_cast<double>(nCJV_20_30);

    n_tot += weight;

    //cuts
    bool pass = true;
    //variable[5]>3.8 //deta cut
    //  && variable[4] > 1100 //Mjj cut
    //  && variable[8] > 100 //MET cut
    //  && variable[10] > 5. //MET significance
    //  ;
    if (pass) {
      for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
	hist[iV]->Fill(variable[iV],weight);
	//if (iV==0 && i<10) std::cout << inputFile->GetName() << " " << i << " jet1_pt = " << variable[iV] << " wt=" << weight << std::endl;
      }//loop on variables
    }
  }//loop on events
}

int CheckTmvaVariables() {//main

  std::string folder = "../output_tmva/nunu/MET130/";

  std::vector<std::string> datafiles;
  datafiles.push_back("Data_MET-2012A-13Jul2012-v1");
  datafiles.push_back("Data_MET-2012A-06Aug2012-v1");
  datafiles.push_back("Data_MET-2012B-13Jul2012-v1");
  datafiles.push_back("Data_MET-2012C-24Aug2012-v1");
  datafiles.push_back("Data_MET-2012C-11Dec2012-v1"); 
  datafiles.push_back("Data_MET-2012C-PromptReco-v2");
  datafiles.push_back("Data_MET-2012D-PromptReco-v1");

  //all files
  std::vector<std::string> qcdfiles;
  qcdfiles.push_back("MC_QCD-Pt-30to50-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-50to80-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-80to120-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-120to170-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-170to300-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-300to470-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-470to600-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-600to800-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-800to1000-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-1000to1400-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-1400to1800-pythia6");
  qcdfiles.push_back("MC_QCD-Pt-1800-pythia6");
  qcdfiles.push_back("MC_GJets-HT-200To400-madgraph");
  qcdfiles.push_back("MC_GJets-HT-400ToInf-madgraph");

  //all bkg without qcd
  std::vector<std::string> files;
  files.push_back("MC_TTJets");
  //powheg samples
  //files.push_back("MC_TT-v1");
  //files.push_back("MC_TT-v2");
  //
  files.push_back("MC_T-tW");
  files.push_back("MC_Tbar-tW");
  files.push_back("MC_SingleT-s-powheg-tauola");
  files.push_back("MC_SingleTBar-s-powheg-tauola");
  files.push_back("MC_SingleT-t-powheg-tauola");
  files.push_back("MC_SingleTBar-t-powheg-tauola");
  files.push_back("MC_WW-pythia6-tauola");
  files.push_back("MC_WZ-pythia6-tauola");
  files.push_back("MC_ZZ-pythia6-tauola");
  files.push_back("MC_W1JetsToLNu_enu");
  files.push_back("MC_W2JetsToLNu_enu");
  files.push_back("MC_W3JetsToLNu_enu");
  files.push_back("MC_W4JetsToLNu_enu");
  files.push_back("MC_WJetsToLNu-v1_enu");
  files.push_back("MC_WJetsToLNu-v2_enu");
  files.push_back("MC_W1JetsToLNu_munu");
  files.push_back("MC_W2JetsToLNu_munu");
  files.push_back("MC_W3JetsToLNu_munu");
  files.push_back("MC_W4JetsToLNu_munu");
  files.push_back("MC_WJetsToLNu-v1_munu");
  files.push_back("MC_WJetsToLNu-v2_munu");
  files.push_back("MC_W1JetsToLNu_taunu");
  files.push_back("MC_W2JetsToLNu_taunu");
  files.push_back("MC_W3JetsToLNu_taunu");
  files.push_back("MC_W4JetsToLNu_taunu");
  files.push_back("MC_WJetsToLNu-v1_taunu");
  files.push_back("MC_WJetsToLNu-v2_taunu");
  files.push_back("MC_DYJetsToLL");
  files.push_back("MC_DY1JetsToLL");
  files.push_back("MC_DY2JetsToLL");
  files.push_back("MC_DY3JetsToLL");
  files.push_back("MC_DY4JetsToLL");
  files.push_back("MC_ZJetsToNuNu_100_HT_200");
  files.push_back("MC_ZJetsToNuNu_200_HT_400");
  files.push_back("MC_ZJetsToNuNu_400_HT_inf");
  files.push_back("MC_ZJetsToNuNu_50_HT_100");
  files.push_back("MC_WGamma");
  files.push_back("MC_EWK-Z2j");
  files.push_back("MC_EWK-Z2jiglep");
  files.push_back("MC_EWK-W2jminus_enu");
  files.push_back("MC_EWK-W2jplus_enu");
  files.push_back("MC_EWK-W2jminus_munu");
  files.push_back("MC_EWK-W2jplus_munu");
  files.push_back("MC_EWK-W2jminus_taunu");
  files.push_back("MC_EWK-W2jplus_taunu");

  const unsigned nData = datafiles.size();
  const unsigned nQCD = qcdfiles.size();
  const unsigned nMC = files.size();
  TFile *fData[nData];
  TFile *fQCD[nQCD];
  TFile *fMC[nMC];

  TFile *fSignal120 = TFile::Open((folder+"/MC_VBF_HToZZTo4Nu_M-120.root").c_str());
  TFile *fSignal400 = TFile::Open((folder+"/MC_VBF_HToZZTo4Nu_M-400.root").c_str());

  //TFile *output_withQCD = new TFile("DataMC_SumBkgWithQCD.root","RECREATE");
  //TFile *output_withoutQCD = new TFile("DataMC_SumBkgWithoutQCD.root","RECREATE");

  const unsigned nVars = 31;

  TCanvas *myc[nVars];
  gStyle->SetOptStat(0);


  bool logy[nVars] = {
    1,1,0,0,1,
    0,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,0,
    1,1,1,1,1,
    1
  };

  double n_dijet_data = 0;
  double n_dijet_bkg = 0;
  double n_dijet_qcd = 0;
  double n_dijet_signal = 0;
  double n_deta_mjj_met_data = 0;
  double n_deta_mjj_met_bkg = 0;
  double n_deta_mjj_met_qcd = 0;
  double n_deta_mjj_met_signal = 0;

  for (unsigned iV(0); iV<nVars;++iV){
    std::ostringstream canvasName;
    canvasName << "myc_" << iV;
    myc[iV] = new TCanvas(canvasName.str().c_str(),"",1);
    myc[iV]->SetLogy(logy[iV]);
  }

  TH1F *hVar[nVars];
  MakeAllHistograms(hVar,"Bkg");

  TH1F *hVarSignal120[nVars];
  MakeAllHistograms(hVarSignal120,"Signal120");
  FillVariables(fSignal120,nVars,hVarSignal120,n_dijet_signal);

  TH1F *hVarSignal400[nVars];
  MakeAllHistograms(hVarSignal400,"Signal400");
  double n_dijet_400 = 0;
  FillVariables(fSignal400,nVars,hVarSignal400,n_dijet_400);

  TH1F *hVarData[nVars];
  MakeAllHistograms(hVarData,"Data");

  TH1F *hVarQCD[nVars];
  MakeAllHistograms(hVarQCD,"QCD");
  

  for (unsigned iBkg = 0; iBkg < files.size(); ++iBkg) {//loop on bkg files
    fMC[iBkg] = TFile::Open((folder+"/"+files[iBkg]+".root").c_str());
    if (!fMC[iBkg]) {
      std::cerr << "Warning, file " << files[iBkg] << " could not be opened." << std::endl;
      continue;//return 1;
    }

    std::cout << " Processing tree " <<  files[iBkg] << std::endl;
    
    FillVariables(fMC[iBkg],nVars,hVar,n_dijet_bkg);

  }//loop on bkg files

  for (unsigned iD(0); iD<nData;++iD){//loop on data files
    fData[iD] = TFile::Open((folder+"/"+datafiles[iD]+".root").c_str());
    if (!fData[iD]) {
      std::cerr << "Warning, file " << datafiles[iD] << " could not be opened." << std::endl;
      continue;//return 1;
    }   
    std::cout << " Processing tree " <<  datafiles[iD] << std::endl;
 
    FillVariables(fData[iD],nVars,hVarData,n_dijet_data);
  }//loop on data files

  for (unsigned iQ(0); iQ<nQCD;++iQ){//loop on qcd files
    fQCD[iQ] = TFile::Open((folder+"/"+qcdfiles[iQ]+".root").c_str());
    if (!fQCD[iQ]) {
      std::cerr << "Warning, file " << qcdfiles[iQ] << " could not be opened." << std::endl;
      continue;//return 1;
    }   
    std::cout << " Processing tree " <<  qcdfiles[iQ] << std::endl;
 
    FillVariables(fQCD[iQ],nVars,hVarQCD,n_dijet_qcd);
  }//loop on qcd files

  n_deta_mjj_met_qcd = hVarQCD[0]->Integral();
  n_deta_mjj_met_bkg = hVar[0]->Integral();
  n_deta_mjj_met_data = hVarData[0]->Integral();
  n_deta_mjj_met_signal = hVarSignal120[0]->Integral();

  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    myc[iV]->cd();
  

    //hVarQCD[iV]->Add(hVar[iV]);
    hVarSignal120[iV]->Scale(hVarData[iV]->Integral()/hVarSignal120[iV]->Integral());
    //hVarSignal400[iV]->Scale(hVarQCD[iV]->Integral()/hVarSignal400[iV]->Integral());
    hVar[iV]->Scale(hVarData[iV]->Integral()/hVar[iV]->Integral());
    hVarQCD[iV]->Scale(hVarData[iV]->Integral()/hVarQCD[iV]->Integral());

    hVarData[iV]->SetMarkerColor(1);
    hVarData[iV]->SetMarkerStyle(2);
    hVarData[iV]->GetYaxis()->SetTitle("Events");
    if (hVarData[iV]->GetMinimum() > hVar[iV]->GetMinimum()) hVarData[iV]->SetMinimum(hVar[iV]->GetMinimum());
    if (hVarData[iV]->GetMaximum() < hVarQCD[iV]->GetMaximum()) hVarData[iV]->SetMaximum(hVarQCD[iV]->GetMaximum());
    if (hVarData[iV]->GetMaximum() < hVarSignal120[iV]->GetMaximum()) hVarData[iV]->SetMaximum(hVarSignal120[iV]->GetMaximum());
    //if (hVarData[iV]->GetMinimum() > hVarSignal120[iV]->GetMinimum()) hVarData[iV]->SetMinimum(hVarSignal120[iV]->GetMinimum());
    if (logy[iV]==true && hVarData[iV]->GetMinimum()==0) hVarData[iV]->SetMinimum(0.01);
    hVarData[iV]->Draw("PE");
    //plot QCD+others
     hVarQCD[iV]->SetFillColor(5);
    hVarQCD[iV]->SetLineColor(5);
    hVarQCD[iV]->GetYaxis()->SetTitle("Events");
    hVarQCD[iV]->Draw("same");
    //plot others
    hVar[iV]->SetFillColor(6);
    hVar[iV]->SetLineColor(6);
    hVar[iV]->SetFillStyle(3005);
    hVar[iV]->GetYaxis()->SetTitle("Events");
    hVar[iV]->Draw("same");
    hVar[iV]->Draw("Psame");

    //draw signal120
    hVarSignal120[iV]->SetFillColor(4);
    hVarSignal120[iV]->SetFillStyle(3004);
    hVarSignal120[iV]->SetLineColor(4);
    hVarSignal120[iV]->SetMarkerColor(4);
    hVarSignal120[iV]->SetMarkerStyle(22);
    hVarSignal120[iV]->GetYaxis()->SetTitle("Events");
    hVarSignal120[iV]->Draw("same");
    hVarSignal120[iV]->Draw("Psame");
    //draw signal400
    //hVarSignal400[iV]->SetFillColor(6);
    hVarSignal400[iV]->SetLineColor(4);
    //hVarSignal400[iV]->SetMarkerColor(4);
    //hVarSignal400[iV]->SetMarkerStyle(23);
    hVarSignal400[iV]->GetYaxis()->SetTitle("Events");
    //hVarSignal400[iV]->Draw("Lsame");

    //replot data on top
    hVarData[iV]->Draw("PEsame");
    hVarData[iV]->Draw("axissame");

    //double proba = hVarData[iV]->KolmogorovTest(hVarQCD[iV],"UON");
    //std::cout << hVar[iV]->GetName() << " = " << proba << std::endl;

    TLegend *leg = new TLegend(0.75,0.8,1.0,1.0);
    leg->SetFillColor(10);
    leg->AddEntry(hVarData[iV],"Data","P");
    leg->AddEntry(hVarSignal120[iV],"VBF m_{H}=120 GeV","F");
    //leg->AddEntry(hVarSignal400[iV],"VBF m_{H}=400 GeV","L");
    leg->AddEntry(hVarQCD[iV],"QCD","F");
    leg->AddEntry(hVar[iV],"V+top+VV","F");
    leg->Draw("same");

    std::ostringstream output;
    //output << "TmvaPlots_deta3p8_Mjj1100_met100_metsig5/" << hVar[iV]->GetName() << ".pdf";
    output << "TmvaPlots_scaleEwkBkg/" << hVar[iV]->GetName() << ".pdf";
    myc[iV]->Print(output.str().c_str() );


  }//loop on variables

  std::cout << "----------------------------------------------" << std::endl
	    << "---- Summary of number of events selected ----" << std::endl
	    << " - N_dijet = " << n_dijet_data << " qcd = " << n_dijet_qcd << " bkg = " << n_dijet_bkg << " signal = " << n_dijet_signal << std::endl
	    <<" - N_deta_mjj_met = " << n_deta_mjj_met_data << " qcd = " << n_deta_mjj_met_qcd << " bkg = " << n_deta_mjj_met_bkg << " signal = " << n_deta_mjj_met_signal 
	    << std::endl;

  return 0;
}//main
