#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH2F.h" 
#include "TGraphErrors.h" 
#include "TDirectory.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TPad.h"

int sfMinDPhi() {//main

  unsigned region = 13;
  std::string filedir = "../output_allmindphi_1_metsig_3_4_mindphi_1_2/";
  if (region==3) filedir = "../output_allmindphi_1_metsig_4_mindphi_1_2/";
  else if (region==12) filedir = "../output_allmindphi_1_metsig_3_4_mindphi_1/";
  else if (region==13) filedir = "../output_allmindphi_1_metsig_3_mindphi_1_2/";
  std::string channel = "nunu";
  std::string histname = "jetmetnomu_mindphi";

  double minX = 1;
  double maxX = 3;


  std::ostringstream suffix;
  suffix << "_norm" << region;

  std::string plotDir = "PLOTS/";

  std::string filename = filedir+channel+".root";
  TFile *inputf = TFile::Open(filename.c_str());
  if (!inputf) {
    std::cout << " Input file " << filename << " not found. Exiting..." << std::endl;
    return 1;
  }

  const unsigned nS = 3;
  TH1F *hist[nS+1];

  const unsigned nDirs = 9;
  std::string dirs[nDirs] = {
    "data_obs","qcd",
    "wmu","wel","wtau",
    "zvv","vv","wg",
    "top"
  };
  
  for (unsigned iD(0); iD<nDirs;++iD){
    inputf->cd(dirs[iD].c_str());
    if (iD==1) hist[iD] = (TH1F*)gDirectory->Get("jetmetnomu_mindphi")->Clone();
    else if (iD<nS) hist[iD] = (TH1F*)gDirectory->Get(histname.c_str())->Clone();
    else hist[nS-1]->Add((TH1F*)gDirectory->Get(histname.c_str()));
  }

  inputf->cd("qqH");
  hist[nS] = (TH1F*)gDirectory->Get(histname.c_str())->Clone();
  inputf->cd("ggH");
  hist[nS]->Add((TH1F*)gDirectory->Get(histname.c_str())->Clone());

  TCanvas *myc = new TCanvas("myc","myc",1500,1000);
  myc->cd();

  for (unsigned iS(0); iS<nS+1;++iS){//loop on samples
    hist[iS]->SetLineColor(iS+1);
    hist[iS]->SetMarkerColor(iS+1);
    hist[iS]->SetMarkerStyle(2);
    if (iS>=nS-1) {
      hist[iS]->SetFillStyle(1001);
      hist[iS]->SetFillColor(iS+1);
    }
    hist[iS]->SetTitle((";"+histname+";entries").c_str());
    //hist[iS]->GetXaxis()->SetRangeUser(0,2);
  }//loop on samples


  gStyle->SetOptStat(0);
  TH1F *dataminusbkg = (TH1F*)hist[0]->Clone("hsubtr");
  dataminusbkg->Add(hist[2],-1.);
  //copy scale factors from job output!
  if (region==1) hist[1]->Scale(1/0.219);
  else if (region==3) hist[1]->Scale(1/1.216);
  else if (region==12) hist[1]->Scale(1/0.153);
  else if (region==13) hist[1]->Scale(1/0.412);

  TH1F *ratio = (TH1F*)dataminusbkg->Clone("ratio");
  ratio->Divide(dataminusbkg,hist[1]);
  TH1F *ratioint = (TH1F*)dataminusbkg->Clone("ratioint");
  int nbins = ratio->GetNbinsX()+1;
  for (int i(1);i<nbins+1;++i){
    if (dataminusbkg->GetBinContent(i)<0) dataminusbkg->SetBinContent(i,0);
  }
  for (int i(1);i<nbins+1;++i){
    double val = 0;
    double valerr = 0;
    double errnum = 0;
    double errden = 0;
    double intnum = dataminusbkg->IntegralAndError(i,nbins,errnum);
    double intden = hist[1]->IntegralAndError(i,nbins,errden);
    if (intnum>0 && intden>0 && hist[1]->GetBinLowEdge(i)<maxX){
      val = intnum/intden;
      valerr = val*sqrt(pow(errnum/intnum,2)+pow(errden/intden,2));
    }
    if (i==1) std::cout << " Total integrals: data-bkg = " << intnum << " qcd-bkg = " << intden << std::endl;
    ratioint->SetBinContent(i,val);
    ratioint->SetBinError(i,valerr);
  }

  myc->Divide(1,2);
  myc->cd(1);
  TPad *top = (TPad*)gPad;
  top->Divide(2,1);
  top->cd(1);
  hist[1]->GetXaxis()->SetRangeUser(minX,maxX);
  hist[1]->SetMinimum(0);
  hist[1]->Draw();
  dataminusbkg->Draw("same");
  top->cd(2);
  ratio->GetXaxis()->SetRangeUser(minX,maxX);
  ratio->GetYaxis()->SetTitle("data/QCD");
  ratio->Draw();
  myc->cd(2);
  gPad->SetGridy(1);
  ratioint->GetYaxis()->SetTitle("QCD SF");
  ratioint->GetXaxis()->SetRangeUser(minX,3.2);
  ratioint->Draw();
  gStyle->SetOptFit(1111);

  TF1 *fit = new TF1("fit","exp([0]+[1]*x)",1,3.2);
  fit->SetParameters(-0.029,-1.8);
  double res=0;
  double reserr=0;
  double res2=0;
  double reserr2=0;
  if (region != 3 && region != 13){
    ratioint->Fit("fit","+","same",minX,maxX);
    //fit->Draw("same");
    res = fit->Integral(2.5,2.6)/0.1;//3.1416);
    reserr = fit->IntegralError(2.5,2.6)/0.1;//3.1416);
    res2 = fit->Integral(2.,2.1)/0.1;//3.1416);
    reserr2 = fit->IntegralError(2.,2.1)/0.1;//3.1416);
  } else {
    ratioint->Fit("pol1","+","same",minX,maxX);
    res = ratioint->GetFunction("pol1")->Integral(2.5,2.6)/0.1;//3.1416);
    reserr = ratioint->GetFunction("pol1")->IntegralError(2.5,2.6)/0.1;//3.1416);
    res2 = ratioint->GetFunction("pol1")->Integral(2.0,2.1)/0.1;//3.1416);
    reserr2 = ratioint->GetFunction("pol1")->IntegralError(2.0,2.1)/0.1;//3.1416);
  }
  std::cout << "Fit integral 2.-Pi: " << res2 << std::endl;
  std::cout << "Fit integral 2.5-Pi: " << res << std::endl;

  TLatex lat;
  lat.SetTextSize(0.05);
  char buf[500];

  sprintf(buf,"Expectation 2.-#pi: SF=%3.3f #pm %3.3f",res2,reserr2);
  lat.DrawLatex(2.2,0.06,buf);

  sprintf(buf,"Expectation 2.5-#pi: SF=%3.3f #pm %3.3f",res,reserr);
  lat.DrawLatex(2.2,res,buf);


  myc->Update();
  myc->Print((plotDir+histname+suffix.str()+"_SF.pdf").c_str());


  return 0;
}//main
