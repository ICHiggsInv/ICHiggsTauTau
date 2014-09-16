#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h" 
#include "TH2F.h" 
#include "TGraphErrors.h" 
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

struct Event {
  unsigned run;
  unsigned event;
  unsigned lumi;

  bool operator==(const Event & rhs)const{
    if(this!=&rhs){
      if(this->run == rhs.run && this->lumi == rhs.lumi && this->event==rhs.event){
        return true;
      }
      return false; 
    }
    return true; 
  }

  bool operator<(const Event & rhs)const {
    if(this==&rhs){
      return false; 
    }

    if(this->run < rhs.run){
      return true; 
    }
    else if (this->run == rhs.run && this->lumi < rhs.lumi) return true;
    else if (this->run == rhs.run && this->lumi == rhs.lumi && this->event<rhs.event) return true;
    return false; 
  }

};

double qcdWeight(double aVar){
  //jet1 pT
  /*
  if (aVar<=70) return 0;
  else if (aVar>70 && aVar<=80) return 0.77139;// 0.120061
  else if (aVar>80 && aVar<=90) return 0.874573;// 0.0820615
  else if (aVar>90 && aVar<=100) return 0.950359;// 0.0653604
  else if (aVar>100 && aVar<=110) return 1.00269;// 0.0540433
  else if (aVar>110 && aVar<=120) return 0.734874;// 0.0400591
  else if (aVar>120 && aVar<=130) return 0.679982;// 0.034953
  else if (aVar>130 && aVar<=140) return 0.480797;// 0.0265533
  else if (aVar>140 && aVar<=150) return 0.420656;// 0.0238441
  else if (aVar>150 && aVar<=160) return 0.310657;// 0.0206916
  else if (aVar>160 && aVar<=170) return 0.243469;// 0.0179453
  else if (aVar>170 && aVar<=180) return 0.204945;// 0.0169434
  else if (aVar>180 && aVar<=190) return 0.182669;// 0.0158033
  else if (aVar>190 && aVar<=200) return 0.13522 ;//0.0145337
  else if (aVar>200 && aVar<=210) return 0.110492;// 0.014156
  else if (aVar>210 && aVar<=220) return 0.0844846;// 0.0129247
  else if (aVar>220 && aVar<=230) return 0.0720299;// 0.0139459
  else if (aVar>230 && aVar<=240) return 0.0703659;// 0.0148048
  else if (aVar>240 && aVar<=250) return 0.0632627;// 0.0141989
  else if (aVar>250 && aVar<=260) return 0.0676825;// 0.0151234
  else if (aVar>260 && aVar<=270) return 0.0398887;// 0.016641
  else if (aVar>270 && aVar<=280) return 0.0213198;// 0.0173002
  else if (aVar>280 && aVar<=290) return 0.0536532;// 0.0168529
  else if (aVar>290) return 0.0483656;// 0.0160575
  else return 1;
  */
  //dphijj
  if (aVar>0 &&  aVar <=0.1) return 0.143877;//// 0.0162334
  else if (aVar>0.1 &&  aVar <=0.2) return 0.142455;//// 0.0162577
  else if (aVar>0.2 &&  aVar <=0.3) return 0.154856;//// 0.0166426
  else if (aVar>0.3 &&  aVar <=0.4) return 0.164215;//// 0.0166573
  else if (aVar>0.4 &&  aVar <=0.5) return 0.151531;//// 0.0162367
  else if (aVar>0.5 &&  aVar <=0.6) return 0.162384;//// 0.0163704
  else if (aVar>0.6 &&  aVar <=0.7) return 0.173669;//// 0.0174333
  else if (aVar>0.7 &&  aVar <=0.8) return 0.155771;//// 0.0163837
  else if (aVar>0.8 &&  aVar <=0.9) return 0.171223;//// 0.0173506
  else if (aVar>0.9 &&  aVar <=1) return 0.186394;//// 0.0175418
  else if (aVar>1 &&  aVar <=1.1) return 0.2224;//// 0.0184566
  else if (aVar>1.1 &&  aVar <=1.2) return 0.203679;//// 0.0180246
  else if (aVar>1.2 &&  aVar <=1.3) return 0.23474;//// 0.0182477
  else if (aVar>1.3 &&  aVar <=1.4) return 0.268086;//// 0.0195426
  else if (aVar>1.4 &&  aVar <=1.5) return 0.276632;//// 0.0207171
  else if (aVar>1.5 &&  aVar <=1.6) return 0.317397;//// 0.02437
  else if (aVar>1.6 &&  aVar <=1.7) return 0.34657;//// 0.0265896
  else if (aVar>1.7 &&  aVar <=1.8) return 0.473792;//// 0.0342175
  else if (aVar>1.8 &&  aVar <=1.9) return 0.641095;//// 0.0491317
  else if (aVar>1.9 &&  aVar <=2) return 0.789086;//// 0.0685297
  else if (aVar>2 &&  aVar <=2.1) return 1.08172;//// 0.0963207
  else if (aVar>2.1 &&  aVar <=2.2) return 1.16063;//// 0.133661
  else if (aVar>2.2 &&  aVar <=2.3) return 1.25094;//// 0.166526
  else if (aVar>2.3 &&  aVar <=2.4) return 1.27074;//// 0.197914
  else if (aVar>2.4 &&  aVar <=2.5) return 1.0452;//// 0.220793
  else if (aVar>2.5 &&  aVar <=2.6) return 0.93385;//// 0.300231
  else if (aVar>2.6 &&  aVar <=2.7) return 0.706125;//// 0.287712
  else if (aVar>2.7 &&  aVar <=2.8) return 1.24828;//// 0.380992
  else if (aVar>2.8 &&  aVar <=2.9) return 0.80586;//// 0.38304
  else if (aVar>2.9 &&  aVar <=3) return 0.267529;//// 0.310501
  

  return 1;

};

int dataDrivenQCD() {

  std::string fileName = "../../output_lighttree/VBFPARKED.root";//_MET-2012A-22Jan2013-v1.root";
  //std::string fileName = "../../output/MC_Powheg-Htoinv-mH125.root";

  std::string type = "DataMC_PARKED";//_rwdphijj";
  //std::string type = "VBFH125";

  unsigned varIdx = 3;

  const unsigned nTrees = 3;
  std::string label[nTrees] = {"j1j2","j1j3","j2j3"};

  TFile *mcFile[nTrees];
  mcFile[0] = TFile::Open("../output/nunu.root");
  mcFile[1] = TFile::Open("../output/nunu_J1J3.root");
  mcFile[2] = TFile::Open("../output/nunu_J2J3.root");

  TFile *data = TFile::Open(fileName.c_str(), "update");
  if (!data) {
    std::cout << " Input file " << fileName << " not found. Exiting..." << std::endl;
    return 1;
  }

  data->cd();


  TTree *tree[nTrees];
  tree[0] = (TTree*)gDirectory->Get("LightTree");
  tree[1] = (TTree*)gDirectory->Get("LightTreeJ1J3");
  tree[2] = (TTree*)gDirectory->Get("LightTreeJ2J3");

  unsigned run;
  unsigned lumi;
  unsigned event;
  double passtrigger;
  double passparkedtrigger1;
  double passparkedtrigger2;
  double l1met;
  double jet1_pt;
  double jet3_pt;
  double dijet_M;
  double dijet_dphi;
  double jet1_eta;
  double jet2_eta;
  double metnomuons;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;
  unsigned n_jets;
  int is_dupl = 0;
  double qcdW = 1.;

  std::set<Event> evtSet;
  std::pair<std::set<Event>::iterator,bool> isInserted;
  unsigned duplicateJ1J2 = 0;
  unsigned duplicateJ1J3 = 0;
  unsigned duplicateJ2J3 = 0;
  unsigned passJ1J2 = 0;
  unsigned passJ1J3 = 0;
  unsigned passJ2J3 = 0;
  TBranch *brDupl[nTrees];
  TBranch *brQCDWeight[nTrees];

  int nEntries[nTrees];
  
  for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
    brDupl[iT] = tree[iT]->Branch("is_dupl", &is_dupl, "is_dupl/I");
    brQCDWeight[iT] = tree[iT]->Branch("qcdW", &qcdW, "qcdW/D");

    nEntries[iT] = tree[iT]->GetEntries();
    tree[iT]->SetBranchAddress("run",&run);
    tree[iT]->SetBranchAddress("lumi",&lumi);
    tree[iT]->SetBranchAddress("event",&event);
    tree[iT]->SetBranchAddress("passtrigger",&passtrigger);
    tree[iT]->SetBranchAddress("passparkedtrigger1",&passparkedtrigger1);
    tree[iT]->SetBranchAddress("passparkedtrigger2",&passparkedtrigger2);
    tree[iT]->SetBranchAddress("l1met",&l1met);
    tree[iT]->SetBranchAddress("jet1_pt",&jet1_pt);
    tree[iT]->SetBranchAddress("jet3_pt",&jet3_pt);
    tree[iT]->SetBranchAddress("dijet_M",&dijet_M);
    tree[iT]->SetBranchAddress("dijet_dphi",&dijet_dphi);
    tree[iT]->SetBranchAddress("jet1_eta",&jet1_eta);
    tree[iT]->SetBranchAddress("jet2_eta",&jet2_eta);
    tree[iT]->SetBranchAddress("metnomuons",&metnomuons);
    tree[iT]->SetBranchAddress("nvetomuons",&nvetomuons);
    tree[iT]->SetBranchAddress("nselmuons",&nselmuons);
    tree[iT]->SetBranchAddress("nvetoelectrons",&nvetoelectrons);
    tree[iT]->SetBranchAddress("nselelectrons",&nselelectrons);
    tree[iT]->SetBranchAddress("n_jets",&n_jets);

    std::cout << " Jet pair " << label[iT] << " has " << nEntries[iT] << " entries in tree." << std::endl;
    
    for (int iE(0); iE<nEntries[iT]; ++iE){//loop on entries
      tree[iT]->GetEntry(iE);
      is_dupl = 0;
      bool passtrig = ((run>=190456 && run<=193621 &&passtrigger==1) || (run>=193833 && run<=196531 && passparkedtrigger1==1) ||(run>=203777 && run<=208686 && passparkedtrigger2==1)) && l1met>40;//parked
      bool passpT = (iT==0 && jet1_pt > 50) || iT>0;
      bool passnj = (iT==0) || (iT>0 && n_jets>2);
      qcdW = qcdWeight(dijet_dphi);

      if (passtrig && passpT && passnj && jet1_eta*jet2_eta < 0 && dijet_M>=600 && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0){//pass sel
	Event lEvt;
	lEvt.run = run;
	lEvt.event = event;
	lEvt.lumi = lumi;
	isInserted = evtSet.insert(lEvt);
	if (!isInserted.second){
	  is_dupl = 1;
	  if (iT==0)
	    duplicateJ1J2++;
	  else if (iT==1)
	    duplicateJ1J3++;
	  else 
	    duplicateJ2J3++;
	}
	else {
	  if (iT==0)
	    passJ1J2++;
	  else if (iT==1)
	    passJ1J3++;
	  else 
	    passJ2J3++;
	}
      }//pass sel
      //if (iE<10){
      //std::cout << iT << " " << iE << " " << is_dupl << " " << qcdW << std::endl;
      //}
      brDupl[iT]->Fill();
      brQCDWeight[iT]->Fill();
    }//loop on entries
  }//loop on trees


  // tree[0]->SetBranchAddress("is_dupl", &is_dupl);
  // tree[0]->SetBranchAddress("qcdW", &qcdW);
  // for (int iE(0); iE<10; ++iE){//loop on entries
  //   tree[0]->GetEntry(iE);
  //   std::cout << 0 << " " << iE << " " << is_dupl << " " << qcdW << std::endl;
  // }
  // //return 1;

  std::cout << " Duplicated/passing events: " << std::endl
	    << " j1j2 & " << duplicateJ1J2+passJ1J2 << " & " << duplicateJ1J2 << " & " << passJ1J2 << "\\\\" << std::endl
	    << " j1j3 & " << duplicateJ1J3+passJ1J3 << " & " << duplicateJ1J3 << " & " << passJ1J3 << "\\\\"<< std::endl
	    << " j2j3 & " << duplicateJ2J3+passJ2J3 << " & " << duplicateJ2J3 << " & " << passJ2J3 << "\\\\"<< std::endl
    ;

  const unsigned nVars = 17;
  std::string vars[nVars] = {
    "dijet_M",
    "metnomuons",	     
    "dijet_deta",
    "dijet_dphi",
    "metnomu_significance",
    "jetmetnomu_mindphi",
    "jet1_pt",
    "jet2_pt",
    "jet3_pt",
    "jet1_csv",
    "jet2_csv",
    "jet3_csv",
    "dijetmetnomu_scalarSum_pt",
    "dijetmetnomu_vectorialSum_pt",
    "dijetmetnomu_ptfraction",
    "n_jets",
    "n_jets_cjv_30"
  };

  std::string xaxis[nVars] = {
    ";M_{jj} (GeV)",
    ";METnoMu (GeV)",
    ";#Delta#eta_{jj}",
    ";#Delta#phi_{jj}",
    ";METnoMu/#sigma(METnoMu)",
    ";min #Delta#phi(j,METnoMu)",
    ";p_{T}^{j1} (GeV)",
    ";p_{T}^{j2} (GeV)",
    ";p_{T}^{j3} (GeV)",
    ";CSV (jet1)",
    ";CSV (jet2)",
    ";CSV (jet3)",
    ";p_{T}^{jet1}+p_{T}^{jet2}+METnoMu",
    ";p_{T}(#vec{j1}+#vec{j2}+#vec{METnoMu})",
    ";p_{T}^{dijet}/(p_{T}^{dijet}+METnoMu)",
    ";n_{jets} (p_{T}>30 GeV)",
    ";CJV jets (30 GeV)"
  };

  //unsigned nbins[nVars] = {75,100,44,30,35,30,50,50,100,50,50,10};
  float min[nVars] = {0,0,3.6,0,3,1.5,30,30,30,0,0,0,0,0,0,0,0};
  float max[nVars] = {3000,500,8,3.1416,10,3.1416,300,300,300,1,1,1,1000,400,1,10,10};  

  TH1F *hist[nVars][nTrees];
  TH1F *histsum[nVars];
  TH2F *weight2D = new TH2F("weight2D",";pT1 (GeV); pT2 (GeV)",
			    27,30,300,
			    27,30,300);

  TGraphErrors *grWeight[nVars];
  TH1F *histMC[nVars][nTrees];
  TH1F *histDataCheck[nVars];
  TH1F *histDataSubtr[nVars][nTrees];
  TH1F *histsumSubtr[nVars];
  TH1F *histSig[nVars];
  TH1F *histQCD[nVars];

  const unsigned nDirs = 7;
  std::string dirs[nDirs] = {
    "wmu","wel","wtau",
    "zvv","vv","wg",
    "top"
  };


  TCanvas *myc[nVars];
  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    std::ostringstream lname;
    lname.str("");
    lname << "myc_" << vars[iV];
    myc[iV] = new TCanvas(lname.str().c_str(),
			  vars[iV].c_str(),
			  1600,800);

    //get MC shapes
    mcFile[0]->cd("qcd");
    histQCD[iV] = (TH1F*)gDirectory->Get(vars[iV].c_str())->Clone();
    histQCD[iV]->SetLineColor(6);
    histQCD[iV]->SetFillColor(6);
    histQCD[iV]->SetFillStyle(3005);
    
    mcFile[0]->cd("data_obs");
    histDataCheck[iV] = (TH1F*)gDirectory->Get(vars[iV].c_str())->Clone();
    histDataCheck[iV]->SetLineColor(7);
    histDataCheck[iV]->SetMarkerColor(7);
    histDataCheck[iV]->SetMarkerStyle(23);
    
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      //get MC
      for (unsigned iD(0); iD<nDirs;++iD){
      	mcFile[iT]->cd(dirs[iD].c_str());
      	if (iD==0) histMC[iV][iT] = (TH1F*)gDirectory->Get(vars[iV].c_str())->Clone();
      	else histMC[iV][iT]->Add((TH1F*)gDirectory->Get(vars[iV].c_str()));
      }
      histMC[iV][iT]->SetLineColor(5);
      histMC[iV][iT]->SetLineWidth(3);
      histMC[iV][iT]->SetFillColor(5);
      histMC[iV][iT]->SetFillStyle(3005);
      
      if (iV==0) std::cout << "tree " << iT 
      			   << " nMC = " 
      			   << histMC[iV][iT]->Integral()
      			   << std::endl;

      mcFile[iT]->cd("qqH");
      histSig[iV] = (TH1F*)gDirectory->Get(vars[iV].c_str())->Clone();
      histSig[iV]->SetLineColor(6);
      
      data->cd();
      // tree[iT]->SetBranchAddress("is_dupl", &is_dupl);
      // tree[iT]->SetBranchAddress("qcdW", &qcdW);
      // for (int iE(0); iE<10; ++iE){//loop on entries
      // 	tree[iT]->GetEntry(iE);
      // 	std::cout << iT << " " << iE << " " << is_dupl << " " << qcdW << std::endl;
      // }
      // return 1;

      std::ostringstream selection;
      std::string passtrig = "(((run>=190456 && run<=193621 && passtrigger==1) || (run>=193833 && run<=196531 &&passparkedtrigger1==1) || (run>=203777 && run<=208686 && passparkedtrigger2==1)) && l1met>40)";
      std::string passpT = "(jet1_pt > 50)";

      selection << "( ";
      if (iT>0) selection << "n_jets>2 && ";
      selection << "is_dupl==0 && ";
      selection << passtrig <<" && ";
      if (iT==0) selection << passpT << " && ";
      selection << "jet1_eta*jet2_eta<0 && dijet_M>=600 && (jet1_pt>50 || jet3_pt>50) && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      selection << " )";
      //if (iT>0) selection << " * (qcdW)";

      lname.str("");
      lname << "p_" << vars[iV] << "_" << label[iT];
      //hist[iV][iT] = new TH1F(lname.str().c_str(),xaxis[iV].c_str(),nbins[iV],min[iV],max[iV]);
      hist[iV][iT] = (TH1F*)histDataCheck[iV]->Clone(lname.str().c_str());
      hist[iV][iT]->Reset();
      hist[iV][iT]->SetTitle(xaxis[iV].c_str());
      hist[iV][iT]->GetXaxis()->SetRangeUser(min[iV],max[iV]);
      //hist[iV][iT]->Sumw2();
      lname.str("");
      lname << vars[iV] << ">>p_" << vars[iV] << "_" << label[iT];
      tree[iT]->Draw(lname.str().c_str(),selection.str().c_str(),"");
      if (iV==6 && iT==0) {
	lname.str("");
	lname << vars[7] << ":" << vars[6] << ">>weight2D";// << label[iT];
	tree[iT]->Draw(lname.str().c_str(),selection.str().c_str(),"");
      }

      hist[iV][iT]->SetLineColor(iT+1);
      if (iT==0) {
	hist[iV][iT]->SetMarkerColor(iT+1);
	hist[iV][iT]->SetMarkerStyle(2);
      }


      //subtract MC
      lname.str("");
      lname << "dataSubtr_" << vars[iV] << "_" << label[iT];
      histDataSubtr[iV][iT] = (TH1F*)hist[iV][iT]->Clone(lname.str().c_str());
      //if (iT==0) 
	histDataSubtr[iV][iT]->Add(histMC[iV][iT],-1);
    }//loop on trees
    

    histsum[iV] = (TH1F*)hist[iV][1]->Clone("histsum");
    //histsum[iV]->Sumw2();
    histsum[iV]->Add(hist[iV][2]);

    TH1F *histTmp = (TH1F*)histDataSubtr[iV][0]->Clone("histTmp");
    //histTmp->Sumw2();
    histTmp->Divide(histsum[iV]);
    grWeight[iV] = new TGraphErrors(histTmp);

    histsumSubtr[iV] = (TH1F*)histDataSubtr[iV][1]->Clone("histsumSubtr");
    //histsumSubtr[iV]->Sumw2();
    histsumSubtr[iV]->Add(histDataSubtr[iV][2]);
    histsumSubtr[iV]->SetLineColor(4);
    histsumSubtr[iV]->SetFillColor(4);
    histsumSubtr[iV]->SetFillStyle(3004);


    myc[iV]->Clear();
    myc[iV]->Divide(2,1);
    gStyle->SetOptStat(0);

    double lmaxY = 0;

    //draw absolute scale
    myc[iV]->cd(1);
    if (hist[iV][0]->GetMaximum() > lmaxY) lmaxY =hist[iV][0]->GetMaximum();
    histsumSubtr[iV]->Scale(histDataSubtr[iV][0]->Integral()/histsumSubtr[iV]->Integral());
    histsumSubtr[iV]->Add(histMC[iV][0]);
    if (histsumSubtr[iV]->GetMaximum() > lmaxY) lmaxY =histsumSubtr[iV]->GetMaximum();
    if (histMC[iV][0]->GetMaximum() > lmaxY) lmaxY =histMC[iV][0]->GetMaximum();
    //if (histQCD[iV]->GetMaximum() > lmaxY) lmaxY =histQCD[iV]->GetMaximum();
    //if (histDataCheck[iV]->GetMaximum() > lmaxY) lmaxY =histDataCheck[iV]->GetMaximum();

    TLegend *leg = new TLegend(0.75,0.75,0.99,0.99);
    leg->SetFillColor(10);
    leg->AddEntry(hist[iV][0],"PARKED","P");
    hist[iV][0]->SetMaximum(lmaxY*1.1);
    hist[iV][0]->Draw("PE");
    histsumSubtr[iV]->Draw("histsame");
    histMC[iV][0]->Draw("histsame");
    //histQCD[iV]->Draw("histsame");
    //histDataCheck[iV]->Draw("PEsame");
    leg->AddEntry(histsumSubtr[iV],"j1j3+j2j3","F");
    leg->AddEntry(histMC[iV][0],"V+Top+VV","F");
    //leg->AddEntry(histDataCheck[iV],"PARKED","F");
    leg->Draw("same");

    //draw normalised
    myc[iV]->cd(2);
    lmaxY = 0;
    for (unsigned iT(0); iT<nTrees; ++iT){
      histDataSubtr[iV][iT]->GetYaxis()->SetTitle("arb. unit");
      histDataSubtr[iV][iT]->Scale(1./histDataSubtr[iV][iT]->Integral());
      if (histDataSubtr[iV][iT]->GetMaximum() > lmaxY) lmaxY =histDataSubtr[iV][iT]->GetMaximum();
    }
    //histsum[iV]->Scale(1./histsum[iV]->Integral());
    //histMC[iV]->Scale(1./histMC[iV]->Integral());
    histQCD[iV]->Scale(1./histQCD[iV]->Integral());
    //histDataCheck[iV]->Scale(1./histDataCheck[iV]->Integral());
    if (histQCD[iV]->GetMaximum() > lmaxY) lmaxY =histQCD[iV]->GetMaximum();
    //if (histDataCheck[iV]->GetMaximum() > lmaxY) lmaxY =histDataCheck[iV]->GetMaximum();

    TLegend *leg2 = new TLegend(0.75,0.75,0.99,0.99);
    leg2->SetFillColor(10);
    leg2->AddEntry(histDataSubtr[iV][0],"Data-MC","P");
    histDataSubtr[iV][0]->SetMaximum(lmaxY*1.1);
    histDataSubtr[iV][0]->Draw("PE");
    for (unsigned iT(1); iT<nTrees; ++iT){//loop on trees
      histDataSubtr[iV][iT]->SetMaximum(lmaxY*1.1);
      histDataSubtr[iV][iT]->Draw("histsame");
      leg2->AddEntry(histDataSubtr[iV][iT],label[iT].c_str(),"L");
    }
    //histsum[iV]->Draw("histsame");
    histQCD[iV]->Draw("histsame");
    //histDataCheck[iV]->Draw("Lsame");
    //leg2->AddEntry(histsum[iV],"j1j3+j2j3","F");
    leg2->AddEntry(histQCD[iV],"MC VBFQCD","F");
    //leg2->AddEntry(histDataCheck[iV],"j1j2 data","F");
    leg2->Draw("same");

    myc[iV]->Update();

    lname.str("");
    lname << "PLOTS/" << type << "_" << vars[iV] << ".pdf";
    myc[iV]->Print(lname.str().c_str());
    

  }//loop on vars


  TCanvas *mycpT1 = new TCanvas("pt1","pt1",1);
  mycpT1->cd();
  grWeight[varIdx]->SetTitle("");
  grWeight[varIdx]->GetXaxis()->SetTitle(vars[varIdx].c_str());
  grWeight[varIdx]->GetYaxis()->SetTitle("data/QCD");
  grWeight[varIdx]->Draw("APE");

  for (unsigned iP(0); iP<grWeight[varIdx]->GetN();++iP){
    double x,y;
    grWeight[varIdx]->GetPoint(iP,x,y);
    if (iP==0) std::cout << "if"; 
    else std::cout << "else if"; 
    std::cout << " (aVar>" << x-grWeight[varIdx]->GetErrorX(iP) << " &&  aVar <=" 
	      << x+grWeight[varIdx]->GetErrorX(iP) << ") return "
	      << y << ";//// " 
	      << grWeight[varIdx]->GetErrorY(iP)
	      << std::endl;
  }
  mycpT1->Print("qcd_reweighting_with_dphijj.pdf");

  TCanvas *mycpT12 = new TCanvas("pt1pt2","pt1pt2",1);
  mycpT12->cd();
  weight2D->Draw("colz");

  return 0;
}
