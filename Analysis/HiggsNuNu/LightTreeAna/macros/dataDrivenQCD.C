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

int dataDrivenQCD() {

  std::string fileName = "../../output_lighttree/VBFPARKED.root";//_MET-2012A-22Jan2013-v1.root";
  //std::string fileName = "../../output/MC_Powheg-Htoinv-mH125.root";

  std::string type = "DataMC_PARKED";
  //std::string type = "VBFH125";

  TFile *data = TFile::Open(fileName.c_str(), "update");
  if (!data) {
    std::cout << " Input file " << fileName << " not found. Exiting..." << std::endl;
    return 1;
  }

  data->cd();

  const unsigned nTrees = 3;
  std::string label[nTrees] = {"j1j2","j1j3","j2j3"};

  TTree *tree[nTrees];
  tree[0] = (TTree*)gDirectory->Get("LightTree");
  tree[1] = (TTree*)gDirectory->Get("LightTreeJ1J3");
  tree[2] = (TTree*)gDirectory->Get("LightTreeJ2J3");
  TFile *mcFile[nTrees];
  mcFile[0] = TFile::Open("../output/nunu.root");
  mcFile[1] = TFile::Open("../output/nunu_J1J3.root");
  mcFile[2] = TFile::Open("../output/nunu_J2J3.root");

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
  double jet1_eta;
  double jet2_eta;
  double metnomuons;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;
  unsigned n_jets;

  std::set<Event> evtSet;
  std::pair<std::set<Event>::iterator,bool> isInserted;
  unsigned duplicateJ1J2 = 0;
  unsigned duplicateJ1J3 = 0;
  unsigned duplicateJ2J3 = 0;
  unsigned passJ1J2 = 0;
  unsigned passJ1J3 = 0;
  unsigned passJ2J3 = 0;
  
  int nEntries[nTrees];
  
  for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
    int is_dupl = 0;
    TBranch *brDupl = tree[iT]->Branch("is_dupl", &is_dupl, "is_dupl/I");
    
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
      brDupl->Fill();
    }//loop on entries
  }//loop on trees


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
  TH1F *histMC[nVars][nTrees];
  TH1F *histDataCheck[nVars];
  TH1F *histDataSubtr[nVars][nTrees];
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
      
    //get MC shapes
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      
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

      std::ostringstream selection;
      std::string passtrig = "((run>=190456 && run<=193621 && passtrigger==1) || (run>=193833 && run<=196531 &&passparkedtrigger1==1) || (run>=203777 && run<=208686 && passparkedtrigger2==1)) && l1met>40";
      std::string passpT = "jet1_pt > 50";

      if (iT>0) selection << " n_jets>2 && ";
      selection << "is_dupl==0 && ";
      selection << passtrig <<" && ";
      if (iT==0) selection << passpT << " && ";
      selection << "jet1_eta*jet2_eta<0 && dijet_M>=600 && (jet1_pt>50 || jet3_pt>50) && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      
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

      hist[iV][iT]->SetLineColor(iT+1);
      if (iT==0) {
	hist[iV][iT]->SetMarkerColor(iT+1);
	hist[iV][iT]->SetMarkerStyle(2);
      }

      //subtract MC
      lname.str("");
      lname << "dataSubtr_" << vars[iV] << "_" << label[iT];
      histDataSubtr[iV][iT] = (TH1F*)hist[iV][iT]->Clone(lname.str().c_str());
      histDataSubtr[iV][iT]->Add(histMC[iV][iT],-1);
    }//loop on trees
    
    histsum[iV] = (TH1F*)histDataSubtr[iV][1]->Clone("histsum");
    //histsum[iV]->Sumw2();
    histsum[iV]->Add(histDataSubtr[iV][2]);
    histsum[iV]->SetLineColor(4);
    histsum[iV]->SetFillColor(4);
    histsum[iV]->SetFillStyle(3004);


    myc[iV]->Clear();
    myc[iV]->Divide(2,1);
    gStyle->SetOptStat(0);

    double lmaxY = 0;

    //draw absolute scale
    myc[iV]->cd(1);
    if (hist[iV][0]->GetMaximum() > lmaxY) lmaxY =hist[iV][0]->GetMaximum();
    histsum[iV]->Scale(histDataSubtr[iV][0]->Integral()/histsum[iV]->Integral());
    histsum[iV]->Add(histMC[iV][0]);
    if (histsum[iV]->GetMaximum() > lmaxY) lmaxY =histsum[iV]->GetMaximum();
    if (histMC[iV][0]->GetMaximum() > lmaxY) lmaxY =histMC[iV][0]->GetMaximum();
    //if (histQCD[iV]->GetMaximum() > lmaxY) lmaxY =histQCD[iV]->GetMaximum();
    //if (histDataCheck[iV]->GetMaximum() > lmaxY) lmaxY =histDataCheck[iV]->GetMaximum();

    TLegend *leg = new TLegend(0.75,0.75,0.99,0.99);
    leg->SetFillColor(10);
    leg->AddEntry(hist[iV][0],"PARKED","P");
    hist[iV][0]->SetMaximum(lmaxY*1.1);
    hist[iV][0]->Draw("PE");
    histsum[iV]->Draw("histsame");
    histMC[iV][0]->Draw("histsame");
    //histQCD[iV]->Draw("histsame");
    //histDataCheck[iV]->Draw("PEsame");
    leg->AddEntry(histsum[iV],"j1j3+j2j3","F");
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

  return 0;
}
