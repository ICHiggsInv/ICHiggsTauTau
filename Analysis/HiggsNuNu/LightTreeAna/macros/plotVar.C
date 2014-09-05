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

int plotVar() {

  std::string fileName = "../../output_lighttree/VBFQCD.root";
  //std::string fileName = "../../output_lighttree/MC_Powheg-Htoinv-mH125.root";

  std::string type = "VBFQCD";
  //std::string type = "VBFH125";

  bool isMC = true;

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

  unsigned run;
  unsigned lumi;
  unsigned event;
  double passtrigger;
  double jet1_pt;
  double dijet_M;
  double jet1_eta;
  double jet2_eta;
  double metnomuons;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;

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
    tree[iT]->SetBranchAddress("jet1_pt",&jet1_pt);
    tree[iT]->SetBranchAddress("dijet_M",&dijet_M);
    tree[iT]->SetBranchAddress("jet1_eta",&jet1_eta);
    tree[iT]->SetBranchAddress("jet2_eta",&jet2_eta);
    tree[iT]->SetBranchAddress("metnomuons",&metnomuons);
    tree[iT]->SetBranchAddress("nvetomuons",&nvetomuons);
    tree[iT]->SetBranchAddress("nselmuons",&nselmuons);
    tree[iT]->SetBranchAddress("nvetoelectrons",&nvetoelectrons);
    tree[iT]->SetBranchAddress("nselelectrons",&nselelectrons);

    std::cout << " Jet pair " << label[iT] << " has " << nEntries[iT] << " entries in tree." << std::endl;
    
    for (int iE(0); iE<nEntries[iT]; ++iE){//loop on entries
      tree[iT]->GetEntry(iE);
      is_dupl = 0;

      if (((!isMC && passtrigger==1) || isMC) && jet1_eta*jet2_eta<0 && dijet_M>=600 && jet1_pt>50 && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0){//pass sel
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

  const unsigned nVars = 19;
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
    "jet1_pdgid",
    "jet2_pdgid",
    "jet3_pdgid"    
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
    ";njets (30 GeV)",
    ";PDGID (jet1)",
    ";PDGID (jet2)",
    ";PDGID (jet3)"
  };

  unsigned nbins[nVars] = {75,100,44,30,35,30,50,50,50,50,50,50,100,50,50,10,44,44,44};
  float min[nVars] = {0,0,3.6,0,3,1.5,30,30,30,0,0,0,0,0,0,0,-22,-22,-22};
  float max[nVars] = {3000,500,8,3.1416,10,3.1416,300,300,300,1,1,1,1000,400,1,10,22,22,22};  

  TH1F *hist[nVars][nTrees];
  TH1F *histsum[nVars];

  TCanvas *myc[nVars];
  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    std::ostringstream lname;
    lname.str("");
    lname << "myc_" << vars[iV];
    myc[iV] = new TCanvas(lname.str().c_str(),
			  vars[iV].c_str(),
			  1000,1000);

    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      
      //std::string selection = "jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      //std::string selection = "is_dupl==0 && jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0"; 
      std::string selection;
      if (!isMC) selection = 
		   "is_dupl==0 && passtrigger==1 && jet1_eta*jet2_eta<0 && dijet_M>=600 && jet1_pt>50 && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      else selection = 
	     "is_dupl==0 && jet1_eta*jet2_eta<0 && dijet_M>=600 && jet1_pt>50 && metnomuons>60 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      
      lname.str("");
      lname << "p_" << vars[iV] << "_" << label[iT];
      hist[iV][iT] = new TH1F(lname.str().c_str(),xaxis[iV].c_str(),nbins[iV],min[iV],max[iV]);
      hist[iV][iT]->Sumw2();
      lname.str("");
      lname << vars[iV] << ">>p_" << vars[iV] << "_" << label[iT];
      tree[iT]->Draw(lname.str().c_str(),selection.c_str(),"");

      hist[iV][iT]->SetLineColor(iT+1);
      if (iT==0) {
	hist[iV][iT]->SetMarkerColor(iT+1);
	hist[iV][iT]->SetMarkerStyle(2);
      }
    }//loop on trees
    
    histsum[iV] = (TH1F*)hist[iV][1]->Clone("histsum");
    //histsum[iV]->Sumw2();
    histsum[iV]->Add(hist[iV][2]);
    histsum[iV]->SetLineColor(4);
    histsum[iV]->SetFillColor(4);
    histsum[iV]->SetFillStyle(3004);

    myc[iV]->Clear();
    myc[iV]->cd();
    gStyle->SetOptStat(0);

    double lmaxY = 0;
    for (unsigned iT(0); iT<nTrees; ++iT){
      hist[iV][iT]->GetYaxis()->SetTitle("arb. unit");
      hist[iV][iT]->Scale(1./hist[iV][iT]->Integral());
      if (hist[iV][iT]->GetMaximum() > lmaxY) lmaxY =hist[iV][iT]->GetMaximum();
    }
    histsum[iV]->Scale(1./histsum[iV]->Integral());
    if (histsum[iV]->GetMaximum() > lmaxY) lmaxY =histsum[iV]->GetMaximum();

    TLegend *leg = new TLegend(0.75,0.75,0.99,0.99);
    leg->SetFillColor(10);
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      hist[iV][iT]->SetMaximum(lmaxY*1.1);
      hist[iV][iT]->Draw(iT==0?"PE":"histsame");
      leg->AddEntry(hist[iV][iT],label[iT].c_str(),"L");
    }
    histsum[iV]->Draw("histsame");
    leg->AddEntry(histsum[iV],"j1j3+j2j3","F");
    leg->Draw("same");

    myc[iV]->Update();

    lname.str("");
    lname << "PLOTS/" << type << "_" << vars[iV] << ".pdf";
    myc[iV]->Print(lname.str().c_str());
    

  }//loop on vars

  return 0;
}
