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

  std::string fileName = "../../output_lighttree/MET_MET-2012A-22Jan2013-v1.root";
  TFile *data = TFile::Open(fileName.c_str());
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
  double jet1_pt;
  double dijet_M;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;

  const unsigned nVars = 8;
  std::string vars[nVars] = {
    "dijet_M",
    "metnomuons",	     
    "dijet_deta",
    "dijet_dphi",
    "metnomu_significance",
    "jetmetnomu_mindphi",
    "jet1_pt",
    "jet2_pt"
  };

  std::string xaxis[nVars] = {
    ";M_{jj} (GeV)",
    ";METnoMu (GeV)",
    ";#Delta#eta_{jj}",
    ";#Delta#phi_{jj}",
    ";METnoMu/#sigma(METnoMu)",
    ";min #Delta#phi(j,METnoMu)",
    ";p_{T}^{j1} (GeV)",
    ";p_{T}^{j2} (GeV)"
  };

  unsigned nbins[nVars] = {150,100,44,60,35,30,100,100};
  float min[nVars] = {0,0,3.6,0,3,1.5,30,30};
  float max[nVars] = {3000,500,8,3.1416,10,3.1416,300,300};  

  int nEntries[nTrees];
  TH1F *hist[nVars][nTrees];

  std::set<Event> evtSet;
  std::pair<std::set<Event>::iterator,bool> isInserted;
  unsigned duplicateJ1J2 = 0;
  unsigned duplicateJ1J3 = 0;
  unsigned duplicateJ2J3 = 0;

  TCanvas *myc[nVars];
  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    myc[iV] = new TCanvas(vars[iV].c_str(),
			  vars[iV].c_str(),
			  1);
    
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      
      nEntries[iT] = tree[iT]->GetEntries();
      tree[iT]->SetBranchAddress("run",&run);
      tree[iT]->SetBranchAddress("lumi",&lumi);
      tree[iT]->SetBranchAddress("event",&event);
      tree[iT]->SetBranchAddress("jet1_pt",&jet1_pt);
      tree[iT]->SetBranchAddress("dijet_M",&dijet_M);
      tree[iT]->SetBranchAddress("nvetomuons",&nvetomuons);
      tree[iT]->SetBranchAddress("nselmuons",&nselmuons);
      tree[iT]->SetBranchAddress("nvetoelectrons",&nvetoelectrons);
      tree[iT]->SetBranchAddress("nselelectrons",&nselelectrons);

      std::cout << " Jet pair " << label[iT] << " has " << nEntries[iT] << " entries in tree." << std::endl;
      
      if (iV==0) {
	for (unsigned iE(0); iE<nEntries[iT]; ++iE){//loop on entries
	  tree[iT]->GetEntry(iE);
	  if (jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0){//pass sel
	    Event lEvt;
	    lEvt.run = run;
	    lEvt.event = event;
	    lEvt.lumi = lumi;
	    isInserted = evtSet.insert(lEvt);
	    if (!isInserted.second){
	      if (iT==0) duplicateJ1J2++;
	      else if (iT==1) duplicateJ1J3++;
	      else duplicateJ2J3++;
	    }
	  }//pass sel
	}//loop on entries
      }
      std::string selection = "jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      
      std::ostringstream lname,laxis;
      lname.str("");
      lname << "p_" << vars[iV] << "_" << label[iT];
      hist[iV][iT] = new TH1F(lname.str().c_str(),xaxis[iV].c_str(),nbins[iV],min[iV],max[iV]);
      lname.str("");
      lname << vars[iV] << ">>p_" << vars[iV] << "_" << label[iT];
      tree[iT]->Draw(lname.str().c_str(),selection.c_str(),"");

      hist[iV][iT]->SetLineColor(iT+1);
      
    }//loop on trees
    
    myc[iV]->cd();
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      hist[iV][iT]->DrawNormalized(iT==0?"":"same");
    }

  }//loop on vars

  std::cout << " Duplicated events: " << std::endl
	    << " -- J1J2 = " << duplicateJ1J2 << std::endl
	    << " -- J1J3 = " << duplicateJ1J3 << std::endl
	    << " -- J2J3 = " << duplicateJ2J3 << std::endl
    ;

  return 0;
}
