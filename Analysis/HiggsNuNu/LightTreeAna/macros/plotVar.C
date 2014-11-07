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

double deltaPhi(double jeta_phi, double jetb_phi){
  double dphi = jeta_phi-jetb_phi;
  if ( dphi > M_PI ) {
    dphi -= 2.0*M_PI;
  } else if ( dphi <= -M_PI ) {
    dphi += 2.0*M_PI;
  }
  return dphi;
}

int plotVar() {

  const unsigned nTrees = 1;

  //std::string fileName = "../../output_lighttree/Wjets_enu.root";
  std::string fileName[nTrees];
  //fileName[0] = "../../output_lighttree/VBFPARKED.root";
  //fileName[1] = "../../output_lighttree/QCDPARKED.root";
  //fileName[2] = "../../output_lighttree/QCDPARKED_0d5.root";
  fileName[0] = "../../output_lighttree/VBFQCD.root";
  //std::string fileName = "../../output_lighttree/MC_Powheg-Htoinv-mH125.root";

  //std::string type = "Wjets_enu";
  //std::string type = "VBFPARKED";
  std::string type = "VBFQCD";
  //std::string type = "VBFH125";

  bool isMC = true;
  //bool isMC = false;

  /*TFile *data = TFile::Open(fileName.c_str(), "update");
  if (!data) {
    std::cout << " Input file " << fileName << " not found. Exiting..." << std::endl;
    return 1;
  }

  data->cd();

  const unsigned nTrees = 3;
  TTree *tree[nTrees];
  tree[0] = (TTree*)gDirectory->Get("LightTree");
  tree[1] = (TTree*)gDirectory->Get("LightTreeJ1J3");
  tree[2] = (TTree*)gDirectory->Get("LightTreeQCD");

  */
  std::string label[nTrees] = {"VBF-QCD"};//data","qcd_1","qcd_0.5"};
  TFile *data[nTrees];
  TTree *tree[nTrees];

  for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
    int is_dupl = 0;
    data[iT] = TFile::Open(fileName[iT].c_str(), "update");
    if (!data[iT]) {
      std::cout << " Input file " << fileName[iT] << " not found. Exiting..." << std::endl;
      return 1;
    }
    
    data[iT]->cd();
    tree[iT] = (TTree*)gDirectory->Get("LightTree");
  }

  
  unsigned run;
  unsigned lumi;
  unsigned event;
  double passtrigger;
  double passparkedtrigger1;
  double passparkedtrigger2;
  double l1met;
  double jet1_pt;
  double jet2_pt;
  double jet3_pt;
  double dijet_M;
  double mindphi;
  double allmindphi;
  double jet1_eta;
  double jet2_eta;
  double jet3_eta;
  double jet1_phi;
  double jet2_phi;
  double jet3_phi;
  double metnomuons;
  double metnomu_x;
  double metnomu_y;
  unsigned n_jets_30;
  unsigned n_jets_15;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;
  float dRjajc;
  float dRjbjc;
  float dphiMETjc;

  /*
  std::set<Event> evtSet;
  std::pair<std::set<Event>::iterator,bool> isInserted;
  unsigned duplicateJ1J2 = 0;
  unsigned duplicateJ1J3 = 0;
  unsigned duplicateJ2J3 = 0;
  unsigned passJ1J2 = 0;
  unsigned passJ1J3 = 0;
  unsigned passJ2J3 = 0;
  */
  int nEntries[nTrees];
  
  for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
    int is_dupl = 0;
    data[iT] = TFile::Open(fileName[iT].c_str(), "update");
    if (!data[iT]) {
      std::cout << " Input file " << fileName[iT] << " not found. Exiting..." << std::endl;
      return 1;
    }
    
    data[iT]->cd();
    tree[iT] = (TTree*)gDirectory->Get("LightTree");

    
    //TBranch *brDupl = tree[iT]->Branch("is_dupl", &is_dupl, "is_dupl/I");
    TBranch *brdRjajc = tree[iT]->Branch("dRjajc", &dRjajc, "dRjajc/F");
    TBranch *brdRjbjc = tree[iT]->Branch("dRjbjc", &dRjbjc, "dRjbjc/F");
    TBranch *brdphiMETjc = tree[iT]->Branch("dphiMETjc", &dphiMETjc, "dphiMETjc/F");


    nEntries[iT] = tree[iT]->GetEntries();
    tree[iT]->SetBranchAddress("run",&run);
    tree[iT]->SetBranchAddress("lumi",&lumi);
    tree[iT]->SetBranchAddress("event",&event);
    tree[iT]->SetBranchAddress("passtrigger",&passtrigger);
    tree[iT]->SetBranchAddress("passparkedtrigger1",&passparkedtrigger1);
    tree[iT]->SetBranchAddress("passparkedtrigger2",&passparkedtrigger2);
    tree[iT]->SetBranchAddress("l1met",&l1met);
    tree[iT]->SetBranchAddress("jet1_pt",&jet1_pt);
    tree[iT]->SetBranchAddress("jet2_pt",&jet2_pt);
    tree[iT]->SetBranchAddress("jet3_pt",&jet3_pt);
    tree[iT]->SetBranchAddress("n_jets_30",&n_jets_30);
    tree[iT]->SetBranchAddress("n_jets_15",&n_jets_15);
    tree[iT]->SetBranchAddress("jetmetnomu_mindphi",&mindphi);
    tree[iT]->SetBranchAddress("alljetsmetnomu_mindphi",&allmindphi);
    tree[iT]->SetBranchAddress("metnomu_x",&metnomu_x);
    tree[iT]->SetBranchAddress("metnomu_y",&metnomu_y);
    tree[iT]->SetBranchAddress("dijet_M",&dijet_M);
    tree[iT]->SetBranchAddress("jet1_eta",&jet1_eta);
    tree[iT]->SetBranchAddress("jet2_eta",&jet2_eta);
    tree[iT]->SetBranchAddress("jet3_eta",&jet3_eta);
    tree[iT]->SetBranchAddress("jet1_phi",&jet1_phi);
    tree[iT]->SetBranchAddress("jet2_phi",&jet2_phi);
    tree[iT]->SetBranchAddress("jet3_phi",&jet3_phi);
    tree[iT]->SetBranchAddress("metnomuons",&metnomuons);
    tree[iT]->SetBranchAddress("nvetomuons",&nvetomuons);
    tree[iT]->SetBranchAddress("nselmuons",&nselmuons);
    tree[iT]->SetBranchAddress("nvetoelectrons",&nvetoelectrons);
    tree[iT]->SetBranchAddress("nselelectrons",&nselelectrons);

    std::cout << " Jet pair " << label[iT] << " has " << nEntries[iT] << " entries in tree." << std::endl;
    
    //std::cout << " Dupl set size: " << evtSet.size() << std::endl;

    for (int iE(0); iE<nEntries[iT]; ++iE){//loop on entries
      tree[iT]->GetEntry(iE);
      //is_dupl = 0;
      dRjajc = -1;
      dRjbjc = -1;
      dphiMETjc = -1;
      if (n_jets_30>2) {
	double deta = jet1_eta-jet3_eta;
	double dphi = deltaPhi(jet1_phi,jet3_phi); 
	dRjajc = sqrt(pow(deta,2)+pow(dphi,2));
	deta = jet2_eta-jet3_eta;
	dphi = deltaPhi(jet2_phi,jet3_phi); 
	dRjbjc = sqrt(pow(deta,2)+pow(dphi,2));
	
	dphiMETjc = deltaPhi(jet3_phi,atan(metnomu_y/metnomu_x));
      }
      /*
      //bool passtrig = (!isMC && passtrigger==1 && l1met>40) || isMC;//MET data
      bool passtrig = (!isMC && ((((run>=190456)&&(run<=193621))&&passtrigger==1)||(((run>=193833)&&(run<=196531))&&passparkedtrigger1==1)||(((run>=203777)&&(run<=208686))&&passparkedtrigger2==1)) && l1met>40) || isMC;//parked

      bool passnj = true;//(iT==0) || (iT>0 && n_jets>2);

      bool passdphi = mindphi>1 && allmindphi<1;

      if (passtrig && passnj && passdphi && jet1_eta*jet2_eta < 0 && dijet_M>=800 && jet1_pt > 50 && jet2_pt > 40  && metnomuons>90 && nvetomuons==0 && nvetoelectrons==0){//pass sel
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
	  else if (iT==1){
	    passJ1J3++;
	    evtSet.erase(lEvt);
	  }
	  else 
	    passJ2J3++;
	}
      
      }//pass sel
      */
      //brDupl->Fill();
      brdRjajc->Fill();
      brdRjbjc->Fill();
      brdphiMETjc->Fill();
    }//loop on entries

    //std::cout << " Dupl set size after: " << evtSet.size() << std::endl;
  
  }//loop on trees

  /*
  std::cout << " Duplicated/passing events: " << std::endl
	    << label[0] << " & " << duplicateJ1J2+passJ1J2 << " & " << duplicateJ1J2 << " & " << passJ1J2 << "\\\\" << std::endl
	    << label[1] << " & " << duplicateJ1J3+passJ1J3 << " & " << duplicateJ1J3 << " & " << passJ1J3 << "\\\\"<< std::endl
	    << label[2] << " & " << duplicateJ2J3+passJ2J3 << " & " << duplicateJ2J3 << " & " << passJ2J3 << "\\\\"<< std::endl
    ;
    */

  const unsigned nVars = 26;
  std::string vars[nVars] = {
    "dijet_M",
    "metnomuons",	     
    "dijet_deta",
    "dijet_dphi",
    "metnomu_significance",
    "jetmetnomu_mindphi",
    "alljetsmetnomu_mindphi",
    "jet1_pt",
    "jet2_pt",
    "jet3_pt",
    "jet4_pt",
    "jet1_csv",
    "jet2_csv",
    "jet3_csv",
    "jet4_csv",
    "dijetmetnomu_scalarSum_pt",
    "dijetmetnomu_vectorialSum_pt",
    "dijetmetnomu_ptfraction",
    "n_jets_30",
    "jet1_pdgid",
    "jet2_pdgid",
    "jet3_pdgid",   
    "jet4_pdgid",
    "dRjajc",
    "dRjbjc",
    "dphiMETjc"
  };

  std::string xaxis[nVars] = {
    ";M_{jj} (GeV)",
    ";METnoMu (GeV)",
    ";#Delta#eta_{jj}",
    ";#Delta#phi_{jj}",
    ";METnoMu/#sigma(METnoMu)",
    ";min #Delta#phi(j_{1,2},METnoMu)",
    ";min #Delta#phi(j,METnoMu)",
    ";p_{T}^{ja} (GeV)",
    ";p_{T}^{jb} (GeV)",
    ";p_{T}^{jc} (GeV)",
    ";p_{T}^{j4} (GeV)",
    ";CSV (jet a)",
    ";CSV (jet b)",
    ";CSV (jet c)",
    ";CSV (jet 4)",
    ";p_{T}^{jeta}+p_{T}^{jetb}+METnoMu",
    ";p_{T}(#vec{ja}+#vec{jb}+#vec{METnoMu})",
    ";p_{T}^{dijet}/(p_{T}^{dijet}+METnoMu)",
    ";njets (30 GeV)",
    ";PDGID (jet a)",
    ";PDGID (jet b)",
    ";PDGID (jet c)",
    ";PDGID (jet 4)",
    ";#Delta R(ja,jc)",
    ";#Delta R(jb,jc)",
    ";#Delta#phi(jc,METnoMu)"
  };

unsigned nbins[nVars] = {75,100,44,30,35,30,30,50,50,50,50,50,50,50,50,100,50,50,10,34,34,34,34,50,50,50};
float min[nVars] = {0,0,3.6,0,3,0,0,30,30,15,15,0,0,0,0,0,0,0,0,-12,-12,-12,-12,0,0,0};
float max[nVars] = {3000,500,8,3.1416,10,3.1416,3.1416,300,300,300,300,1,1,1,1,1000,400,1,10,22,22,22,22,10,10,3.1416};  

  TH1F *hist[nVars][nTrees];
  //TH1F *histsum[nVars];

  std::ostringstream lname;
  TCanvas *myc[nVars+4];
  /*for (unsigned iV(0); iV<nVars+4; ++iV){//loop on variables
    lname.str("");
    if (iV<nVars) lname << "myc_" << vars[iV];
    else lname << "myc_" << iV-nVars;
    myc[iV] = new TCanvas(lname.str().c_str(),
			  lname.str().c_str(),
			  1000,1000);
			  }*/

  for (unsigned iV(0); iV<nVars; ++iV){//loop on variables
    lname.str("");
    lname << "myc_" << vars[iV];
    myc[iV] = new TCanvas(lname.str().c_str(),
			  lname.str().c_str(),
			  1000,1000);
    myc[iV]->cd();
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      
      //std::string selection = "jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0";
      //std::string selection = "is_dupl==0 && jet1_pt > 50 && dijet_M > 600 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0"; 
      std::ostringstream selection;
      //bool passtrig = (!isMC && passtrigger==1 && l1met>40) || isMC;//MET data
      std::string passtrig = "((run>=190456 && run<=193621 && passtrigger==1) || (run>=193833 && run<=196531 &&passparkedtrigger1==1) || (run>=203777 && run<=208686 && passparkedtrigger2==1)) && l1met>40";

      //std::string signal = "metnomuons<130 && alljetsmetnomu_mindphi>1.0 && alljetsmetnomu_mindphi<2.0 && metnomu_significance < 4 && dijet_dphi > 1.5";
      std::string signal = "jetmetnomu_mindphi>1.5 && metnomu_significance > 3";
      if (iT==0) selection << signal << " && " ;
      //if (iT>0) selection << " n_jets_30>2 && ";
      //selection << "is_dupl==0 && ";
      if (!isMC) selection << passtrig <<" && ";
      selection << "jet1_eta*jet2_eta<0 && dijet_deta>3.6 && jet1_pt>50 && jet2_pt>40 && metnomuons>90 && dijet_M > 800 && nvetomuons==0 && nvetoelectrons==0";
      //if (iT>0) selection << " && alljetsmet_mindphi<1.0 && jetmetnomu_mindphi>1.";
      
      lname.str("");
      lname << "p_" << vars[iV] << "_" << label[iT];
      hist[iV][iT] = new TH1F(lname.str().c_str(),xaxis[iV].c_str(),nbins[iV],min[iV],max[iV]);
      hist[iV][iT]->Sumw2();
      lname.str("");
      lname << vars[iV] << ">>p_" << vars[iV] << "_" << label[iT];
      tree[iT]->Draw(lname.str().c_str(),selection.str().c_str(),"");

      hist[iV][iT]->SetLineColor(6);//iT+1);
      hist[iV][iT]->SetFillColor(6);//iT+1);
      //if (iT==0) {
      //hist[iV][iT]->SetMarkerColor(iT+1);
      //hist[iV][iT]->SetMarkerStyle(2);
      //}
    }//loop on trees
    
    //histsum[iV] = (TH1F*)hist[iV][1]->Clone("histsum");
    //histsum[iV]->Sumw2();
    //histsum[iV]->Add(hist[iV][2]);
    //histsum[iV]->SetLineColor(4);
    //histsum[iV]->SetFillColor(4);
    //histsum[iV]->SetFillStyle(3004);

    myc[iV]->Clear();
    myc[iV]->cd();
    gStyle->SetOptStat(0);

    double lmaxY = 0;
    for (unsigned iT(0); iT<nTrees; ++iT){
      hist[iV][iT]->GetYaxis()->SetTitle("Events");//arb. unit");
      //hist[iV][iT]->Scale(1./hist[iV][iT]->Integral());
      if (hist[iV][iT]->GetMaximum() > lmaxY) lmaxY =hist[iV][iT]->GetMaximum();
    }
    //histsum[iV]->Scale(1./histsum[iV]->Integral());
    //if (histsum[iV]->GetMaximum() > lmaxY) lmaxY =histsum[iV]->GetMaximum();

    TLegend *leg = new TLegend(0.62,0.75,0.86,0.99);
    leg->SetFillColor(10);
    for (unsigned iT(0); iT<nTrees; ++iT){//loop on trees
      hist[iV][iT]->SetMaximum(lmaxY*1.1);
      //hist[iV][iT]->Draw(iT==0?"PE":"histsame");
      hist[iV][iT]->Draw("hist");
      leg->AddEntry(hist[iV][iT],label[iT].c_str(),"L");
    }
    //histsum[iV]->Draw("histsame");
    //leg->AddEntry(histsum[iV],"j1j3+j2j3","F");
    //leg->Draw("same");

    myc[iV]->Update();

    lname.str("");
    lname << "PLOTS/" << type << "_" << vars[iV] << ".pdf";
    myc[iV]->Print(lname.str().c_str());
    

  }//loop on vars

  /*
  myc[nVars]->cd();
  hist[6][0]->SetMaximum(hist[6][0]->GetMaximum()*1.1);
  hist[6][0]->Draw("PE");
  hist[6][1]->Draw("histsame");
  hist[8][2]->Draw("histsame");
  //histsum[6]->Add(hist[6][1],hist[8][2]);
  //histsum[6]->Scale(1./histsum[6]->Integral());
  //histsum[6]->Draw("histsame");
  myc[nVars]->Update();
  lname.str("");
  lname << "PLOTS/" << type << "_pTlead.pdf";
  myc[nVars]->Print(lname.str().c_str());

  myc[nVars+1]->cd();
  hist[7][0]->SetMaximum(0.1);
  hist[7][0]->Draw("PE");
  hist[8][1]->Draw("histsame");
  hist[6][2]->Draw("histsame");
  //histsum[7]->Add(hist[8][1],hist[6][2]);
  //histsum[7]->Scale(1./histsum[7]->Integral());
  //histsum[7]->Draw("histsame");
  myc[nVars+1]->Update();
  lname.str("");
  lname << "PLOTS/" << type << "_pTsublead.pdf";
  myc[nVars+1]->Print(lname.str().c_str());

  myc[nVars+2]->cd();
  //histsum[8]->Add(hist[7][1],hist[7][2]);
  //histsum[8]->Scale(1./histsum[8]->Integral());
  //hist[8][0]->SetMaximum(histsum[8]->GetMaximum()*1.1);
  hist[8][0]->Draw("PE");
  hist[7][1]->Draw("histsame");
  hist[7][2]->Draw("histsame");
  //histsum[8]->Draw("histsame");
  myc[nVars+2]->Update();
  lname.str("");
  lname << "PLOTS/" << type << "_pTthird.pdf";
  myc[nVars+2]->Print(lname.str().c_str());

  myc[nVars+3]->cd();
  //histsum[9]->Add(hist[9][1],hist[9][2]);
  //histsum[9]->Scale(1./histsum[9]->Integral());
  //hist[9][0]->SetMaximum(histsum[9]->GetMaximum()*1.1);
  hist[9][0]->Draw("PE");
  hist[9][1]->Draw("histsame");
  hist[9][2]->Draw("histsame");
  //histsum[9]->Draw("histsame");
  myc[nVars+3]->Update();
  lname.str("");
  lname << "PLOTS/" << type << "_pTfourth.pdf";
  myc[nVars+3]->Print(lname.str().c_str());
  */
    


  return 0;
}
