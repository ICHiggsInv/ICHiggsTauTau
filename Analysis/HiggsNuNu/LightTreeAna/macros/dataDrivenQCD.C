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
  /*
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
  */
  //mindphi_j1j2met
if (aVar>0 &&  aVar <=0.05236) return 1.77429;//// 0.0285863
else if (aVar>0.05236 &&  aVar <=0.10472) return 1.81757;//// 0.0306516
else if (aVar>0.10472 &&  aVar <=0.15708) return 1.90692;//// 0.0365168
else if (aVar>0.15708 &&  aVar <=0.20944) return 1.97956;//// 0.0453416
else if (aVar>0.20944 &&  aVar <=0.2618) return 2.0237;//// 0.055573
else if (aVar>0.2618 &&  aVar <=0.31416) return 2.05141;//// 0.0680555
else if (aVar>0.31416 &&  aVar <=0.36652) return 2.08484;//// 0.0869034
else if (aVar>0.36652 &&  aVar <=0.41888) return 1.99372;//// 0.0981167
else if (aVar>0.41888 &&  aVar <=0.47124) return 1.92691;//// 0.12008
else if (aVar>0.47124 &&  aVar <=0.5236) return 1.6321;//// 0.108482
else if (aVar>0.5236 &&  aVar <=0.57596) return 1.50407;//// 0.109121
else if (aVar>0.57596 &&  aVar <=0.62832) return 1.35802;//// 0.110702
else if (aVar>0.62832 &&  aVar <=0.68068) return 1.27123;//// 0.113009
else if (aVar>0.68068 &&  aVar <=0.73304) return 1.01642;//// 0.0944162
else if (aVar>0.73304 &&  aVar <=0.7854) return 1.07142;//// 0.123253
else if (aVar>0.7854 &&  aVar <=0.83776) return 0.813294;//// 0.0869153
else if (aVar>0.83776 &&  aVar <=0.89012) return 0.649379;//// 0.0676864
else if (aVar>0.89012 &&  aVar <=0.94248) return 0.689959;//// 0.0842654
else if (aVar>0.94248 &&  aVar <=0.99484) return 0.650578;//// 0.0801251
else if (aVar>0.99484 &&  aVar <=1.0472) return 0.602928;//// 0.0771214
else if (aVar>1.0472 &&  aVar <=1.09956) return 0.529747;//// 0.0663426
else if (aVar>1.09956 &&  aVar <=1.15192) return 0.478993;//// 0.0700237
else if (aVar>1.15192 &&  aVar <=1.20428) return 0.464871;//// 0.0649029
else if (aVar>1.20428 &&  aVar <=1.25664) return 0.496775;//// 0.0777062
else if (aVar>1.25664 &&  aVar <=1.309) return 0.420256;//// 0.06702
else if (aVar>1.309 &&  aVar <=1.36136) return 0.41736;//// 0.0757746
else if (aVar>1.36136 &&  aVar <=1.41372) return 0.318631;//// 0.0573408
else if (aVar>1.41372 &&  aVar <=1.46608) return 0.416459;//// 0.0671769
else if (aVar>1.46608 &&  aVar <=1.51844) return 0.446135;//// 0.0781064
else if (aVar>1.51844 &&  aVar <=1.5708) return 0.432387;//// 0.0760104
else if (aVar>1.5708 &&  aVar <=1.62316) return 0.341891;//// 0.0708654
else if (aVar>1.62316 &&  aVar <=1.67552) return 0.411269;//// 0.0828041
else if (aVar>1.67552 &&  aVar <=1.72788) return 0.538317;//// 0.085166
else if (aVar>1.72788 &&  aVar <=1.78024) return 0.395874;//// 0.0821543
else if (aVar>1.78024 &&  aVar <=1.8326) return 0.500289;//// 0.105725
else if (aVar>1.8326 &&  aVar <=1.88496) return 0.595132;//// 0.10996
else if (aVar>1.88496 &&  aVar <=1.93732) return 0.492405;//// 0.0935388
else if (aVar>1.93732 &&  aVar <=1.98968) return 0.427878;//// 0.0887358
else if (aVar>1.98968 &&  aVar <=2.04204) return 0.413447;//// 0.0882058
else if (aVar>2.04204 &&  aVar <=2.0944) return 0.434929;//// 0.0845145
else if (aVar>2.0944 &&  aVar <=2.14676) return 0.368647;//// 0.0867877
else if (aVar>2.14676 &&  aVar <=2.19912) return 0.442155;//// 0.0774579
else if (aVar>2.19912 &&  aVar <=2.25148) return 0.471863;//// 0.0837436
else if (aVar>2.25148 &&  aVar <=2.30384) return 0.324949;//// 0.0654281
else if (aVar>2.30384 &&  aVar <=2.3562) return 0.464887;//// 0.081463
else if (aVar>2.3562 &&  aVar <=2.40856) return 0.339325;//// 0.0725929
else if (aVar>2.40856 &&  aVar <=2.46092) return 0.314148;//// 0.0705466
else if (aVar>2.46092 &&  aVar <=2.51328) return 0.378439;//// 0.0671156
else if (aVar>2.51328 &&  aVar <=2.56564) return 0.282832;//// 0.062221
else if (aVar>2.56564 &&  aVar <=2.618) return 0.175896;//// 0.0552767
else if (aVar>2.618 &&  aVar <=2.67036) return 0.1828;//// 0.0521285
else if (aVar>2.67036 &&  aVar <=2.72272) return 0.199943;//// 0.0522814
else if (aVar>2.72272 &&  aVar <=2.77508) return 0.213227;//// 0.0517293
else if (aVar>2.77508 &&  aVar <=2.82744) return 0.2273;//// 0.0489419
else if (aVar>2.82744 &&  aVar <=2.8798) return 0.156259;//// 0.0498699
else if (aVar>2.8798 &&  aVar <=2.93216) return 0.098822;//// 0.0449051
else if (aVar>2.93216 &&  aVar <=2.98452) return 0.200001;//// 0.0516363
else if (aVar>2.98452 &&  aVar <=3.03688) return 0.215024;//// 0.0557044
else if (aVar>3.03688 &&  aVar <=3.08924) return 0.239018;//// 0.0746044
else if (aVar>3.08924 &&  aVar <=3.1416) return 0.162946;//// 0.106139

 return 1;

};

int dataDrivenQCD() {

  std::string fileName = "../../output_lighttree/VBFPARKED.root";//_MET-2012A-22Jan2013-v1.root";
  //std::string fileName = "../../output/MC_Powheg-Htoinv-mH125.root";

  bool doReweighting = false;
  std::string type = "DataMC_PARKED";//_rwmindphi";
  //std::string type = "VBFH125";
  unsigned varIdx = 5;

  const unsigned nTrees = 3;
  std::string label[nTrees] = {"j1j2","j1j3","qcd"};

  TFile *mcFile[nTrees];
  mcFile[0] = TFile::Open("../output/nunu.root");
  mcFile[1] = TFile::Open("../output/nunu_J1J3.root");
  mcFile[2] = TFile::Open("../output/nunu_QCD.root");

  TFile *data = TFile::Open(fileName.c_str(), "update");
  if (!data) {
    std::cout << " Input file " << fileName << " not found. Exiting..." << std::endl;
    return 1;
  }

  data->cd();


  TTree *tree[nTrees];
  tree[0] = (TTree*)gDirectory->Get("LightTree");
  tree[1] = (TTree*)gDirectory->Get("LightTreeJ1J3");
  tree[2] = (TTree*)gDirectory->Get("LightTreeQCD");

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
  double dijet_dphi;
  double mindphi;
  double allmindphi;
  double jet1_eta;
  double jet2_eta;
  double metnomuons;
  int nvetomuons;
  int nselmuons;
  int nvetoelectrons;
  int nselelectrons;
  unsigned n_jets_30;
  unsigned n_jets_15;
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
    tree[iT]->SetBranchAddress("jet2_pt",&jet2_pt);
    tree[iT]->SetBranchAddress("jet3_pt",&jet3_pt);
    tree[iT]->SetBranchAddress("dijet_M",&dijet_M);
    tree[iT]->SetBranchAddress("dijet_dphi",&dijet_dphi);
    tree[iT]->SetBranchAddress("jetmetnomu_mindphi",&mindphi);
    tree[iT]->SetBranchAddress("alljetsmetnomu_mindphi",&allmindphi);
    tree[iT]->SetBranchAddress("jet1_eta",&jet1_eta);
    tree[iT]->SetBranchAddress("jet2_eta",&jet2_eta);
    tree[iT]->SetBranchAddress("metnomuons",&metnomuons);
    tree[iT]->SetBranchAddress("nvetomuons",&nvetomuons);
    tree[iT]->SetBranchAddress("nselmuons",&nselmuons);
    tree[iT]->SetBranchAddress("nvetoelectrons",&nvetoelectrons);
    tree[iT]->SetBranchAddress("nselelectrons",&nselelectrons);
    tree[iT]->SetBranchAddress("n_jets_30",&n_jets_30);
    tree[iT]->SetBranchAddress("n_jets_15",&n_jets_15);

    std::cout << " Jet pair " << label[iT] << " has " << nEntries[iT] << " entries in tree." << std::endl;
    
    for (int iE(0); iE<nEntries[iT]; ++iE){//loop on entries
      tree[iT]->GetEntry(iE);
      is_dupl = 0;
      bool passtrig = ((run>=190456 && run<=193621 &&passtrigger==1) || (run>=193833 && run<=196531 && passparkedtrigger1==1) ||(run>=203777 && run<=208686 && passparkedtrigger2==1)) && l1met>40;//parked
      //bool passpT = (iT==0 && jet1_pt > 50) || iT>0;
      bool passpT = jet1_pt > 50 && jet2_pt > 40;
      bool passnj = (iT==0) || (iT>0 && n_jets_30>2);
      qcdW = qcdWeight(mindphi);

      if (passtrig && passpT && passnj && jet1_eta*jet2_eta < 0 && dijet_M>=800 && metnomuons>90 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0 && mindphi > 1.5){//pass sel
	Event lEvt;
	lEvt.run = run;
	lEvt.event = event;
	lEvt.lumi = lumi;
	isInserted = evtSet.insert(lEvt);
	if (!isInserted.second){
	  is_dupl = 1;
	  if (iT==0)
	    duplicateJ1J2 += 1;
	  else if (iT==1)
	    duplicateJ1J3 += 1;
	  else 
	    duplicateJ2J3 += 1;
	}
	else {
	  if (iT==0)
	    passJ1J2 += 1;
	  else if (iT==1){
	    passJ1J3 += 1;
	    evtSet.erase(lEvt);
	  }
	  else 
	    passJ2J3 += 1;
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
	    << label[0]<< " & " << duplicateJ1J2+passJ1J2 << " & " << duplicateJ1J2 << " & " << passJ1J2 << "\\\\" << std::endl
	    << label[1]<< " & " << duplicateJ1J3+passJ1J3 << " & " << duplicateJ1J3 << " & " << passJ1J3 << "\\\\"<< std::endl
	    << label[2]<< " & " << duplicateJ2J3+passJ2J3 << " & " << duplicateJ2J3 << " & " << passJ2J3 << "\\\\"<< std::endl
    ;

  const unsigned nVars = 19;
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
    "jet1_csv",
    "jet2_csv",
    "jet3_csv",
    "dijetmetnomu_scalarSum_pt",
    "dijetmetnomu_vectorialSum_pt",
    "dijetmetnomu_ptfraction",
    "n_jets_30",
    "n_jets_15",
    "n_jets_cjv_30"
  };

  std::string xaxis[nVars] = {
    ";M_{jj} (GeV)",
    ";METnoMu (GeV)",
    ";#Delta#eta_{jj}",
    ";#Delta#phi_{jj}",
    ";METnoMu/#sigma(METnoMu)",
    ";min #Delta#phi(ja-jb,METnoMu)",
    ";min #Delta#phi(all jets,METnoMu)",
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
    ";n_{jets} (p_{T}>15 GeV)",
    ";CJV jets (30 GeV)"
  };

  //unsigned nbins[nVars] = {75,100,44,30,35,30,50,50,100,50,50,10};
  float min[nVars] = {0,0,3.6,0,3,0,0,30,30,30,0,0,0,0,0,0,0,0,0};
  float max[nVars] = {3000,500,8,3.1416,10,3.1416,3.1416,300,300,300,1,1,1,1000,400,1,10,10,10};  

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
      std::string passpT = "(jet1_pt > 50 && jet2_pt > 40)";

      selection << "( ";
      if (iT>0) selection << "n_jets_30>2 && ";
      selection << "is_dupl==0 && ";
      selection << passtrig <<" && ";
      //if (iT==0) 
	selection << passpT << " && ";
	selection << "jet1_eta*jet2_eta<0 && dijet_M>=800 && metnomuons>90 && nvetomuons==0 && nselmuons==0 && nvetoelectrons==0 && nselelectrons==0 && jetmetnomu_mindphi>1.5";
      selection << " )";
      if (iT==2 && doReweighting) selection << " * (qcdW)";
      //selection << " * (1)";
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
      if (iV==7 && iT==0) {
	lname.str("");
	lname << vars[8] << ":" << vars[7] << ">>weight2D";// << label[iT];
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
    

    histsum[iV] = (TH1F*)hist[iV][2]->Clone("histsum");
    //histsum[iV]->Sumw2();
    //histsum[iV]->Add(hist[iV][2]);

    TH1F *histTmp = (TH1F*)histDataSubtr[iV][0]->Clone("histTmp");
    //histTmp->Sumw2();
    histTmp->Divide(histDataSubtr[iV][2]);
    grWeight[iV] = new TGraphErrors(histTmp);

    histsumSubtr[iV] = (TH1F*)histDataSubtr[iV][2]->Clone("histsumSubtr");
    //histsumSubtr[iV]->Sumw2();
    //histsumSubtr[iV]->Add(histDataSubtr[iV][2]);
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
    leg->AddEntry(histsumSubtr[iV],"Data-QCD","F");
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
  mycpT1->Print("qcd_reweighting_with_mindphi.pdf");

  TCanvas *mycpT12 = new TCanvas("pt1pt2","pt1pt2",1);
  mycpT12->cd();
  weight2D->Draw("colz");

  return 0;
}
