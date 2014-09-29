#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeAnalyser.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/NormPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/SimplePlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/MVATrain.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/AddFriends.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/Plotter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/HistPlotter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/SummaryTable.h"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "TColor.h"

namespace po=boost::program_options;
using namespace ic;

int main(int argc, char* argv[]){
  /*##########################################
  #                                          #
  #            SET UP OPTIONS                #
  #                                          #
  ##########################################*/

  //Input output and config options
  std::string cfg;
  std::string outputname;
  std::string inputfolder;
  std::string inputparams;
  std::string filelist;
  std::string basesel;
  std::string channel;
  std::string treename;
  double signalf = 1.0;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options
    ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))
    ("input_folder,i",           po::value<std::string>(&inputfolder)->default_value("../output_lighttree_withrle2/"))
    ("input_params,p",           po::value<std::string>(&inputparams)->default_value("../filelists/Dec18/ParamsDec18test.dat"))
    ("filelist,f",               po::value<std::string>(&filelist)->default_value("filelists/filelist.dat"))
    ("basesel",                  po::value<std::string>(&basesel)->default_value("jet1_eta*jet2_eta<0 && jet1_eta<4.7 && jet2_eta<4.7 && dijet_M>=600&&jet1_pt>50&&dijet_deta>3.6&& jet2_pt>60&&metnomuons>60&&metnomu_significance>3&&jetmetnomu_mindphi>1.5"))
    ("treename",                 po::value<std::string>(&treename)->default_value("LightTree"))
    ("signalFactor",             po::value<double>(&signalf)->default_value(1.))
    ("channel",                  po::value<std::string>(&channel)->default_value("nunu"));

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  /*##########################################
  #                                          #
  #          INSTANTIATE ANALYSER            #
  #                                          #
  ##########################################*/

  LTAnalyser* analysis = new LTAnalyser(outputname);
  analysis->SetTreeName(treename);

  analysis->AddFiles(filelist);
  std::cout<<"Taking input from: "<<inputfolder<<std::endl;
  analysis->SetInFolder(inputfolder);
  analysis->SetInputParams(inputparams);

  //Set selection step common to all categories
  //analysis->set_baseselection("passtrigger==1&&jet1_eta<4.7&&jet2_eta<4.7&& jet1_pt>"+boost::lexical_cast<std::string>(jet1ptcut)+"&& jet2_pt>"+boost::lexical_cast<std::string>(jet2ptcut)+" && dijet_M >"+boost::lexical_cast<std::string>(mjjcut)+"&& jet1_eta*jet2_eta<"+boost::lexical_cast<std::string>(etaprodcut)+"&& dijet_dphi<"+boost::lexical_cast<std::string>(dphicut)+"&& dijet_deta >"+boost::lexical_cast<std::string>(detacut));
  //analysis->set_baseselection("jet1_pt>50&&jet2_pt>50&&dijet_deta>3.6&&metnomu_significance>3&&jetmetnomu_mindphi>1.5");
  //analysis->set_baseselection("dijet_M>=1100&&l1met>=40&&jet1_pt>50&&jet2_pt>50&&metnomuons>60&&dijet_deta>3.6&&metnomu_significance>4.0&&jetmetnomu_mindphi>2.0");//BEST SO FAR
  //analysis->set_baseselection("jet1_eta*jet2_eta<0 && jet1_eta<4.7 && jet2_eta<4.7 && dijet_M>=1100&&jet1_pt>50&&dijet_deta>3.5&& jet2_pt>40&&metnomuons>90&&metnomu_significance>3.5&&jetmetnomu_mindphi>1.5");//&&jetunclet_mindphi<2.5&&jetunclet_mindphi>0.8");
  analysis->set_baseselection(basesel);
  
  /*##########################################
  #                                          #
  #            DEFINE MODULES                #
  #                                          #
  ##########################################*/

  //ADD FRIENDS TO GET EXTRA VARIABLES
  std::vector<std::string> setswithfriends;
  setswithfriends.push_back("WJets_taunu");
  setswithfriends.push_back("WJets_munu");
  setswithfriends.push_back("WJets_enu");
  setswithfriends.push_back("VBF-QCD");
  setswithfriends.push_back("MET");
  setswithfriends.push_back("PARKED");
  setswithfriends.push_back("PARKEDPLUSA");
  setswithfriends.push_back("VV");
  setswithfriends.push_back("Top");
  setswithfriends.push_back("ZJets_ll");
  setswithfriends.push_back("ZJets_ll_vbf");
  setswithfriends.push_back("ZJets_nunu");
  setswithfriends.push_back("sig125");
  
  AddFriends addfriends("addfriends");
  addfriends.set_frienddir("output/")
    .set_friendtreename("mvafriend")
    .set_sets(setswithfriends);

  std::vector<std::string> shape;
  std::vector<std::string> histTitle;
  //shape.push_back("BDT(12,-1.,0.2)");

  shape.push_back("jet1_pt(28,20.,300.)");  histTitle.push_back(";p_{T}^{j1} (GeV);entries");
  shape.push_back("jet2_pt(28,20.,300.)");  histTitle.push_back(";p_{T}^{j2} (GeV);entries");
  shape.push_back("jet3_pt(29,5.,295.)");  histTitle.push_back(";p_{T}^{j3} (GeV);entries");
  shape.push_back("metnomuons(25,80.,300.)"); histTitle.push_back(";METnoMu (GeV);entries");
  //shape.push_back("l1met(20,00.,200.)");histTitle.push_back(";L1MET (GeV);entries");
  shape.push_back("dijet_M(14,700.,2000.)");histTitle.push_back(";M_{jj} (GeV);entries");
  shape.push_back("jetmetnomu_mindphi(60,0,3.1416)");histTitle.push_back(";min #Delta#phi(j,METnoMu);entries");
  shape.push_back("alljetsmetnomu_mindphi(60,0,3.1416)");histTitle.push_back(";min #Delta#phi(j,METnoMu);entries");
  shape.push_back("metnomu_significance(101,2.9,13.)");histTitle.push_back(";METnoMu/#sigma(METnoMu);entries");
  //shape.push_back("dijet_sumeta(50,-10,10)");histTitle.push_back(";#eta_{j1}+#eta_{j2};entries");
  shape.push_back("ht(50,0,1000)");histTitle.push_back(";H_{T} (GeV);entries");
  //shape.push_back("jetunclet_mindphi(32,0,3.2)");histTitle.push_back(";min #Delta#phi(j,E_{T}^{uncl});entries");
  //shape.push_back("metnomuunclet_dphi(32,0,3.2)");histTitle.push_back(";#Delta#phi(METnoMu,E_{T}^{uncl};entries");
  shape.push_back("dijetmetnomu_scalarSum_pt(70,0,1400)");histTitle.push_back(";p_{T}^{jeta}+p_{T}^{jetb}+METnoMu;entries");
  shape.push_back("dijetmetnomu_vectorialSum_pt(20,0,400)");histTitle.push_back(";p_{T}(#vec{ja}+#vec{jb}+#vec{METnoMu});entries");
  shape.push_back("dijetmetnomu_ptfraction(20,0.,1.)");histTitle.push_back(";p_{T}^{dijet}/(p_{T}^{dijet}+METnoMu);entries");
  shape.push_back("n_jets_30(10,0,10)");histTitle.push_back(";njets (p_{T}>30 GeV);entries");
  shape.push_back("n_jets_15(10,0,10)");histTitle.push_back(";njets (p_{T}>15 GeV);entries");
  shape.push_back("n_jets_cjv_30(10,0,10)");histTitle.push_back(";CJV jets (30 GeV);entries");
  //shape.push_back("n_jets_cjv_20EB_30EE(5,0,5)");histTitle.push_back(";CJV jets (20 GeV EB, 30 GeV EE);entries");
  shape.push_back("dijet_dphi(30,0.,3.1416)");histTitle.push_back(";#Delta#phi_{jj};entries");
  shape.push_back("dijet_deta(18,3.4,7.)");histTitle.push_back(";#Delta#eta_{jj};entries");
  shape.push_back("jet1_csv(51,0,1.02)");histTitle.push_back(";CSV (jet 1);entries");
  shape.push_back("jet2_csv(51,0,1.02)");histTitle.push_back(";CSV (jet 2);entries");
  shape.push_back("jet3_csv(51,0,1.02)");histTitle.push_back(";CSV (jet 3);entries");
  shape.push_back("m_mumu(60,60,120)");histTitle.push_back(";M_{#mu#mu} (GeV);entries");
  shape.push_back("lep_mt(50,0.,200.)");histTitle.push_back(";m_{T}(lepton+MET (GeV);entries");
  shape.push_back("thrust(50,-10,5)");histTitle.push_back(";ln(#tau_{#perp}); entries");
  shape.push_back("thrust_minor(50,-10,5)");histTitle.push_back(";ln(T_{m}); entries");
  shape.push_back("thrust_met(50,-10,5)");histTitle.push_back(";ln(#tau_{#perp}); entries");
  shape.push_back("thrust_met_minor(50,-10,5)");histTitle.push_back(";ln(T_{m}); entries");

  assert(shape.size() == histTitle.size());

  /*
  if(channel!="taunu"){
    shape.push_back("jet2_pt(27,30.,300.)");histTitle.push_back(";p_{T}^{j1} (GeV);entries");
    shape.push_back("jet1_pt(27,30.,300.)");histTitle.push_back(";p_{T}^{j2} (GeV);entries");
    shape.push_back("metnomuons(25,50.,300.)");histTitle.push_back(";METnoMu (GeV);entries");
    shape.push_back("l1met(20,00.,200.)");histTitle.push_back(";L1MET (GeV);entries");
    shape.push_back("dijet_M(14,600.,2000.)");histTitle.push_back(";M_{jj} (GeV);entries");
    shape.push_back("jetmetnomu_mindphi(32,0.,3.2)");histTitle.push_back(";min #Delta#phi(j1/j2,METnoMu);entries");
    shape.push_back("alljetsmetnomu_mindphi(32,0.,3.2)");histTitle.push_back(";min #Delta#phi(all jets,METnoMu);entries");
    shape.push_back("metnomu_significance(25,3.,8.)");histTitle.push_back(";METnoMu/#sigma(METnoMu);entries");
    shape.push_back("dijet_sumeta(50,-10,10)");histTitle.push_back(";#eta_{j1}+#eta_{j2};entries");
    shape.push_back("ht(50,0,1000)");histTitle.push_back(";H_{T} (GeV);entries");
    shape.push_back("jetunclet_mindphi(32,0,3.2)");histTitle.push_back(";min #Delta#phi(j,E_{T}^{uncl});entries");
    shape.push_back("metnomuunclet_dphi(32,0,3.2)");histTitle.push_back(";#Delta#phi(METnoMu,E_{T}^{uncl};entries");
    shape.push_back("dijetmetnomu_scalarSum_pt(70,0,1400)");histTitle.push_back(";p_{T}^{jeta}+p_{T}^{jetb}+METnoMu;entries");
    shape.push_back("dijetmetnomu_vectorialSum_pt(20,0,400)");histTitle.push_back(";p_{T}(#vec{ja}+#vec{jb}+#vec{METnoMu});entries");
    shape.push_back("n_jets_cjv_30(5,0,5)");histTitle.push_back(";CJV jets (30 GeV);entries");
    shape.push_back("n_jets_cjv_20EB_30EE(5,0,5)");histTitle.push_back(";CJV jets (20 GeV EB, 30 GeV EE);entries");
    shape.push_back("dijet_dphi(30,0.,3.)");histTitle.push_back(";#Delta#phi_{jj};entries");
    shape.push_back("dijet_deta(18,3.4,7.)");histTitle.push_back(";#Delta#eta_{jj};entries");
    shape.push_back("lep_mt(20,0.,100.)");histTitle.push_back(";m_{T}(lepton+MET (GeV);entries");
    shape.push_back("dijetmetnomu_ptfraction(20,0.,1.)");histTitle.push_back(";p_{T}^{dijet}/(p_{T}^{dijet}+METnoMu);entries");
  }
  else{
shape.push_back("jet2_pt(14,30.,300.)");histTitle.push_back(";p_{T}^{j1} (GeV);entries");
    shape.push_back("jet1_pt(14,30.,300.)");histTitle.push_back(";p_{T}^{j2} (GeV);entries");
    shape.push_back("metnomuons(12,50.,300.)");histTitle.push_back(";METnoMu (GeV);entries");
    shape.push_back("l1met(10,00.,200.)");histTitle.push_back(";L1MET (GeV);entries");
    shape.push_back("dijet_M(7,600.,2000.)");histTitle.push_back(";M_{jj} (GeV);entries");
    shape.push_back("jetmetnomu_mindphi(16,0.,3.2)");histTitle.push_back(";min #Delta#phi(j1/j2,METnoMu);entries");
    shape.push_back("alljetsmetnomu_mindphi(16,0.,3.2)");histTitle.push_back(";min #Delta#phi(all jets,METnoMu);entries");
    shape.push_back("metnomu_significance(12,3.,8.)");histTitle.push_back(";METnoMu/#sigma(METnoMu);entries");
    shape.push_back("dijet_sumeta(25,-10,10)");histTitle.push_back(";#eta_{j1}+#eta_{j2};entries");
    shape.push_back("ht(25,0,1000)");histTitle.push_back(";H_{T} (GeV);entries");
    shape.push_back("jetunclet_mindphi(16,0,3.2)");histTitle.push_back(";min #Delta#phi(j,E_{T}^{uncl});entries");
    shape.push_back("metnomuunclet_dphi(16,0,3.2)");histTitle.push_back(";#Delta#phi(METnoMu,E_{T}^{uncl};entries");
    shape.push_back("dijetmetnomu_scalarSum_pt(35,0,1400)");histTitle.push_back(";p_{T}^{jeta}+p_{T}^{jetb}+METnoMu;entries");
    shape.push_back("dijetmetnomu_vectorialSum_pt(10,0,400)");histTitle.push_back(";p_{T}(#vec{ja}+#vec{jb}+#vec{METnoMu});entries");
    shape.push_back("n_jets_cjv_30(5,0,5)");histTitle.push_back(";CJV jets (30 GeV);entries");
    shape.push_back("n_jets_cjv_20EB_30EE(5,0,5)");histTitle.push_back(";CJV jets (20 GeV EB, 30 GeV EE);entries");
    shape.push_back("dijet_dphi(15,0.,3.)");histTitle.push_back(";#Delta#phi_{jj};entries");
    shape.push_back("dijet_deta(9,3.4,7.)");histTitle.push_back(";#Delta#eta_{jj};entries");
    shape.push_back("lep_mt(10,0.,100.)");histTitle.push_back(";m_{T}(lepton+MET (GeV);entries");
    shape.push_back("dijetmetnomu_ptfraction(10,0.,1.)");histTitle.push_back(";p_{T}^{dijet}/(p_{T}^{dijet}+METnoMu);entries");
  }
  */

  std::string dataset="PARKEDPLUSA";
  std::string dataextrasel="&&((((run>=190456)&&(run<=193621))&&passtrigger==1)||(((run>=193833)&&(run<=196531))&&passparkedtrigger1==1)||(((run>=203777)&&(run<=208686))&&passparkedtrigger2==1))&&l1met>40";
  std::string sigcat;
  std::string zextrasigcat;

  std::string nunucat="nvetomuons==0&&nvetoelectrons==0 && thrust > -5";
  std::string nunuzcat=" && thrust > -5";
  
  std::string mumucat="nselmuons==2&&nvetomuons==2&&nvetoelectrons==0&&m_mumu>60&&m_mumu<120";
  std::string mumuzcat="&&nselmuons==2&&nvetomuons==2&&m_mumu>60&&m_mumu<120";//zmumu

  std::string munucat="nselmuons==1&&nvetomuons==1&&nvetoelectrons==0";
  std::string munuzcat="&&nselmuons==1&&nvetomuons==1&&nvetoelectrons==0&&m_mumu>60&&m_mumu<120";//wmu

  std::string enucat="nselelectrons==1&&nvetomuons==0&&nvetoelectrons==1";
  std::string enuzcat="&&nselelectrons==1&&nvetoelectrons==1";//wel

  std::string taunucat="ntaus==1&&nvetomuons==0&&nvetoelectrons==0&&lep_mt>20";
  std::string taunuzcat="&&ntaus==1&&nvetoelectrons==0";//wtau

  //std::string qcdcat="nvetoelectrons==0&&nvetomuons==0&&dijetmetnomu_ptfraction>0.6";
  //std::string qcdzcat="&&dijetmetnomu_ptfraction>0.6";//QCD
  std::string qcdcat="nvetoelectrons==0&&nvetomuons==0&&thrust<-5";
  std::string qcdzcat="&&thrust<-5";//QCD

  if(channel=="nunu"){//nunu
    sigcat=nunucat;
    zextrasigcat=nunuzcat;
  }
  else if(channel=="mumu"){//zmumu
    sigcat=mumucat;
    zextrasigcat=mumuzcat;
  }
  else if(channel=="munu"){//wmu
    sigcat=munucat;
    zextrasigcat=munuzcat;
  }
  else if(channel=="enu"){//wel
    sigcat=enucat;
    zextrasigcat=enuzcat;
  }
  else if(channel=="taunu"){//wtau
    sigcat=taunucat;
    zextrasigcat=taunuzcat;
  }
  else if(channel=="qcd"){//QCD
    sigcat=qcdcat;
    zextrasigcat=qcdzcat;
  }
  else{
    std::cout<<"Error: Channel "<<channel<<" not recognised, exiting"<<std::endl;
    return 1;
  }




    
  //DATA SHAPE GENERATION
  DataShape data("data");
  data.set_dataset(dataset)
    .set_dirname("data_obs")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+dataextrasel);

  DataShape signal("signal");
  signal.set_dataset("sig125")
    .set_dirname("qqH")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape vv("vv");
  vv.set_dataset("VV")
    .set_dirname("vv")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape top("top");
  top.set_dataset("Top")
    .set_dirname("top")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape wgamma("wgamma");
  wgamma.set_dataset("WGamma")
    .set_dirname("wg")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape znunuraw("znunuraw");
  znunuraw.set_dataset("ZJets_nunu")
    .set_dirname("zvv")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape zmumuraw("zmumuraw");
  zmumuraw.set_dataset("ZJets_ll_all")
    .set_dirname("zmumu")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);
  
  DataShape wmunuraw("wmunuraw");
  wmunuraw.set_dataset("WJets_munu")
    .set_dirname("wmuraw")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape wenuraw("wenuraw");
  wenuraw.set_dataset("WJets_enu")
    .set_dirname("welraw")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape wtaunuraw("wtaunuraw");
  wtaunuraw.set_dataset("WJets_taunu")
    .set_dirname("wtauraw")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  DataShape QCDraw("QCDraw");
  QCDraw.set_dataset("VBF-QCD")
    .set_dirname("qcdraw")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat);

  //ZBKG SHAPE GENERATION
  std::vector<std::string> Zcontbkgsets;
  Zcontbkgsets.push_back("VV");
  Zcontbkgsets.push_back("Top");
  Zcontbkgsets.push_back("WGamma");
  Zcontbkgsets.push_back("WJets_enu");
  Zcontbkgsets.push_back("WJets_munu");
  Zcontbkgsets.push_back("WJets_taunu");
  
  DataZNormShape zmumu("zmumu");
  zmumu.set_sigmcewkset("ZJets_ll_vbf")
    .set_shape(shape)
    .set_dirname("zvv")
    .set_sigmcqcdset("ZJets_ll")
    .set_contmcewkset("ZJets_ll_vbf")
    .set_contmcqcdset("ZJets_ll")
    .set_contbkgset(Zcontbkgsets)
    .set_contdataset(dataset)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat("m_mumu_gen>80&&m_mumu_gen<100"+zextrasigcat)
    .set_contcat(mumucat)//"nvetoelectrons==0 && nvetomuons==2 && nselmuons==2&&m_mumu>60&&m_mumu<120")
    .set_contdataweight("weight_nolep")
    .set_sigmainccontewk(303)
    .set_sigmainccontqcd(3503700./3)
    .set_sigmaincsigewk(460*3)
    .set_sigmaincsigqcd(3503700./3*5.651)
    .set_ngenincewk(5781.91)
    .set_ngenincqcd(22789300)
    .set_ngenmassfilteredewk(4226.53)
    .set_ngenmassfilteredqcd(20334900);
    
  DataNormShape zmumuinzcont("zmumuinzcont");
  zmumuinzcont.set_sigmcset("ZJets_ll_all")
    .set_shape(shape)
    .set_dirname("zvv")
    .set_contmcset("ZJets_ll_all")
    .set_contbkgset(Zcontbkgsets)
    .set_contdataset(dataset)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigmcweight("total_weight_leptight")
    .set_sigcat("m_mumu_gen>80&&m_mumu_gen<100"+zextrasigcat)
    .set_contcat(mumucat)
    .set_contdataweight("weight_nolep");
;//"nvetoelectrons==0 && nvetomuons==2 && nselmuons==2&&m_mumu>60&&m_mumu<120");

  //WBKG SHAPE GENERATION
  std::vector<std::string> Wcontbkgsets; //List of sets for ncbkg
  Wcontbkgsets.push_back("VV");
  Wcontbkgsets.push_back("Top");
  Wcontbkgsets.push_back("WGamma");
//   Wcontbkgsets.push_back("ZJets_ll");
//   Wcontbkgsets.push_back("ZJets_ll_vbf");
//   Wcontbkgsets.push_back("ZJets_nunu");

  DataNormShape wmunu("wmunu");
  wmunu.set_sigmcset("WJets_munu")
    .set_shape(shape)
    .set_dirname("wmu")
    .set_contmcset("WJets_munu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat(sigcat)
    .set_contcat(munucat)
    .set_contdataweight("weight_nolep");
  if((channel=="enu")||(channel=="munu")){
    wmunu.set_sigmcweight("total_weight_leptight");
  }
  

  DataNormShape wenu("wenu");
  wenu.set_sigmcset("WJets_enu")
    .set_shape(shape)
    .set_dirname("wel")
    .set_contmcset("WJets_enu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat(sigcat)
    .set_contcat(enucat)
    .set_contdataweight("weight_nolep");
  if((channel=="enu")||(channel=="munu")){
    wenu.set_sigmcweight("total_weight_leptight");
  }

  DataNormShape wtaunu("wtaunu");
  wtaunu.set_sigmcset("WJets_taunu")
    .set_shape(shape)
    .set_dirname("wtau")
    .set_contmcset("WJets_taunu")
    .set_contdataset(dataset)
    .set_contbkgset(Wcontbkgsets)
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat(sigcat)
    .set_contcat(taunucat)//"ntaus>=1&&nvetoelectrons ==0 && nvetomuons==0&&lep_mt>20")
    .set_sigmcweight("total_weight_lepveto")
    .set_contmcweight("total_weight_lepveto")
    .set_contdataweight("weight_nolep");
  if((channel=="enu")||(channel=="munu")){
    wtaunu.set_sigmcweight("total_weight_leptight");
  }

  //QCDBKG
  std::vector<std::string> QCDcontbkgsets; //list of sets for ncbkg
  QCDcontbkgsets.push_back("VV");
  QCDcontbkgsets.push_back("Top");
  QCDcontbkgsets.push_back("WGamma");
  QCDcontbkgsets.push_back("ZJets_ll");
  QCDcontbkgsets.push_back("ZJets_ll_vbf");
  QCDcontbkgsets.push_back("WJets_enu");
  QCDcontbkgsets.push_back("WJets_munu");
  QCDcontbkgsets.push_back("WJets_taunu");

  std::vector<std::string> QCDcontbkgextrafactordir;//list of dirs with data driven weights for above backgrounds
  QCDcontbkgextrafactordir.push_back("");
  QCDcontbkgextrafactordir.push_back("");
  QCDcontbkgextrafactordir.push_back("");
  QCDcontbkgextrafactordir.push_back("zvv");
  QCDcontbkgextrafactordir.push_back("zvv");
  QCDcontbkgextrafactordir.push_back("wel");
  QCDcontbkgextrafactordir.push_back("wmu");
  QCDcontbkgextrafactordir.push_back("wtau");
  
  std::vector<int> QCDcontbkgisz;
  QCDcontbkgisz.push_back(0);
  QCDcontbkgisz.push_back(0);
  QCDcontbkgisz.push_back(0);
  if(channel!="mumu"){
    QCDcontbkgisz.push_back(2);
    QCDcontbkgisz.push_back(1);
  }
  else{
    QCDcontbkgisz.push_back(0);
    QCDcontbkgisz.push_back(0);
  }
  QCDcontbkgisz.push_back(0);
  QCDcontbkgisz.push_back(0);
  QCDcontbkgisz.push_back(0);

  /*
  DataNormShape QCD("QCD");
  QCD.set_sigmcset("VBF-QCD")//VBF-QCD")
    .set_shape(shape)
    .set_dirname("qcd")
    .set_contmcset("VBF-QCD")//VBF-QCD")
    .set_contdataset(dataset)
    .set_contbkgset(QCDcontbkgsets)
    .set_contbkgextrafactordir(QCDcontbkgextrafactordir)
    .set_contbkgisz(QCDcontbkgisz)
    .set_sigmcweight("total_weight_lepveto")
    .set_contmcweight("total_weight_lepveto")
    .set_contdataweight("weight_nolep")
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat(sigcat)
    .set_zcontcat("m_mumu_gen>80&&m_mumu_gen<100")
    .set_contcat("nvetoelectrons==0&&nvetomuons==0");
  if(channel=="enu"||channel=="munu"){
    QCD.set_sigmcweight("total_weight_leptight");
  }
  */

  //QCD from DATA SHAPE GENERATION
  DataNormShape QCD("QCD");
  QCD.set_sigmcset(dataset)//VBF-QCD")
    .set_shape(shape)
    .set_dirname("qcd")
    .set_contmcset(dataset)//VBF-QCD")
    .set_contdataset(dataset)
    .set_contbkgset(QCDcontbkgsets)
    .set_contbkgextrafactordir(QCDcontbkgextrafactordir)
    .set_contbkgisz(QCDcontbkgisz)
    .set_sigmcweight("weight_nolep")
    .set_contmcweight("weight_nolep")
    .set_contdataweight("weight_nolep")
    .set_basesel(analysis->baseselection())
    .set_contdataextrasel(dataextrasel)
    .set_sigcat(qcdcat+dataextrasel)
    .set_zcontcat("m_mumu_gen>80&&m_mumu_gen<100")
    .set_contcat(nunucat);


  //NORMALISED PLOTS FOR REFEREE
  std::vector<std::string> ewksets; //List of sets for ewk
  ewksets.push_back("VV");
  ewksets.push_back("Top");
  ewksets.push_back("ZJets_ll");
  ewksets.push_back("ZJets_ll_vbf");
  ewksets.push_back("ZJets_nunu");
  ewksets.push_back("WJets_enu");
  ewksets.push_back("WJets_munu");
  ewksets.push_back("WJets_taunu");

  std::vector<std::string> shapes; //List of shapes to draw
  shapes.push_back("dijet_M(370,150.,2000.)");
  shapes.push_back("dijet_deta(160,0.,8.)");
  shapes.push_back("dijet_dphi(310,0.,3.1)");
  shapes.push_back("met(80,0.,400.)");
  shapes.push_back("cjvjetpt(100,0.,100.)");
  shapes.push_back("met(50,0.,500.)");
  shapes.push_back("met_x(50,0.,500.)");
  shapes.push_back("met_y(50,0.,500.)");

  NormPlots normplots("normplots");
  normplots.set_qcdset("QCD")
    .set_sigset("sig125")
    .set_ewkset(ewksets)
    .set_cat("")
    .set_basesel("jet1_eta<4.7&&jet2_eta<4.7&&jet1_pt>50&&jet2_pt>50&&nvetoelectrons==0 && nvetomuons==0&&dijet_M>150&&met>130")
    .set_shapes(shapes);


  //MVA TRAIN
  std::vector<std::string> sigsets;
  sigsets.push_back("sig125");
  std::vector<std::string> bkgsets;
  bkgsets.push_back("QCD");
  std::vector<std::string> variables;
  variables.push_back("jet2_eta");
  variables.push_back("met");
  variables.push_back("met_significance");
  variables.push_back("mht");
  variables.push_back("jet2met_dphi");
  variables.push_back("jetmet_mindphi");
  variables.push_back("jetunclet_mindphi");
  variables.push_back("metunclet_dphi");
  variables.push_back("dijetmet_scalarSum_pt");
  variables.push_back("dijetmet_vectorialSum_pt");
  variables.push_back("dijetmet_ptfraction");
  variables.push_back("jet1met_scalarprod");
  variables.push_back("jet2met_scalarprod");
  variables.push_back("jet1met_scalarprod_frac := jet1met_scalarprod/met");
  variables.push_back("jet2met_scalarprod_frac := jet2met_scalarprod/met");
  std::vector<std::string> specvariables;

  MVATrain mvatrainer("mvatrainer");
  mvatrainer.set_sigsets(sigsets)
    .set_bkgsets(bkgsets)
    .set_variables(variables)
    .set_specvariables(specvariables)
    .set_basesel("passtrigger==1&&nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7")
    .set_sigcat("")
    .set_bkgcat("");

  //PLOTTER
//   std::vector<std::vector<std::string> > plottersets;
//   std::vector<std::string> sets1;
//   sets1.push_back("sig125");
//   std::vector<std::string> sets2;
//   sets2.push_back("WJets_enu");
//   plottersets.push_back(sets1);
//   plottersets.push_back(sets2);
 
//   std::vector<std::string> plottervariables;
//   plottervariables.push_back("BDT(40,-1.,1.)");

//   Plotter plotter("plotter");
//   plotter.set_sets(plottersets)
//     .set_shapes(plottervariables)
//     .set_basesel("passtrigger==1&&nvetomuons==0&&nvetoelectrons==0&&jet1_pt>50&&jet2_pt>50&&jet1_eta<4.7&&jet2_eta<4.7&&met_significance>3&&dijet_deta>3.6&&jetmet_mindphi>1.5")
//     .set_cat("");
  
  //HISTPLOTTER
  std::vector<std::string> shapevec;
  for(unsigned ishape=0;ishape<shape.size();ishape++){
    std::vector<std::string> strs;
    boost::split(strs, shape[ishape], boost::is_any_of("("));
    shapevec.push_back(strs[0]);
  }

  std::vector<LTPlotElement> elementvec;
  LTPlotElement dataele;
  dataele.set_is_data(true)
    .set_scale(1)
    .set_legname("Data")
    .set_is_inrationum(true)
    .set_sample("data_obs");

  LTPlotElement wmunuele;
  wmunuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange-4)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W#rightarrow#mu#nu")
    .set_sample("wmu");

  LTPlotElement wenuele;
  wenuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange  + 2)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W#rightarrow e#nu")
    .set_sample("wel");

  LTPlotElement wtaunuele;
  wtaunuele.set_is_data(false)
    .set_scale(1)
    .set_color(kOrange + 4)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W#rightarrow#tau#nu")
    .set_sample("wtau");

  LTPlotElement znunuele;
  znunuele.set_is_data(false)
    .set_scale(1)
    .set_color(kAzure  + 2)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z#rightarrow#nu#nu")
    .set_sample("zvv");

  LTPlotElement zmumuele;
  zmumuele.set_is_data(false)
    .set_scale(1)
    .set_color(kAzure  + 2)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z#rightarrow#mu#mu")
    .set_sample("zmumu");

  LTPlotElement qcdele;
  qcdele.set_is_data(false)
    .set_scale((64903-16011.1)/13016.8)
    .set_color(kMagenta-10)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("QCD")
    .set_sample("qcd");

  LTPlotElement vvele;
  vvele.set_is_data(false)
    .set_scale(1)
    .set_color(kGreen-5)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("VV")
    .set_sample("vv");

  LTPlotElement wgele;
  wgele.set_is_data(false)
    .set_scale(1)
    .set_color(kGreen-3)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("WGamma")
    .set_sample("wg");

  LTPlotElement topele;
  topele.set_is_data(false)
    .set_scale(1)
    .set_color(kBlue-8)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Top")
    .set_sample("top");

  std::ostringstream legname;
  legname << "Signal #times " << signalf; 

  LTPlotElement sigele;
  sigele.set_is_data(false)
    .set_scale(signalf)
    .set_color(kRed)
    .set_in_stack(false)
    .set_legname(legname.str())
    .set_sample("qqH");

  elementvec.push_back(dataele);
  elementvec.push_back(wmunuele);
  elementvec.push_back(wenuele);
  elementvec.push_back(wtaunuele);
  //  elementvec.push_back(zmumuele);
  elementvec.push_back(znunuele);
  elementvec.push_back(qcdele);
  elementvec.push_back(vvele);
  elementvec.push_back(wgele);
  elementvec.push_back(topele);
  elementvec.push_back(sigele);

  HistPlotter plotter("plotter");
  plotter.set_dirname("ControlPlots")
    .set_do_ratio(true)
    .set_elements(elementvec)
    .set_histTitles(histTitle)
    .set_shapes(shapevec);
  
  std::vector<std::string> dirvec;
  dirvec.push_back("wel");
  dirvec.push_back("wmu");
  dirvec.push_back("wtau");
  dirvec.push_back("zvv");
  dirvec.push_back("qcd");
  dirvec.push_back("vv");
  dirvec.push_back("wg");  
  dirvec.push_back("top");
  dirvec.push_back("qqH");
  dirvec.push_back("data_obs");

  SummaryTable summary("summary");
  summary.set_shape(shapevec)
    .set_dirs(dirvec);
  

  /*##########################################
  #                                          #
  #   SET UP ANALYSIS SEQUENCE AND RUN       #
  #                                          #
  ##########################################*/
  
  //analysis->AddModule(&addfriends);
  //analysis->AddModule(&mvatrainer);
  //analysis->AddModule(&normplots);
  analysis->AddModule(&wmunu);
  analysis->AddModule(&wenu);
  analysis->AddModule(&wtaunu);
  if(channel!="mumu"){
    analysis->AddModule(&zmumu);
  }
  else analysis->AddModule(&zmumuinzcont);
  analysis->AddModule(&QCD);
  //analysis->AddModule(&wmunuraw);
  //analysis->AddModule(&wenuraw);
  //analysis->AddModule(&wtaunuraw);  
  //analysis->AddModule(&QCDraw);
   //analysis->AddModule(&zmumuraw);
   //analysis->AddModule(&znunuraw);
  analysis->AddModule(&vv);
  analysis->AddModule(&wgamma);
  analysis->AddModule(&top);
  analysis->AddModule(&data);
  analysis->AddModule(&signal);
  analysis->AddModule(&plotter);
  analysis->AddModule(&summary);

  analysis->RunAnalysis();

  return 0;

}
