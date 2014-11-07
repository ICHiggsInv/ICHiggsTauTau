#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/QCDNormShape.h"
#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TVectorD.h"
#include <map>
#include "boost/lexical_cast.hpp"

namespace ic{

  QCDNormShape::QCDNormShape(std::string name) : LTModule(name){
    sigqcdweight_="total_weight_lepveto";
    contqcdweight_="total_weight_lepveto";
    contdataweight_="weight_nolep";
    contdataextrasel_="";
    contbkgextrasel_="";
    contqcdextrasel_="";
    contqcdbkgextrasel_="";
    sigqcdextrasel_="";
    sigbkgextrasel_="";
    sigcontextrafactor_=1.;
    std::vector<std::string> shapes;
    shapes.push_back("jet2_pt(200,0.,1000.)");
    shape_=shapes;
    dirname_="";
  };

  QCDNormShape::~QCDNormShape(){ ;};

  int QCDNormShape::Init(TFile* fs){
    fs_=fs;
    std::cout<<"Initialisation info for "<<module_name_<<":"<<std::endl;
    std::cout<<"  Signal QCD set is: "<<sigqcdset_<<std::endl;
    std::cout<<"  Control QCD set is: "<<contqcdset_<<std::endl;
    std::cout<<"  Control data set is: "<<contdataset_<<std::endl;
    std::cout<<"  Base selection is: "<<basesel_<<std::endl;
    std::cout<<"  Signal extra selection is: "<<sigcat_<<std::endl;
    std::cout<<"  Control extra selection is: "<<contcat_<<std::endl;
    std::cout<<"  Signal QCD weight is: "<<sigqcdweight_<<std::endl;
    std::cout<<"  Control QCD weight is: "<<contqcdweight_<<std::endl;
    std::cout<<"  Control data weight is: "<<contdataweight_<<std::endl;
    if((shapename_.size()!=shape_.size())&&shapename_.size()!=0){
      std::cout<<"WARNING: different numbers of shape names and shapes expect errors!"<<std::endl;
    }
    return 0;
  };

  int QCDNormShape::Run(LTFiles* filemanager){
    std::cout<<module_name_<<":"<<std::endl;

    TFile *file=fs_;
    TDirectory* dir;
    if(dirname_==""){
      dir=file->mkdir(sigqcdset_.c_str());
    }
    else if(!fs_->GetDirectory(dirname_.c_str())){
      dir=file->mkdir(dirname_.c_str());
    }
    else{
      dir=fs_->GetDirectory(dirname_.c_str());
    }
    dir->cd();
    //Get Shapes for NSQCD, NCQCD, NCQCDBkg, NCData and NCBkg
    std::cout<<"  Getting control QCD shape"<<std::endl;
    TH1F  contqcdshape = filemanager->GetSetShape(contqcdset_,"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contqcdextrasel_),contqcdweight_,false);

    std::cout<<"  Getting control QCD Backgrounds shape"<<std::endl;
    TH1F contqcdbkgshape;
    bool firstbkg=true;
    //IF NO DIRS TO GET DATA DRIVEN WEIGHTS, GET SETS SHAPE
    if(contqcdbkgextrafactordir_.size()==0){
      contqcdbkgshape= filemanager->GetSetsShape(contqcdbkgset_,"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contqcdbkgextrasel_),contqcdweight_,false);
    }
    //WEIGHT INDIVIDUAL SHAPES BY DATA DRIVEN WEIGHTS
    else{
      if(contqcdbkgextrafactordir_.size()==contqcdbkgset_.size()){
	for(unsigned iBkg=0;iBkg<contqcdbkgset_.size();iBkg++){
	  double extrafactor=1;
	  TDirectory* extrafactordir;
	  TVectorD *ddweight;
	  if(contqcdbkgextrafactordir_[iBkg]==""){
	    extrafactor=1;
	  }
	  else if(!fs_->GetDirectory(contqcdbkgextrafactordir_[iBkg].c_str())){
	    std::cout<<"Directory "<<contqcdbkgextrafactordir_[iBkg]<<" does not exist, exiting"<<std::endl;
	    return 1;
	  }
	  else{
	    extrafactordir=fs_->GetDirectory(contqcdbkgextrafactordir_[iBkg].c_str());
	    ddweight = (TVectorD*)extrafactordir->Get("ddweight");
	    if(contqcdbkgisz_.size()==contqcdbkgset_.size()){
	      if(contqcdbkgisz_[iBkg]==0){
		extrafactor=(*ddweight)[0];
	      }
	      if(contqcdbkgisz_[iBkg]==1){
		extrafactor=(*ddweight)[0];//EWK WEIGHT
	      }
	      if(contqcdbkgisz_[iBkg]==2){
		extrafactor=(*ddweight)[1];//QCD WEIGHT
	      }
	    }
	  }
	  TH1F nextbkg;
	  //SPECIAL CASE FOR Z BKG
	  if(contqcdbkgisz_.size()==contqcdbkgset_.size()){
	    if(contqcdbkgisz_[iBkg]!=0){
	      //GET Z SHAPE AND WEIGHT IT
	      nextbkg=filemanager->GetSetShape(contqcdbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,zcontqcdcat_,"weight_nolep*"+boost::lexical_cast<std::string>(extrafactor),false);
	    }
	    else nextbkg=filemanager->GetSetShape(contqcdbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contqcdbkgextrasel_),contqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  }
	  else nextbkg=filemanager->GetSetShape(contqcdbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contqcdbkgextrasel_),contqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  
	  //ADD HISTOGRAMS TOGETHER
	  if(firstbkg){
	    contqcdbkgshape=nextbkg;
	    firstbkg=false;
	  }
	  else{
	    contqcdbkgshape.Add(&nextbkg);
	  }
	}
      }
      else{
	std::cout<<"Extra factor dir vector must be same size as bkg set vector expect errors!"<<std::endl;
	return 1;
      }
    }

    std::cout<<"  Getting control Data Backgrounds shape"<<std::endl;
    TH1F contbkgshape;
    firstbkg=true;
    //IF NO DIRS TO GET DATA DRIVEN WEIGHTS, GET SETS SHAPE
    if(contbkgextrafactordir_.size()==0){
      contbkgshape= filemanager->GetSetsShape(contbkgset_,"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contbkgextrasel_),contqcdweight_,false);
    }
    //WEIGHT INDIVIDUAL SHAPES BY DATA DRIVEN WEIGHTS
    else{
      if(contbkgextrafactordir_.size()==contbkgset_.size()){
	for(unsigned iBkg=0;iBkg<contbkgset_.size();iBkg++){
	  double extrafactor=1;
	  TDirectory* extrafactordir;
	  TVectorD *ddweight;
	  if(contbkgextrafactordir_[iBkg]==""){
	    extrafactor=1;
	  }
	  else if(!fs_->GetDirectory(contbkgextrafactordir_[iBkg].c_str())){
	    std::cout<<"Directory "<<contbkgextrafactordir_[iBkg]<<" does not exist, exiting"<<std::endl;
	    return 1;
	  }
	  else{
	    extrafactordir=fs_->GetDirectory(contbkgextrafactordir_[iBkg].c_str());
	    ddweight = (TVectorD*)extrafactordir->Get("ddweight");
	    if(contbkgisz_.size()==contbkgset_.size()){
	      if(contbkgisz_[iBkg]==0){
		extrafactor=(*ddweight)[0];
	      }
	      if(contbkgisz_[iBkg]==1){
		extrafactor=(*ddweight)[0];//EWK WEIGHT
	      }
	      if(contbkgisz_[iBkg]==2){
		extrafactor=(*ddweight)[1];//QCD WEIGHT
	      }
	    }
	  }
	  TH1F nextbkg;
	  //SPECIAL CASE FOR Z BKG
	  if(contbkgisz_.size()==contbkgset_.size()){
	    if(contbkgisz_[iBkg]!=0){
	      //GET Z SHAPE AND WEIGHT IT
	      nextbkg=filemanager->GetSetShape(contbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,zcontcat_,"weight_nolep*"+boost::lexical_cast<std::string>(extrafactor),false);
	    }
	    else nextbkg=filemanager->GetSetShape(contbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contbkgextrasel_),contqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  }
	  else nextbkg=filemanager->GetSetShape(contbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(contcat_+contbkgextrasel_),contqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  
	  //ADD HISTOGRAMS TOGETHER
	  if(firstbkg){
	    contbkgshape=nextbkg;
	    firstbkg=false;
	  }
	  else{
	    contbkgshape.Add(&nextbkg);
	  }
	}
      }
      else{
	std::cout<<"Extra factor dir vector must be same size as bkg set vector expect errors!"<<std::endl;
	return 1;
      }
    }

    std::cout<<"  Getting control Data shape"<<std::endl;
    TH1F  contdatashape = filemanager->GetSetShape(contdataset_,"jet2_pt(200,0.,1000.)",basesel_,contcat_+contdataextrasel_,contdataweight_,false);

    std::cout<<"  Getting signal QCD Backgrounds shape"<<std::endl;
    TH1F sigbkgshape;
    firstbkg=true;
    //IF NO DIRS TO GET DATA DRIVEN WEIGHTS, GET SETS SHAPE
    if(sigbkgextrafactordir_.size()==0){
      sigbkgshape= filemanager->GetSetsShape(sigbkgset_,"jet2_pt(200,0.,1000.)",basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_,false);
    }
    //WEIGHT INDIVIDUAL SHAPES BY DATA DRIVEN WEIGHTS
    else{
      if(sigbkgextrafactordir_.size()==sigbkgset_.size()){
	for(unsigned iBkg=0;iBkg<sigbkgset_.size();iBkg++){
	  double extrafactor=1;
	  TDirectory* extrafactordir;
	  TVectorD *ddweight;
	  if(sigbkgextrafactordir_[iBkg]==""){
	    extrafactor=1;
	  }
	  else if(!fs_->GetDirectory(sigbkgextrafactordir_[iBkg].c_str())){
	    std::cout<<"Directory "<<sigbkgextrafactordir_[iBkg]<<" does not exist, exiting"<<std::endl;
	    return 1;
	  }
	  else{
	    extrafactordir=fs_->GetDirectory(sigbkgextrafactordir_[iBkg].c_str());
	    ddweight = (TVectorD*)extrafactordir->Get("ddweight");
	    if(sigbkgisz_.size()==sigbkgset_.size()){
	      if(sigbkgisz_[iBkg]==0){
		extrafactor=(*ddweight)[0];
	      }
	      if(sigbkgisz_[iBkg]==1){
		extrafactor=(*ddweight)[0];//EWK WEIGHT
	      }
	      if(sigbkgisz_[iBkg]==2){
		extrafactor=(*ddweight)[1];//QCD WEIGHT
	      }
	    }
	  }
	  TH1F nextbkg;
	  //SPECIAL CASE FOR Z BKG
	  if(sigbkgisz_.size()==sigbkgset_.size()){
	    if(sigbkgisz_[iBkg]!=0){
	      //GET Z SHAPE AND WEIGHT IT
	      nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,zsigcat_,"weight_nolep*"+boost::lexical_cast<std::string>(extrafactor),false);
	    }
	    else nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  }
	  else nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],"jet2_pt(200,0.,1000.)",basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	  
	  //ADD HISTOGRAMS TOGETHER
	  if(firstbkg){
	    sigbkgshape=nextbkg;
	    firstbkg=false;
	  }
	  else{
	    sigbkgshape.Add(&nextbkg);
	  }
	}
      }
      else{
	std::cout<<"Extra factor dir vector must be same size as bkg set vector expect errors!"<<std::endl;
	return 1;
      }
    }
    //Integrate over shape to get number in each region
    double ncqcd = Integral(&contqcdshape);
    double ncqcdbkg = Integral(&contqcdbkgshape);
    double ncdata = Integral(&contdatashape);
    double ncbkg = Integral(&contbkgshape);
    double ncqcderr = Error(&contqcdshape);
    double ncqcdbkgerr = Error(&contqcdbkgshape);
    double ncdataerr = Error(&contdatashape);
    double ncbkgerr = Error(&contbkgshape);
    double nsbkg = Integral(&sigbkgshape);
    double nsbkgerr = Error(&sigbkgshape);

    std::cout<<"  ncqcd: "<<ncqcd<<"+-"<<ncqcderr<<", ncqcdbkg: " <<ncqcdbkg<<"+/-"<<ncqcdbkgerr<<", ncdata: "<<ncdata<<"+-"<<ncdataerr<<", ncbkg: "<<ncbkg<<"+-"<<ncbkgerr<<std::endl;
    std::cout<<"  nsbkg: "<<nsbkg<<"+-"<<nsbkgerr<<std::endl;

    //Calculate weight
    double weight=(ncdata-ncbkg)/(ncqcd-ncqcdbkg);
    double weighterrdatastatfrac=ncdataerr/(ncdata-ncbkg);
    double weighterrmcstatfrac=sqrt(((ncbkgerr/(ncdata-ncbkg))*(ncbkgerr/(ncdata-ncbkg)))+((ncqcderr/(ncqcd-ncqcdbkg))*(ncqcderr/(ncqcd-ncqcdbkg)))+((ncqcdbkgerr/(ncqcd-ncqcdbkg))*(ncqcdbkgerr/(ncqcd-ncqcdbkg))));
    std::cout<<"  weight: "<<weight<<"+-"<<weight*weighterrdatastatfrac<<"(data stat.)+-"<<weight*weighterrmcstatfrac<<"(MC stat.)"<<std::endl;
    TVectorD weightvec(1);
    TVectorD errvec(2);
    weightvec[0]=weight*sigcontextrafactor_;

    std::cout<<"  Getting signal QCD shapes"<<std::endl;

    for(unsigned iShape=0;iShape<shape_.size();iShape++){
      TH1F  sigqcdshape = filemanager->GetSetShape(sigqcdset_,shape_[iShape],basesel_,(sigcat_+sigqcdextrasel_),sigqcdweight_,false);


      //std::cout<<"  Getting signal QCD Backgrounds shape"<<std::endl;
      TH1F sigbkgshape;
      firstbkg=true;
      //IF NO DIRS TO GET DATA DRIVEN WEIGHTS, GET SETS SHAPE
      if(sigbkgextrafactordir_.size()==0){
	sigbkgshape= filemanager->GetSetsShape(sigbkgset_,shape_[iShape],basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_,false);
      }
      //WEIGHT INDIVIDUAL SHAPES BY DATA DRIVEN WEIGHTS
      else{
	if(sigbkgextrafactordir_.size()==sigbkgset_.size()){
	  for(unsigned iBkg=0;iBkg<sigbkgset_.size();iBkg++){
	    double extrafactor=1;
	    TDirectory* extrafactordir;
	    TVectorD *ddweight;
	    if(sigbkgextrafactordir_[iBkg]==""){
	      extrafactor=1;
	    }
	    else if(!fs_->GetDirectory(sigbkgextrafactordir_[iBkg].c_str())){
	      std::cout<<"Directory "<<sigbkgextrafactordir_[iBkg]<<" does not exist, exiting"<<std::endl;
	      return 1;
	    }
	    else{
	      extrafactordir=fs_->GetDirectory(sigbkgextrafactordir_[iBkg].c_str());
	      ddweight = (TVectorD*)extrafactordir->Get("ddweight");
	      if(sigbkgisz_.size()==sigbkgset_.size()){
		if(sigbkgisz_[iBkg]==0){
		  extrafactor=(*ddweight)[0];
		}
		if(sigbkgisz_[iBkg]==1){
		  extrafactor=(*ddweight)[0];//EWK WEIGHT
		}
		if(sigbkgisz_[iBkg]==2){
		  extrafactor=(*ddweight)[1];//QCD WEIGHT
		}
	      }
	    }
	    TH1F nextbkg;
	    //SPECIAL CASE FOR Z BKG
	    if(sigbkgisz_.size()==sigbkgset_.size()){
	      if(sigbkgisz_[iBkg]!=0){
		//GET Z SHAPE AND WEIGHT IT
		nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],shape_[iShape],basesel_,zsigcat_,"weight_nolep*"+boost::lexical_cast<std::string>(extrafactor),false);
	      }
	      else nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],shape_[iShape],basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	    }
	    else nextbkg=filemanager->GetSetShape(sigbkgset_[iBkg],shape_[iShape],basesel_,(sigcat_+sigbkgextrasel_),sigqcdweight_+"*"+boost::lexical_cast<std::string>(extrafactor),false);
	    
	    //ADD HISTOGRAMS TOGETHER
	    if(firstbkg){
	      sigbkgshape=nextbkg;
	      firstbkg=false;
	    }
	    else{
	      sigbkgshape.Add(&nextbkg);
	    }
	  }
	}
	else{
	  std::cout<<"Extra factor dir vector must be same size as bkg set vector expect errors!"<<std::endl;
	  return 1;
	}
      }



      sigqcdshape.Add(&sigbkgshape,-1.);
      dir->cd();
      std::string histname;
      if(shapename_.size()==0){
	std::vector<std::string> strs;
	boost::split(strs, shape_[iShape], boost::is_any_of("("));
	histname=strs[0];
      }
      else{
	histname=shapename_[iShape];
      }
      sigqcdshape.SetName(histname.c_str());
      if(iShape==0) std::cout << " Raw NSMC: " << Integral(&sigqcdshape) << "+-" << Error(&sigqcdshape) << std::endl;

      sigqcdshape.Scale(weight*sigcontextrafactor_);
      double nsqcd=Integral(&sigqcdshape);
      double nsqcderr=Error(&sigqcdshape);
      double weightednsmcstatfrac=sqrt(((nsqcderr/nsqcd)*(nsqcderr/nsqcd))+(weighterrmcstatfrac*weighterrmcstatfrac));
      //!!MAKE SURE  TO TAKE NSQCD ERROR OUT OF NORMALISATION ERROR AS SHOULD BE DONE BIN BY BIN
      if(iShape==0){
	std::cout<<"  Normalised NSQCD: "<<nsqcd<<"+-"<<nsqcd*weighterrdatastatfrac<<"(data stat.)+-"<<nsqcd*weightednsmcstatfrac<<"(MC stat.)"<<std::endl;
	errvec[0]=weighterrdatastatfrac;
	errvec[1]=weightednsmcstatfrac;
      }
      sigqcdshape.Write();
    }
    dir->cd();
    errvec.Write("normerrs");
    weightvec.Write("ddweight");
  
    return 0;
  };

}
