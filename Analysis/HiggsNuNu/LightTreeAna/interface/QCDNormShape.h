#ifndef ICHiggsTauTau_HiggsNuNu_QCDNormShape_h                                                                       
#define ICHiggsTauTau_HiggsNuNu_QCDNormShape_h
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include <string>

namespace ic {

  class QCDNormShape : public LTModule{
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,shape)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,shapename)
    CLASS_MEMBER(QCDNormShape,std::string,sigqcdset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,sigbkgset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,sigbkgextrafactordir)
    CLASS_MEMBER(QCDNormShape,std::vector<int>,sigbkgisz)
    CLASS_MEMBER(QCDNormShape,std::string,contqcdset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,contqcdbkgset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,contqcdbkgextrafactordir)
    CLASS_MEMBER(QCDNormShape,std::vector<int>,contqcdbkgisz)
    CLASS_MEMBER(QCDNormShape,std::string,contdataset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,contbkgset)
    CLASS_MEMBER(QCDNormShape,std::vector<std::string>,contbkgextrafactordir)
    CLASS_MEMBER(QCDNormShape,std::vector<int>,contbkgisz)
    CLASS_MEMBER(QCDNormShape,std::string,sigcat)
    CLASS_MEMBER(QCDNormShape,std::string,zsigcat)
    CLASS_MEMBER(QCDNormShape,std::string,zcontcat)
    CLASS_MEMBER(QCDNormShape,std::string,zcontqcdcat)
    CLASS_MEMBER(QCDNormShape,std::string,contcat)
    CLASS_MEMBER(QCDNormShape,std::string,contdataextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,contqcdextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,contqcdbkgextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,contbkgextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,sigbkgextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,sigqcdextrasel)
    CLASS_MEMBER(QCDNormShape,std::string,basesel)
    CLASS_MEMBER(QCDNormShape,std::string,sigqcdweight)
    CLASS_MEMBER(QCDNormShape,std::string,contqcdweight)
    CLASS_MEMBER(QCDNormShape,std::string,contdataweight)
    CLASS_MEMBER(QCDNormShape,double,sigcontextrafactor)
    CLASS_MEMBER(QCDNormShape,std::string,dirname)
  public:
    virtual QCDNormShape & set_shape(std::string const& shape) {
      std::vector<std::string> shapes;
      shapes.push_back(shape);
      shape_ = shapes;
      return *this;
    }
    QCDNormShape(std::string);
    virtual ~QCDNormShape();
    virtual int Init(TFile*);
    virtual int Run(LTFiles*);
  };

}
#endif
