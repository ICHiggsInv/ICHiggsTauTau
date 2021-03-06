#ifndef ICHiggsTauTau_HiggsNuNu_HistPlotter_h                                                                       
#define ICHiggsTauTau_HiggsNuNu_HistPlotter_h
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include <string>
#include <vector>
#include <utility>
#include "TH1F.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"


namespace ic {

  class LTPlotElement{
  public:
    LTPlotElement();
    ~LTPlotElement();
    void ApplyStyle();
  protected:
    CLASS_MEMBER(LTPlotElement,int,color)
    CLASS_MEMBER(LTPlotElement,int,marker_color)
    CLASS_MEMBER(LTPlotElement,int,line_color)
    CLASS_MEMBER(LTPlotElement,int,fill_color)
    CLASS_MEMBER(LTPlotElement,int,fill_style)
    CLASS_MEMBER(LTPlotElement,int,line_style)
    CLASS_MEMBER(LTPlotElement,int,marker_style)
    CLASS_MEMBER(LTPlotElement,bool,draw_fill)
    CLASS_MEMBER(LTPlotElement,bool,draw_fill_in_legend)
    CLASS_MEMBER(LTPlotElement,bool,draw_marker)    
    CLASS_MEMBER(LTPlotElement,bool,draw_line)    
    CLASS_MEMBER(LTPlotElement,bool,draw_stat_error_y)    
    CLASS_MEMBER(LTPlotElement,int,line_width)    
    CLASS_MEMBER(LTPlotElement,double,marker_size)    
    CLASS_MEMBER(LTPlotElement,std::string,legname)
    CLASS_MEMBER(LTPlotElement,std::string,unit)
    CLASS_MEMBER(LTPlotElement,double,scale)
    CLASS_MEMBER(LTPlotElement,bool,in_stack)//NOTE DATA CANNOT BE STACKED AT THE MOMENT
    CLASS_MEMBER(LTPlotElement,bool,is_data)
    CLASS_MEMBER(LTPlotElement,bool,is_inrationum)
    CLASS_MEMBER(LTPlotElement,bool,is_inratioden)
    CLASS_MEMBER(LTPlotElement,int,has_dderrors)
    CLASS_MEMBER(LTPlotElement,std::string,drawopts)
    CLASS_MEMBER(LTPlotElement,std::string,legopts)
    CLASS_MEMBER(LTPlotElement,std::string,sample)
    CLASS_MEMBER(LTPlotElement,TH1F*,hist_ptr)
    CLASS_MEMBER(LTPlotElement,std::vector<std::string>,blindvar)
    std::vector<std::pair<double,double> > blindrange_;
    LTPlotElement set_blindrange(std::vector<std::pair<double,double> > const& blindrange) {blindrange_ = blindrange; return *this; }
    std::vector<std::pair<double,double> > blindrange() {return blindrange_; }
    
  };

  class LTShapeElement{
  public:
    LTShapeElement();
    ~LTShapeElement();
  private:
    CLASS_MEMBER(LTShapeElement,std::string,name)
    CLASS_MEMBER(LTShapeElement,std::string,histtitle)
    CLASS_MEMBER(LTShapeElement,bool,dology)
    CLASS_MEMBER(LTShapeElement,double,axisrangemultiplier)
    CLASS_MEMBER(LTShapeElement,double,legleft)
    CLASS_MEMBER(LTShapeElement,double,legright)
  };

  class HistPlotter : public LTModule{ 
    CLASS_MEMBER(HistPlotter,std::string,dirname)
    CLASS_MEMBER(HistPlotter,std::vector<LTPlotElement>,elements)   
    CLASS_MEMBER(HistPlotter,std::vector<LTShapeElement>,shapes)
      //CLASS_MEMBER(HistPlotter,std::vector<std::string>,shapes)   
      //CLASS_MEMBER(HistPlotter,std::vector<std::string>,histTitles)   
    CLASS_MEMBER(HistPlotter,bool,do_ratio)
    CLASS_MEMBER(HistPlotter,bool,do_ratio_line)
    CLASS_MEMBER(HistPlotter,bool,do_ratio_fitline)
    CLASS_MEMBER(HistPlotter,bool,add_underflows)
    CLASS_MEMBER(HistPlotter,bool,add_overflows)

  public:	
    HistPlotter(std::string);
    virtual ~HistPlotter();
    virtual int Init(TFile*);
    virtual int Run(LTFiles*);
    static void SetMCStackStyle(LTPlotElement*);
    static void SetSignalStyle(LTPlotElement*);
    static void SetDataStyle(LTPlotElement*);
  };

}
#endif
