#include "CombineTools/interface/HelperFunctions.h"
#include "Utilities/interface/FnRootTools.h"
#include <iostream>
#include <vector>
#include "RooFitResult.h"
#include "RooRealVar.h"

namespace ch {

std::vector<ch::Parameter> ExtractFitParameters(RooFitResult const& res) {
  std::vector<ch::Parameter> params;
  params.resize(res.floatParsFinal().getSize());
  for (int i = 0; i < res.floatParsFinal().getSize(); ++i) {
    RooRealVar const* var = dynamic_cast<RooRealVar const*>(res.floatParsFinal().at(i));
    params[i].set_name(std::string(var->GetName()));
    params[i].set_val(var->getVal());
    params[i].set_err_d(var->getErrorLo());
    params[i].set_err_u(var->getErrorHi());
  }
  return params;
}

std::vector<ch::Parameter> ExtractSampledFitParameters(RooFitResult const& res) {
  std::vector<ch::Parameter> params;
  params.resize(res.floatParsFinal().getSize());
  RooArgList const& rands = res.randomizePars();
  for (int i = 0; i < res.floatParsFinal().getSize(); ++i) {
    RooRealVar const* var = dynamic_cast<RooRealVar const*>(rands.at(i));
    params[i].set_name(std::string(var->GetName()));
    params[i].set_val(var->getVal());
    params[i].set_err_d(var->getErrorLo());
    params[i].set_err_u(var->getErrorHi());
  }
  return params;
}

SOverBInfo::SOverBInfo(TH1F const* sig, TH1F const* bkg, unsigned steps, double frac) {
  double xmin = sig->GetXaxis()->GetXmin();
  double xmax = sig->GetXaxis()->GetXmax();
  double step_size = (xmax-xmin)/double(steps);
  double sig_tot = sig->Integral();
  double lower_limit = 0;
  double upper_limit = 0;
  double ofrac = (1.-frac)/2.;
  for (unsigned j = 0; j < steps; ++j) {
    double integral = ch::IntegrateFloatRange(sig, xmin, xmin+(step_size*double(j)));
    if (integral/sig_tot > ofrac) {
      lower_limit = xmin+(step_size*double(j));
      break;
    }
  }
  for (unsigned j = 0; j < steps; ++j) {
    double integral = ch::IntegrateFloatRange(sig, xmax - (step_size*double(j)), xmax);
    if (integral/sig_tot > ofrac) {
      upper_limit = xmax - (step_size*double(j));
      break;
    }
  }
  x_lo = lower_limit;
  x_hi = upper_limit;
  s = ch::IntegrateFloatRange(sig, lower_limit, upper_limit);
  b = ch::IntegrateFloatRange(bkg, lower_limit, upper_limit);
}

double IntegrateFloatRange(TH1F const* hist, double xmin, double xmax) {
    TAxis *axis = hist->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = hist->Integral(bmin, bmax);
    integral -= hist->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/
              axis->GetBinWidth(bmin);
    integral -= hist->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/
              axis->GetBinWidth(bmax);
    return integral;
}

void ParseTable(std::map<std::string, TGraph>* graphs,
                        std::string const& file,
                        std::vector<std::string> const& fields) {
  auto lines = ic::ParseFileLines(file);
  for (auto const& f : fields) {
    (*graphs)[f] = TGraph(lines.size());
  }
  for (unsigned i = 0; i < lines.size(); ++i) {
    std::vector<std::string> words;
    boost::split(words, lines[i], boost::is_any_of("\t "),
        boost::token_compress_on);
    // std::cout << "\"" << words.at(0) << "\"\n";
    double m = boost::lexical_cast<double>(words.at(0));
    for (unsigned f = 0; f < fields.size(); ++f) {
    // std::cout << "\"" << words.at(f+1) << "\"\n";
      (*graphs)[fields[f]]
          .SetPoint(i, m, boost::lexical_cast<double>(words.at(f+1)));
    }
  }
}

void ScaleProcessRate(ch::Process* p,
                      std::map<std::string, TGraph> const* graphs,
                      std::string const& prod, std::string const& decay) {
  double mass = boost::lexical_cast<double>(p->mass());
  double scale = 1.0;
  if (prod != "" && graphs->count(prod)) {
    scale *= graphs->find(prod)->second.Eval(mass);
  }
  if (decay != "" && graphs->count(decay)) {
    scale *= graphs->find(decay)->second.Eval(mass);
  }
  p->set_rate(p->rate() * scale);
}

std::vector<std::string> JoinStr(std::vector<vector<std::string>> const& in) {
  return Join<std::string>(in);
}
}
