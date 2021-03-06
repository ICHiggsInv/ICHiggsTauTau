#ifndef ICHiggsTauTau_CombineTools_CombineHarvester_h
#define ICHiggsTauTau_CombineTools_CombineHarvester_h
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <map>
#include <list>
#include <cmath>
#include <set>
#include <functional>
#include "TFile.h"
#include "TH1.h"
#include "RooWorkspace.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/HelperFunctions.h"

namespace ch {

class CombineHarvester {
 public:
  // Low-level
  CombineHarvester();
  ~CombineHarvester();
  CombineHarvester(CombineHarvester const& other);
  CombineHarvester(CombineHarvester&& other);
  CombineHarvester& operator=(CombineHarvester other);
  CombineHarvester cp();
  CombineHarvester& PrintAll();

  /**
   * @name Datacards
   * @brief Methods for the reading and writing of datacards.
   * @details As well as reading or writing the plain-text datacard
   * files these methods also handle the automatic loading and saving
   * of any TH1 objects or RooWorkspaces.
   */
  /**@{*/
  int ParseDatacard(std::string const& filename,
      std::string const& analysis,
      std::string const& era,
      std::string const& channel,
      int bin_id,
      std::string const& mass);
  int ParseDatacard(std::string const& filename,
      std::string parse_rule);

  void WriteDatacard(std::string const& name, std::string const& root_file);
  void WriteDatacard(std::string const& name, TFile & root_file);
  /**@}*/

  /**
   * @name Filters
   * @brief A collection of methods to filter the Observation, Process and Nuisance objects
   * @details Each method here will remove objects from the Observation, Process and Nuisance
   * collections where the value of the relevant property is not equal to one of those specified
   * in the list vec. Each methods returns a reference to the class instance, meaning
   * it is simple to chain multiple filters. Note that filters only act on a set of objects
   * where that property is defined. For example, the method nus_name will only filter the Nuisance
   * collection, and process will only filter the Process and Nuisance collections.
   *
   * @param vec A vector containing the possible values of the given property.
   */
  /**@{*/
  CombineHarvester& bin(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& bin_id(std::vector<int> const& vec, bool cond = true);
  CombineHarvester& process(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& analysis(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& era(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& channel(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& mass(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& nus_name(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& nus_type(std::vector<std::string> const& vec, bool cond = true);

  CombineHarvester& signals();
  CombineHarvester& backgrounds();
  CombineHarvester& histograms();
  CombineHarvester& pdfs();
  CombineHarvester& data();
  /**@}*/


  // Set generation
  template<typename T>
  std::set<T> GenerateSetFromProcs(std::function<T(ch::Process const*)> func);

  template<typename T>
  std::set<T> GenerateSetFromObs(std::function<T(ch::Observation const*)> func);

  template<typename T>
  std::set<T> GenerateSetFromNus(std::function<T(ch::Nuisance const*)> func);

  /**
   * @name Modification
   * @brief Methods to modify existing objects.
   * @details See the documentation of each method for details
   */
  /**@{*/
  void SetParameters(std::vector<ch::Parameter> params);
  void UpdateParameters(std::vector<ch::Parameter> params);
  std::vector<ch::Parameter> GetParameters() const;

  template<typename Function>
  void ForEachProc(Function func);

  template<typename Function>
  void ForEachObs(Function func);

  template<typename Function>
  void ForEachNus(Function func);

  void VariableRebin(std::vector<double> bins);
  /**@}*/

  /**
   * @name Rates and shapes
   * @brief Methods to calculate total yields and shapes
   * @details All of these methods are greedy, meaning they will sum
   * over all available objects and evaluate the effect of all uncertainties.
   * They should be used at the end of a chain of filter methods to give the
   * desired yield, shape or uncertainty.
   */
  /**@{*/
  double GetRate();
  double GetObservedRate();
  double GetUncertainty();
  double GetUncertainty(RooFitResult const* fit, unsigned n_samples);
  TH1F GetShape();
  TH1F GetShapeWithUncertainty();
  TH1F GetShapeWithUncertainty(RooFitResult const* fit, unsigned n_samples);
  TH1F GetObservedShape();
  /**@}*/

  /**
   * @name Creating new entries
   * @brief Methods to create new Observation, Process, Nuisance and Parameter objects
   * @details The general order in which these are run is given below.
   *
   *    1.  The CombineHarvester::AddObservations and CombineHarvester::AddProcesses methods
   *        build the desired ch::Observation and ch::Process entries. The rate entries will
   *        default to zero at this stage, and are determined automatically in step 3.
   *    2.  The CombineHarvester::AddSyst method will iterate through each ch::Process and
   *        add a matching ch::Nuisance entry, with the name, type and value specfied in the
   *        method arguments. This method can be used at the end of a chain of filters to only
   *        create ch::Nuisance objects for a subset of processes.
   *    3.  The method CombineHarvester::ExtractShapes opens a specified ROOT file and uses the
   *        provided pattern rules to load the TH1 objects for all entries. The TH1 integrals
   *        are used to set the event rates at this stage.
   *    4.  Optionally, methods such as  CombineHarvester::AddBinByBin may be used to prepare
   *        additional entries.
   */
  /**@{*/
  void AddObservations(std::vector<std::string> mass,
                       std::vector<std::string> analysis,
                       std::vector<std::string> era,
                       std::vector<std::string> channel,
                       std::vector<std::pair<int, std::string>> bin);

  void AddProcesses(std::vector<std::string> mass,
                    std::vector<std::string> analysis,
                    std::vector<std::string> era,
                    std::vector<std::string> channel,
                    std::vector<std::string> procs,
                    std::vector<std::pair<int, std::string>> bin, bool signal);

  template <class Map>
  void AddSyst(CombineHarvester* other, std::string const& name,
               std::string const& type, Map const& valmap);

  void ExtractShapes(std::string const& file, std::string const& rule,
                     std::string const& syst_rule);

  void AddWorkspace(RooWorkspace const* ws);
  void ExtractPdfs(std::string const& ws_name, std::string const& rule);
  void ExtractData(std::string const& ws_name, std::string const& rule);
  void AddBinByBin(double threshold, bool fixed_norm, CombineHarvester* other);
  void MergeBinErrors(double bbb_threshold, double merge_threshold);
  /**@}*/

 private:
  std::vector<std::shared_ptr<Observation>> obs_;
  std::vector<std::shared_ptr<Process>> procs_;
  std::vector<std::shared_ptr<Nuisance>> nus_;
  std::map<std::string, std::shared_ptr<Parameter>> params_;
  std::map<std::string, std::shared_ptr<RooWorkspace>> wspaces_;

  struct HistMapping {
    std::string process;
    std::string category;
    std::shared_ptr<TFile> file;
    std::string pattern;
    std::string syst_pattern;
    bool IsHist();
    bool IsPdf();
    bool IsData();
  };

  typedef std::vector<std::pair<std::string, std::string>> StrPairVec;
  typedef std::vector<std::string> StrVec;
  StrPairVec GenerateShapeMapAttempts(std::string process,
      std::string category);
  HistMapping ResolveMapping(std::vector<HistMapping> const& mappings,
      std::string const& bin,
      std::string const& process,
      std::string const& mass,
      std::string const& nuisance);
  TH1 * GetHistFromFile(HistMapping const& mapping,
      unsigned type);
  RooAbsPdf * GetPdfFromFile(HistMapping const& mapping);
  RooAbsData * GetDataFromFile(HistMapping const& mapping);

  void WriteHistToFile(
      TH1 const* hist,
      TFile * file,
      std::vector<HistMapping> const& mappings,
      std::string const& bin,
      std::string const& process,
      std::string const& mass,
      std::string const& nuisance,
      unsigned type);

  friend void swap(CombineHarvester& first, CombineHarvester& second);
  typedef std::vector<std::vector<Nuisance const*>>
    ProcNusMap;
  ProcNusMap GenerateProcNusMap();

  double GetRateInternal(ProcNusMap const& lookup,
    std::string const& single_nus = "");

  TH1F GetShapeInternal(ProcNusMap const& lookup,
    std::string const& single_nus = "");

  inline double smoothStepFunc(double x) const {
    if (std::fabs(x) >= 1.0/*_smoothRegion*/) return x > 0 ? +1 : -1;
    double xnorm = x/1.0;/*_smoothRegion*/
    double xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }

  TH1F ShapeDiff(double x, TH1 const* nom, TH1 const* low, TH1 const* high);

  void CreateParameterIfEmpty(CombineHarvester *cmb, std::string const& name);
};

template<typename T>
std::set<T> CombineHarvester::GenerateSetFromProcs(
    std::function<T (ch::Process const*)> func) {
  std::set<T> ret;
  for (auto const& item : procs_) ret.insert(func(item.get()));
  return ret;
}

template<typename T>
std::set<T> CombineHarvester::GenerateSetFromObs(
    std::function<T (ch::Observation const*)> func) {
  std::set<T> ret;
  for (auto const& item : obs_) ret.insert(func(item.get()));
  return ret;
}

template<typename T>
std::set<T> CombineHarvester::GenerateSetFromNus(
    std::function<T (ch::Nuisance const*)> func) {
  std::set<T> ret;
  for (auto const& item : nus_) ret.insert(func(item.get()));
  return ret;
}

template<typename Function>
void CombineHarvester::ForEachProc(Function func) {
  for (auto & item: procs_) func(item.get());
}

template<typename Function>
void CombineHarvester::ForEachObs(Function func) {
  for (auto & item: obs_) func(item.get());
}

template<typename Function>
void CombineHarvester::ForEachNus(Function func) {
  for (auto & item: nus_) func(item.get());
}

template <class Map>
void CombineHarvester::AddSyst(CombineHarvester* other, std::string const& name,
                               std::string const& type,
                               Map const& valmap) {
  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (!valmap.Contains(procs_[i].get())) continue;
    std::string subbed_name = name;
    boost::replace_all(subbed_name, "$BIN", procs_[i]->bin());
    boost::replace_all(subbed_name, "$PROCESS", procs_[i]->process());
    boost::replace_all(subbed_name, "$MASS", procs_[i]->mass());
    boost::replace_all(subbed_name, "$ERA", procs_[i]->era());
    boost::replace_all(subbed_name, "$CHANNEL", procs_[i]->channel());
    boost::replace_all(subbed_name, "$ANALYSIS", procs_[i]->analysis());
    auto nus = std::make_shared<Nuisance>();
    ch::SetProperties(nus.get(), procs_[i].get());
    nus->set_name(subbed_name);
    nus->set_type(type);
    if (type == "lnN") {
      nus->set_asymm(false);
      nus->set_value_u(valmap.Val(procs_[i].get()));
    } else if (type == "shape") {
      nus->set_asymm(true);
      nus->set_value_u(1.0);
      nus->set_value_d(1.0);
    }
    CombineHarvester::CreateParameterIfEmpty(other ? other : this,
                                             nus->name());
    if (other) {
      other->nus_.push_back(nus);
    } else {
      nus_.push_back(nus);
    }
  }
}
}

#endif
