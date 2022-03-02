#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"

#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"


class HOTVRJetCorrectionModule : public uhh2::AnalysisModule {
 public:
  HOTVRJetCorrectionModule(uhh2::Context & ctx);
  virtual bool process(uhh2::Event &event) override;
 private:
  void set_subjet_handles(uhh2::Event &event);
  void rebuild_jets(uhh2::Event &event);
  bool is_mc;
  std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;
  std::shared_ptr<RunSwitcher> jec_switcher_UL16preVFP,jec_switcher_UL16postVFP, jec_switcher_UL17, jec_switcher_UL18;
  std::unique_ptr<GenericJetResolutionSmearer> jer_module;
  uhh2::Event::Handle<std::vector<Jet> > h_subj;
  uhh2::Event::Handle<std::vector<std::pair<int, int> > > h_subjmap;
  uhh2::Event::Handle<std::vector<GenJet> > h_gensubj;
};
