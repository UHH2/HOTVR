#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"

class HOTVRJetCorrectionModule : public uhh2::AnalysisModule {
 public:
  HOTVRJetCorrectionModule(uhh2::Context & ctx);
  virtual bool process(uhh2::Event &event) override;
 private:
  bool is_mc;
  std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;
  std::unique_ptr<AnalysisModule> jlc_module;
};
