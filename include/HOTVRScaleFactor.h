#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/HOTVR/include/HadronicTop.h"

#include <vector>
#include <TH1F.h>
#include <TFile.h>

/** \brief Top tagging scale factors for HOTVR
 */
class HOTVRScaleFactor : public uhh2::AnalysisModule {
 public:
  explicit HOTVRScaleFactor(uhh2::Context &ctx, 
			    const TopJetId &id_topjet,
			    const std::string &sys_direction = "nominal",
			    const std::string &gen_handle_name = "HadronicTop",
			    const std::string &param_name = "TopTagSF",
			    const std::string &xmlpathname = "HOTVRTopTagSFs");

  virtual bool process(uhh2::Event &event) override;

 private:
  void get_sf(double pt, int category);

  const TopJetId m_id_topjet;
  const TString m_sys_direction;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
  uhh2::Event::Handle<std::vector<GenTop>> h_tophad;

  double m_weight = 1.;
  double m_weight_up = 1.;
  double m_weight_down = 1.;

  TH1F *sf_merged, *sf_merged_up, *sf_merged_down;
  TH1F *sf_semi, *sf_semi_up, *sf_semi_down;
  TH1F *sf_not, *sf_not_up, *sf_not_down;

  uhh2::Event::Handle<double> h_toptag_weight;
  uhh2::Event::Handle<double> h_toptag_weight_up;
  uhh2::Event::Handle<double> h_toptag_weight_down;

};
