#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/HOTVR/include/HadronicTop.h"

#include <vector>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
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
  void get_sf(double jetPt, int category);

  const TopJetId m_id_topjet;
  const TString m_sys_direction;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
  uhh2::Event::Handle<std::vector<GenTop>> h_tophad;

  double m_weight = 1.;
  double m_weight_up = 1.;
  double m_weight_down = 1.;
  double m_weight_merged_up = 1.;
  double m_weight_merged_down = 1.;
  double m_weight_semi_up = 1.;
  double m_weight_semi_down = 1.;
  double m_weight_non_up = 1.;
  double m_weight_non_down = 1.;

  TGraphAsymmErrors *graph_merged_tot, *graph_merged_stat, *graph_merged_syst, *graph_non_tot, *graph_non_stat, *graph_non_syst;

  uhh2::Event::Handle<double> h_toptag_weight;
  uhh2::Event::Handle<double> h_toptag_weight_up, h_toptag_weight_merged_up, h_toptag_weight_semi_up, h_toptag_weight_non_up;
  uhh2::Event::Handle<double> h_toptag_weight_down, h_toptag_weight_merged_down, h_toptag_weight_semi_down, h_toptag_weight_non_down;

};
