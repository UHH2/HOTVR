#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"


using namespace uhh2;
using namespace std;


HOTVRScaleFactor::HOTVRScaleFactor(uhh2::Context &ctx, const TopJetId &id_topjet, const string &sys_direction, const string &gen_handle_name, const string &param_name, const string &xmlpathname):
  m_id_topjet(id_topjet),
  m_sys_direction(sys_direction),
  h_tophad(ctx.get_handle<vector<GenTop> >(gen_handle_name))
{
  h_toptag_weight = ctx.declare_event_output<double>(param_name);
  h_toptag_weight_up = ctx.declare_event_output<double>(param_name+"_up");
  h_toptag_weight_down = ctx.declare_event_output<double>(param_name+"_down");
  h_toptag_weight_merged_up = ctx.declare_event_output<double>(param_name+"_merged_up");
  h_toptag_weight_merged_down = ctx.declare_event_output<double>(param_name+"_merged_down");
  h_toptag_weight_non_up = ctx.declare_event_output<double>(param_name+"_non_up");
  h_toptag_weight_non_down = ctx.declare_event_output<double>(param_name+"_non_down");
  // open file and get scale factor histograms
  string path = ctx.get(xmlpathname);
  TFile *f = TFile::Open(path.c_str(), "READ");

  //For merged category
  TString graph_basename_merged = "HOTVR_PUPPI_default/FullyMerged_";
  graph_merged_tot  = (TGraphAsymmErrors*)f->Get((graph_basename_merged+"tot").Data());
  graph_merged_stat = (TGraphAsymmErrors*)f->Get((graph_basename_merged+"stat").Data());
  graph_merged_syst = (TGraphAsymmErrors*)f->Get((graph_basename_merged+"syst").Data());

  //For non-merged category
  TString graph_basename_non = "HOTVR_PUPPI_default/NotMerged_";
  graph_non_tot  = (TGraphAsymmErrors*)f->Get((graph_basename_non+"tot").Data());
  graph_non_stat = (TGraphAsymmErrors*)f->Get((graph_basename_non+"stat").Data());
  graph_non_syst = (TGraphAsymmErrors*)f->Get((graph_basename_non+"syst").Data());

  f->Close();
}

void HOTVRScaleFactor::get_sf(double jetPt, int category) {

//merged category
if (category == 3){

  UInt_t ibin = 0;
  Double_t x;
  Double_t y;
  graph_merged_tot->GetPoint(ibin, x, y);

  const Double_t minimum_jetPt = x - graph_merged_tot->GetErrorXlow(ibin);
  if(jetPt < minimum_jetPt) {
    std::cerr << "No scale factor for such low jet pT available. Will set scale factor to 1 with zero uncertainty" << std::endl;
  }

  while(!(jetPt >= x - graph_merged_tot->GetErrorXlow(ibin) && jetPt < x + graph_merged_tot->GetErrorXhigh(ibin))) { // check if jetPt within pT interval of bin with bin number "ibin"
    ++ibin;
    graph_merged_tot->GetPoint(ibin, x, y);
    if(ibin + 1 >= (UInt_t)graph_merged_tot->GetN()) break; // use scale factor of last bin if given jetPt exceeds highest bin edge of the graph
  }

  m_weight *= y;
  m_weight_up *= graph_merged_tot->GetErrorYhigh(ibin);
  m_weight_down *= graph_merged_tot->GetErrorYlow(ibin);
  m_weight_merged_up *= graph_merged_tot->GetErrorYhigh(ibin);
  m_weight_merged_down *= graph_merged_tot->GetErrorYlow(ibin);
  m_weight_non_up *= graph_merged_tot->GetErrorYhigh(ibin);
  m_weight_non_down *= graph_merged_tot->GetErrorYlow(ibin);

}
//non-merged category
else
{

  UInt_t ibin = 0;
  Double_t x;
  Double_t y;
  graph_non_tot->GetPoint(ibin, x, y);

  const Double_t minimum_jetPt = x - graph_non_tot->GetErrorXlow(ibin);
  if(jetPt < minimum_jetPt) {
    std::cerr << "No scale factor for such low jet pT available. Will set scale factor to 1 with zero uncertainty" << std::endl;
  }

  while(!(jetPt >= x - graph_non_tot->GetErrorXlow(ibin) && jetPt < x + graph_non_tot->GetErrorXhigh(ibin))) { // check if jetPt within pT interval of bin with bin number "ibin"
    ++ibin;
    graph_non_tot->GetPoint(ibin, x, y);
    if(ibin + 1 >= (UInt_t)graph_non_tot->GetN()) break; // use scale factor of last bin if given jetPt exceeds highest bin edge of the graph
  }

  m_weight *= y;
  m_weight_up *= graph_non_tot->GetErrorYhigh(ibin);
  m_weight_down *= graph_non_tot->GetErrorYlow(ibin);
  m_weight_merged_up *= graph_non_tot->GetErrorYhigh(ibin);
  m_weight_merged_down *= graph_non_tot->GetErrorYlow(ibin);
  m_weight_non_up *= graph_non_tot->GetErrorYhigh(ibin);
  m_weight_non_down *= graph_non_tot->GetErrorYlow(ibin);

}

}

bool HOTVRScaleFactor::process(Event &event) {

  m_weight = 1.;
  m_weight_up = 1.;
  m_weight_down = 1.;
  m_weight_merged_up = 1.;
  m_weight_merged_down = 1.;
  m_weight_non_up = 1.;
  m_weight_non_down = 1.;
  if (event.isRealData || event.get(h_tophad).size() == 0)
    {
      event.set(h_toptag_weight, 1.);
      event.set(h_toptag_weight_up, 1.);
      event.set(h_toptag_weight_down, 1.);
      event.set(h_toptag_weight_merged_up, 1.);
      event.set(h_toptag_weight_merged_down, 1.);
      event.set(h_toptag_weight_non_up, 1.);
      event.set(h_toptag_weight_non_down, 1.);
      return false;
    }

  vector<GenTop> gentops = event.get(h_tophad);
  for (const auto &topjet : *event.topjets)
    {
      int nMatched = 0;
      bool bMatched = false;
      bool q1Matched = false;
      bool q2Matched = false;
      if (m_id_topjet(topjet, event))
	{
	  double dRmatch = min(1.5, max(0.1, 600.0 / (topjet.pt() * topjet.JEC_factor_raw()) )); // calculate distance using clustering distance parameter
	  for (auto top : gentops)
	    {
	      if (deltaR(top.get_b(), topjet.v4()) < dRmatch) bMatched = true;
	      if (deltaR(top.get_q1(), topjet.v4()) < dRmatch) q1Matched = true;
	      if (deltaR(top.get_q2(), topjet.v4()) < dRmatch) q2Matched = true;
	    }
	  if (bMatched) ++nMatched;
	  if (q1Matched) ++nMatched;
	  if (q2Matched) ++nMatched;
	  get_sf(topjet.pt(), nMatched);
	}
    }


  event.set(h_toptag_weight, m_weight);
  event.set(h_toptag_weight_up, m_weight_up);
  event.set(h_toptag_weight_down, m_weight_down);
  event.set(h_toptag_weight_merged_up, m_weight_merged_up);
  event.set(h_toptag_weight_merged_down, m_weight_merged_down);
  event.set(h_toptag_weight_non_up, m_weight_non_up);
  event.set(h_toptag_weight_non_down, m_weight_non_down);

  if (m_sys_direction == "nominal") event.weight *= m_weight;
  else if (m_sys_direction == "up") event.weight *= m_weight_up;
  else if (m_sys_direction == "down") event.weight *= m_weight_down;
  else if (m_sys_direction == "up_merged") event.weight *= m_weight_merged_up;
  else if (m_sys_direction == "down_merged") event.weight *= m_weight_merged_down;
  else if (m_sys_direction == "up_non") event.weight *= m_weight_non_up;
  else if (m_sys_direction == "down_non") event.weight *= m_weight_non_down;
  return true;
}
