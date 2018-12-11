#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"


using namespace uhh2;
using namespace std;

HOTVRScaleFactor::HOTVRScaleFactor(uhh2::Context &ctx, const TopJetId &id_topjet, const string &sys_direction, const string &gen_handle_name, const string &param_name, const string &xmlpathname):
  m_id_topjet(id_topjet),
  m_sys_direction(sys_direction),
  h_tophad(ctx.get_handle<vector<GenTop> >(gen_handle_name)),
  h_toptag_weight(ctx.declare_event_output<float>(param_name)),
  h_toptag_weight_up(ctx.declare_event_output<float>(param_name+"_up")),
  h_toptag_weight_down(ctx.declare_event_output<float>(param_name+"_down"))
{
  // open file and get scale factor histograms
  string path = ctx.get(xmlpathname); 
  TFile *f = new TFile(path.c_str(), "READ");
  sf_merged = (TH1F*) f->Get("HOTVR/sf_mergedTop_nominal");
  sf_merged->SetDirectory(0);
  sf_merged_up = (TH1F*) f->Get("HOTVR/sf_mergedTop_up");
  sf_merged_up->SetDirectory(0);
  sf_merged_down = (TH1F*) f->Get("HOTVR/sf_mergedTop_down");
  sf_merged_down->SetDirectory(0);

  sf_semi = (TH1F*) f->Get("HOTVR/sf_semimerged_nominal");
  sf_semi->SetDirectory(0);
  sf_semi_up = (TH1F*) f->Get("HOTVR/sf_semimerged_up");
  sf_semi_up->SetDirectory(0);
  sf_semi_down = (TH1F*) f->Get("HOTVR/sf_semimerged_down");
  sf_semi_down->SetDirectory(0);

  sf_not = (TH1F*) f->Get("HOTVR/sf_notmerged_nominal");
  sf_not->SetDirectory(0);
  sf_not_up = (TH1F*) f->Get("HOTVR/sf_notmerged_up");
  sf_not_up->SetDirectory(0);
  sf_not_down = (TH1F*) f->Get("HOTVR/sf_notmerged_down");
  sf_not_down->SetDirectory(0);
  f->Close();
}

void HOTVRScaleFactor::get_sf(double pt, int category) {

  if (category == 3) 
    {
      int bin = sf_merged->FindFixBin(pt);
      if(pt >= 5000.) bin = sf_merged->GetNbinsX();
      
      m_weight *= sf_merged->GetBinContent(bin);
      m_weight_up *= sf_merged_up->GetBinContent(bin);
      m_weight_down *= sf_merged_down->GetBinContent(bin);
    }
  else if (category == 2)
    {
      int bin = sf_semi->FindFixBin(pt);
      if(pt >= 5000.) bin = sf_semi->GetNbinsX();
      
      m_weight *= sf_semi->GetBinContent(bin);
      m_weight_up *= sf_semi_up->GetBinContent(bin);
      m_weight_down *= sf_semi_down->GetBinContent(bin);
    }
  else
    {
      int bin = sf_not->FindFixBin(pt);
      if(pt >= 5000.) bin = sf_not->GetNbinsX();
      
      m_weight *= sf_not->GetBinContent(bin);
      m_weight_up *= sf_not_up->GetBinContent(bin);
      m_weight_down *= sf_not_down->GetBinContent(bin);
    }  
}

bool HOTVRScaleFactor::process(Event &event) {
  
  vector<GenTop> gentops = event.get(h_tophad);
  if (gentops.size() == 0) return false;
  for (const auto &topjet : *event.topjets)
    {
      int nMatched = 0;
      if (m_id_topjet(topjet, event))
	{
	  for (const auto &subjet : topjet.subjets())
	    {
	      double dRmatch = sqrt(subjet.jetArea()/3.14);
	      for (auto top : gentops)
		{
		  if (deltaR(top.get_b(), subjet.v4()) < dRmatch) ++nMatched;
		  if (deltaR(top.get_q1(), subjet.v4()) < dRmatch) ++nMatched;
		  if (deltaR(top.get_q2(), subjet.v4()) < dRmatch) ++nMatched;
		}
	    }
	}
      get_sf(topjet.pt(), nMatched);
    }
  event.set(h_toptag_weight, m_weight);
  event.set(h_toptag_weight_up, m_weight_up);
  event.set(h_toptag_weight_down, m_weight_down);
  if (m_sys_direction == "nominal") event.weight *= m_weight;
  else if (m_sys_direction == "up") event.weight *= m_weight_up;
  else if (m_sys_direction == "down") event.weight *= m_weight_down;
  return true;   
}

