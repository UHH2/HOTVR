#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/Utils.h"

using namespace uhh2;
using namespace std;

HOTVRJetCorrectionModule::HOTVRJetCorrectionModule(Context & ctx) {

  is_mc = ctx.get("dataset_type") == "MC";

  string jec_tag_2016 = "Summer16_07Aug2017";
  string jec_ver_2016 = "11";

  string jec_tag_2017 = "Fall17_17Nov2017";
  string jec_ver_2017 = "32";

  string jec_tag_2018 = "Autumn18";
  string jec_ver_2018 = "19";

  string jec_tag_UL16preVFP = "Summer19UL16APV"; 
  string jec_ver_UL16preVFP = "7"; 

  string jec_tag_UL16postVFP = "Summer19UL16"; 
  string jec_ver_UL16postVFP = "7"; 

  string jec_tag_UL17 = "Summer19UL17";
  string jec_ver_UL17 = "5";

  string jec_tag_UL18 = "Summer19UL18";
  string jec_ver_UL18 = "5";

  string jec_jet_coll = "AK4PFPuppi";

  h_subj = ctx.get_handle<vector<Jet> >("hotvr_subjets");
  h_subjmap = ctx.get_handle<vector<pair<int, int> > >("hotvr_subjmap");
  h_gensubj = ctx.get_handle<vector<GenJet> >("hotvr_gensubjets");

  // setup GenericJetCorrector for different years / runs
  if (is_mc)
    {
      jet_corrector_MC.reset(new YearSwitcher(ctx));
      jet_corrector_MC->setup2016(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setup2017(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setup2018(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setupUL16preVFP(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setupUL16postVFP(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setupUL17(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_jet_coll),"hotvr_subjets"));
      jet_corrector_MC->setupUL18(std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_jet_coll),"hotvr_subjets"));
    }
  else
    {
      jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
        jec_switcher_16->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
      for (const auto & runItr : runPeriods2017) {
        jec_switcher_17->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
      for (const auto & runItr : runPeriods2018) {
        jec_switcher_18->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_UL16preVFP.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriodsUL16preVFP) { 
        jec_switcher_UL16preVFP->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_UL16postVFP.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriodsUL16postVFP) { 
        jec_switcher_UL16postVFP->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_UL17.reset(new RunSwitcher(ctx, "2017"));
      for (const auto & runItr : runPeriods2017) {
        jec_switcher_UL17->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jec_switcher_UL18.reset(new RunSwitcher(ctx, "2018"));
      for (const auto & runItr : runPeriods2018) {
        jec_switcher_UL18->setupRun(runItr, std::make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_jet_coll, runItr),"hotvr_subjets"));
      }

      jet_corrector_data.reset(new YearSwitcher(ctx));
      jet_corrector_data->setup2016(jec_switcher_16);
      jet_corrector_data->setup2017(jec_switcher_17);
      jet_corrector_data->setup2018(jec_switcher_18);
      jet_corrector_data->setupUL16preVFP(jec_switcher_UL16preVFP);
      jet_corrector_data->setupUL16postVFP(jec_switcher_UL16postVFP);
      jet_corrector_data->setupUL17(jec_switcher_UL17);
      jet_corrector_data->setupUL18(jec_switcher_UL18);
    }

  // finding JER file as in JetCorrections.cxx JetResolutionSmearer constructor
  std::string jetCollection = "AK4PFPuppi";
  const Year & year = extract_year(ctx);
  std::string jer_tag = "";
  if (year == Year::is2016v2 || year == Year::is2016v3) {
    jer_tag = "Summer16_25nsV1";
  } else if (year == Year::is2017v1 || year == Year::is2017v2) {
    jer_tag = "Fall17_V3";
  } else if (year == Year::is2018) {
    jer_tag = "Autumn18_V7";
  } else if (year == Year::isUL16preVFP) {
    jer_tag = "Summer20UL16APV_JRV3";
  } else if (year == Year::isUL16postVFP) {
    jer_tag = "Summer20UL16_JRV3";
  } else if (year == Year::isUL17) {
    jer_tag = "Summer19UL17_JRV2";
  } else if (year == Year::isUL18) {
    jer_tag = "Summer19UL18_JRV2";
  } else {
    throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
  }

  jer_module.reset(new GenericJetResolutionSmearer(ctx, "hotvr_subjets", "hotvr_gensubjets", JERFiles::JERPathStringMC(jer_tag, jetCollection, "SF"), JERFiles::JERPathStringMC(jer_tag, jetCollection, "PtResolution")));

}

void HOTVRJetCorrectionModule::set_subjet_handles(Event &event) {
  vector<GenJet> gensubjs;
  if (!event.isRealData)
    {
      for (const GenTopJet &genjet : *event.gentopjets)
	{
	  for (GenJet gensubj : genjet.subjets())
	    {
	      gensubjs.push_back(gensubj);
	    }
	}
    }
  event.set(h_gensubj, gensubjs);

  vector<Jet> subjs;
  vector<TopJet> *topjets = event.topjets;
  int j = 0;
  vector<pair<int, int> > subj_map;
  for (unsigned int i = 0; i < topjets->size(); ++i)
    {
      for (Jet subj : topjets->at(i).subjets())
	{
	  subjs.push_back(subj);
	  subj_map.push_back(pair<int, int>({j,i}));
	  ++j;
	}
    }
  event.set(h_subj, subjs);
  event.set(h_subjmap, subj_map);
}

void HOTVRJetCorrectionModule::rebuild_jets(Event &event) {
  vector<Jet> subjets = event.get(h_subj);
  vector<pair<int, int> > subj_map = event.get(h_subjmap);
  int i = 0;
  for (auto &topjet : *event.topjets)
    {
      vector<Jet> new_subjets;
      LorentzVector v4;
      for (auto & p : subj_map)
	{
	  if (p.second == i)
	    {
	      new_subjets.push_back(subjets.at(p.first));
	      v4 += subjets.at(p.first).v4();
	    }
	}
      topjet.set_subjets(new_subjets);
      double jec_factor_raw = topjet.pt() / v4.pt(); // calculate JEC_factor raw to use raw pt later
      topjet.set_JEC_factor_raw(jec_factor_raw);
      topjet.set_v4(v4);
      ++i;
    }
}

bool HOTVRJetCorrectionModule::process(Event &event) {

  // fill all (gen)subjets into single vector and create map to identify subjets with their mother jet
  set_subjet_handles(event);

  // run JEC and JER
  if (is_mc)
    {
       jet_corrector_MC->process(event);
       jer_module->process(event);
    }
  else
    {
      jet_corrector_data->process(event);
    }

  // rebuild HOTVR jets from subjets using the mapping
  rebuild_jets(event);

  return true;
}
