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

  string jec_jet_coll = "AK4PFPuppi";

  // setup HOTVRJetCorrector for different years / runs
  jlc_module.reset(new HOTVRJetLeptonCleaner());

  if (is_mc)
    {
      jet_corrector_MC.reset(new YearSwitcher(ctx));
      jet_corrector_MC->setup2016(std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll)));
      jet_corrector_MC->setup2017(std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll)));
      jet_corrector_MC->setup2018(std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll)));
    }
  else
    {
      jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
	jec_switcher_16->setupRun(runItr, std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr)));
      }

      jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
      for (const auto & runItr : runPeriods2017) {
	jec_switcher_17->setupRun(runItr, std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr)));
      }

      jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
      for (const auto & runItr : runPeriods2018) {
	jec_switcher_18->setupRun(runItr, std::make_shared<HOTVRJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr)));
      }

      jet_corrector_data.reset(new YearSwitcher(ctx));
      jet_corrector_data->setup2016(jec_switcher_16);
      jet_corrector_data->setup2017(jec_switcher_17);
      jet_corrector_data->setup2018(jec_switcher_18);
    }
}

bool HOTVRJetCorrectionModule::process(Event &event) {

  jlc_module->process(event);

  if (is_mc)
    {
      jet_corrector_MC->process(event);
    }
  else
    {
      jet_corrector_data->process(event);
    }
  return true;
}
