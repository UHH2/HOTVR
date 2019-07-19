#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace uhh2;
using namespace std;


// shamelessly stolen from JetCorrections.h, deleted L1 corrections
std::unique_ptr<FactorizedJetCorrector> build_corrector(const std::vector<std::string> & filenames){
  std::vector<JetCorrectorParameters> pars;
  for(const auto & filename : filenames){
    pars.emplace_back(locate_file(filename));
  }
  return uhh2::make_unique<FactorizedJetCorrector>(pars);
}

void correct_jet(FactorizedJetCorrector & corrector, Jet & jet, const Event & event, JetCorrectionUncertainty* jec_unc, int jec_unc_direction){
  auto factor_raw = jet.JEC_factor_raw();
  corrector.setJetPt(jet.pt() * factor_raw);
  corrector.setJetEta(jet.eta());
  corrector.setJetE(jet.energy() * factor_raw);
  corrector.setJetA(jet.jetArea());
  corrector.setRho(event.rho);
  auto correctionfactors = corrector.getSubCorrections();
  auto correctionfactor = correctionfactors.back();

  LorentzVector jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);
   
  if(jec_unc_direction!=0){
    if (jec_unc==NULL){
      std::cerr << "JEC variation should be applied, but JEC uncertainty object is NULL! Abort." << std::endl;
      exit(EXIT_FAILURE);
    }
    // ignore jets with very low pt or high eta, avoiding a crash from the JESUncertainty tool
    double pt = jet_v4_corrected.Pt();
    double eta = jet_v4_corrected.Eta();
    if (!(pt<5. || fabs(eta)>5.)) {
      
      jec_unc->setJetEta(eta);
      jec_unc->setJetPt(pt);
	
      double unc = 0.;	  
      if (jec_unc_direction == 1){
	unc = jec_unc->getUncertainty(1);
	correctionfactor *= (1 + fabs(unc));
      } else if (jec_unc_direction == -1){
	unc = jec_unc->getUncertainty(-1);
	correctionfactor *= (1 - fabs(unc));
      }
      jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);
    }
  }
  jet.set_v4(jet_v4_corrected);
  jet.set_JEC_factor_raw(1. / correctionfactor);

}

JetCorrectionUncertainty* corrector_uncertainty(uhh2::Context & ctx, const std::vector<std::string> & filenames, int &direction){
    
  auto dir = ctx.get("jecsmear_direction", "nominal");
  if(dir == "up"){
    direction = 1;
  }
  else if(dir == "down"){
    direction = -1;
  }
  else if(dir != "nominal"){
    // direction = 0 is default
    throw runtime_error("JetCorrector: invalid value jecsmear_direction='" + dir + "' (valid: 'nominal', 'up', 'down')");
  }

  //initialize JetCorrectionUncertainty if shift direction is not "nominal", else return NULL pointer
  if(direction!=0){
    //take name from the L1FastJet correction (0th element of filenames) and replace "L1FastJet" by "Uncertainty" to get the proper name of the uncertainty file
    TString unc_file = locate_file(filenames[0]);
    if (unc_file.Contains("L1FastJet")) {
      unc_file.ReplaceAll("L1FastJet","Uncertainty");
    }
    else if (unc_file.Contains("L2Relative")) {
      unc_file.ReplaceAll("L2Relative","Uncertainty");
    }
    else {
      throw runtime_error("WARNING No JEC Uncertainty File found!");
    }
    JetCorrectionUncertainty* jec_uncertainty = new JetCorrectionUncertainty(unc_file.Data());
    return jec_uncertainty;
  }
  return NULL;
}



HOTVRJetCorrector::HOTVRJetCorrector(Context &ctx, const std::vector<std::string> & filenames) {
  corrector = build_corrector(filenames);
  jec_uncertainty = corrector_uncertainty(ctx, filenames, direction) ;


  std::ifstream in("/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/HOTVR/data/CorrectionFactors.txt");
  std::string line;

  int i = 0;

  while (std::getline(in, line))
    {
      float value;
      int k = 0;
      std::stringstream ss(line);

      while (ss >> value)
        {
	  par[i][k] = value;
	  // cout << "par[" << i << "][" << k << "] = " << par[i][k] << endl;
	  ++k;
        }
      ++i;
    }
}

/** Apply ak4-corrections to all subjets, then combine all subjets
 * back to hotvr jet.
 * To Do: After ak4-corrections apply pt_rec/pt_gen corrections.
 */
bool HOTVRJetCorrector::process(Event &event) {
  assert(event.topjets);
  for(auto & topjet : *event.topjets)
    {
      vector<Jet> new_subjets;
      LorentzVector temp_jet;
      auto subjets = topjet.subjets();
      for (auto & subjet : subjets) 
	{ 
	  subjet.set_JEC_factor_raw(1); // initialize jec factor to 1 (is not set in ntuplewriter)

	  // do ak4 corrections (! USE FILES WITHOUT L1 CORRECTIONS !)
	  correct_jet(*corrector, subjet, event, jec_uncertainty, direction);

	  // // apply additional corrections
	  // double corr_factor = get_factor(subjet.pt(), subjet.eta());
	  // subjet.set_v4(subjet.v4() * corr_factor);
	  // subjet.set_JEC_factor_raw(subjet.JEC_factor_raw() / corr_factor / L1_corr);

	  temp_jet += subjet.v4();
	  new_subjets.push_back(subjet);

	}
      sort_by_pt(new_subjets);
      topjet.set_subjets(move(new_subjets));
      topjet.set_v4(temp_jet);
    }
  return true;
}

double HOTVRJetCorrector::get_factor(double pt, double eta) {
// get eta bin
  int etabin = 0;
  for(unsigned int i = 0; i < 6; i++)
    {
      if(eta > eta_bins[i])
	{
	  etabin = i;
	  break;
	}
    }
  // get factor from function
  // if(pt > 425) pt = 425;
  double factor = par[etabin][0] + exp(par[etabin][2]*pt + par[etabin][1]);
  // double factor = par[etabin][0] + par[etabin][1] * pt + par[etabin][2] * pt * pt;
  return factor;
}

HOTVRJetCorrector::~HOTVRJetCorrector() {}

HOTVRJetLeptonCleaner::HOTVRJetLeptonCleaner(const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection) : 
h_topjetcollection(topjetcollection) {}

bool HOTVRJetLeptonCleaner::process(Event &event) {
  vector<TopJet> &topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.topjets;
  if (event.muons)
    {
      for (const Muon &muo : *event.muons)
	{
	  for (TopJet &topjet : topjets)
	    {
	      vector<Jet> new_subjets;
	      LorentzVector temp_jet;
	      auto subjets = topjet.subjets();
	      for (const Jet &subjet : topjet.subjets())
		{
		  double dRmatch = sqrt(subjet.jetArea()/M_PI);
		  LorentzVector subjet_v4 = subjet.v4();
		  if (deltaR(subjet_v4,muo) < dRmatch)
		      subjet_v4 -= muo.v4(); // remove muon from subjet four-vector
		  temp_jet += subjet.v4();
		  new_subjets.push_back(subjet);
		} // end subjet loop
	      sort_by_pt(new_subjets);
	      topjet.set_subjets(move(new_subjets));
	      topjet.set_v4(temp_jet);
	    } // end topjet loop
	} // end muon loop
    } // endif muons

    if (event.electrons)
    {
      for (const Electron &ele : *event.electrons)
	{
	  for (TopJet &topjet : topjets)
	    {
	      vector<Jet> new_subjets;
	      LorentzVector temp_jet;
	      auto subjets = topjet.subjets();
	      for (const Jet &subjet : topjet.subjets())
		{
		  double dRmatch = sqrt(subjet.jetArea()/M_PI);
		  LorentzVector subjet_v4 = subjet.v4(); 
		  if (deltaR(subjet_v4,ele) < dRmatch)
		      subjet_v4 -= ele.v4(); // remove muon from subjet four-vector
		  temp_jet += subjet.v4();
		  new_subjets.push_back(subjet);
		} // end subjet loop
	      sort_by_pt(new_subjets);
	      topjet.set_subjets(move(new_subjets));
	      topjet.set_v4(temp_jet);
	    } // end topjet loop
	} // end electron loop
    } // endif electrons
    return true;
}
