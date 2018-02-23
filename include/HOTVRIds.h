#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"


/** \brief HOTVR top tagger
 * 
 * Cuts on the pt fraction of the leading subjet, the minimum pairwise
 * subjet mass, and on the top jet mass.  A cut on n_subjets >= 3 is
 * always applied, as the calculation of the minumum pairwise mass
 * requires three subjets to be present.
 * 
 * The default values of the constructor correspond to the default values used in HOTVR.
 */
class HOTVRTopTag {
 public:
  explicit HOTVRTopTag(double fpt_upper = 0.8, double mjet_lower = 140., double mjet_upper = 220., double mmin_lower = 50.);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_fpt_upper;
  double m_mjet_lower;
  double m_mjet_upper;
  double m_mmin_lower;

};

/** \brief minimum delta phi between HOTVR jet and primary lepton
 *
 * Selects HOTVR jets with deltaPhi > deltaphi_lower, where deltaPhi
 * is the distance between the HOTVR jet and the primary lepton in
 * phi.
 *
 */
class DeltaPhiCut {
 public: 
  explicit DeltaPhiCut(uhh2::Context &ctx, double deltaphi_lower = M_PI/2, std::string h_name = "PrimaryLepton");
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_deltaphi_lower;
  uhh2::Event::Handle<FlavorParticle> h_primlep;

};

/** \breif Tau32 Cut
 *
 * Cut on N-subjettiness ration tau3/tau2. tau32_groomed is used, as
 * HOTVR jets are assumed to be groomed after clustering.
 *
 */
class Tau32Groomed {
 public: 
  explicit Tau32Groomed(double tau32_upper);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_tau32_upper;
};
  

/** \breif AntiTau32 Cut
 *
 * Cut on N-subjettiness ration tau3/tau2. tau32_groomed is used, as
 * HOTVR jets are assumed to be groomed after clustering.
 *
 */
class AntiTau32Groomed {
 public: 
  explicit AntiTau32Groomed(double tau32_lower);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_tau32_lower;
};
  
