#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenTopJet.h"
#include "UHH2/common/include/ObjectIdUtils.h"


/*
 * Implementation of GenTopJetId
 */
typedef std::function<bool (const GenTopJet &, const uhh2::Event &)> GenTopJetId;

/** \brief HOTVR top tagger
 * 
 * Cuts on the pt fraction of the leading subjet, the minimum pairwise
 * subjet mass, and on the top jet mass.  A cut on n_subjets >= 3 is
 * always applied, as the calculation of the minumum pairwise mass
 * requires three subjets to be present.
 * 
 * The default values of the constructor correspond to the default values used in HOTVR.
 */
class GenHOTVRTopTag {
 public:
  explicit GenHOTVRTopTag(double fpt_upper = 0.8, double mjet_lower = 140., double mjet_upper = 220., double mmin_lower = 50.);
  bool operator()(const GenTopJet &gentopjet, const uhh2::Event &event) const;

 private:
  double m_fpt_upper;
  double m_mjet_lower;
  double m_mjet_upper;
  double m_mmin_lower;

};

