#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace uhh2;
using namespace std;

HOTVRTopTag::HOTVRTopTag(double fpt_upper, double mjet_lower, double mjet_upper, double mmin_lower):
  m_fpt_upper(fpt_upper), m_mjet_lower(mjet_lower), m_mjet_upper(mjet_upper), m_mmin_lower(mmin_lower) {}

bool HOTVRTopTag::operator()(const TopJet &topjet, const Event &event) const {
  vector<Jet> subjets = topjet.subjets();
  if(subjets.size() < 3) return false;

  double fpt = subjets.at(0).v4().pt() / topjet.v4().pt();
  if (fpt > m_fpt_upper) return false;

  double mjet = topjet.v4().M();
  if(mjet < m_mjet_lower) return false;
  if(mjet > m_mjet_upper) return false;

  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  double mmin = min(min(m12, m13), m23);

  if(mmin < m_mmin_lower) return false;

  return true;
}

HOTVRInvertedMassCut::HOTVRInvertedMassCut(double fpt_upper, double mjet_lower, double mjet_upper, double mmin_lower):
  m_fpt_upper(fpt_upper), m_mjet_lower(mjet_lower), m_mjet_upper(mjet_upper), m_mmin_lower(mmin_lower) {}

bool HOTVRInvertedMassCut::operator()(const TopJet &topjet, const Event &event) const {
  vector<Jet> subjets = topjet.subjets();
  if(subjets.size() < 3) return false;

  double fpt = subjets.at(0).v4().pt() / topjet.v4().pt();
  if (fpt > m_fpt_upper) return false;

  double mjet = topjet.v4().M();
  // if(mjet < m_mjet_lower) return false;
  // if(mjet > m_mjet_upper) return false;

  if (mjet < m_mjet_upper && mjet > m_mjet_lower) return false;

  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  double mmin = min(min(m12, m13), m23);

  if(mmin < m_mmin_lower) return false;

  return true;
}

DeltaPhiCut::DeltaPhiCut(Context &ctx, double deltaphi_lower, std::string h_name):
  m_deltaphi_lower(deltaphi_lower)
{
  h_primlep = ctx.get_handle<FlavorParticle>(h_name);
}

bool DeltaPhiCut::operator()(const TopJet &topjet, const Event &event) const {
  const Particle &lepton = event.get(h_primlep);
  return ( deltaPhi(lepton.v4(), topjet.v4()) > m_deltaphi_lower );
}

DeltaRCut::DeltaRCut(Context &ctx, double deltaR_lower, std::string h_name):
  m_deltaR_lower(deltaR_lower)
{
  h_primlep = ctx.get_handle<FlavorParticle>(h_name);
}

bool DeltaRCut::operator()(const TopJet &topjet, const Event &event) const {
  const Particle &lepton = event.get(h_primlep);
  return ( deltaR(lepton.v4(), topjet.v4()) > m_deltaR_lower );
}

Tau32Groomed::Tau32Groomed(double tau32_upper):
  m_tau32_upper(tau32_upper) {}

bool Tau32Groomed::operator()(const TopJet &topjet, const Event &event) const {
  double tau32 = topjet.tau3_groomed() / topjet.tau2_groomed();
  return ( tau32 < m_tau32_upper );
}

AntiTau32Groomed::AntiTau32Groomed(double tau32_lower):
  m_tau32_lower(tau32_lower) {}

bool AntiTau32Groomed::operator()(const TopJet &topjet, const Event &event) const {
  double tau32 = topjet.tau3_groomed() / topjet.tau2_groomed();
  return ( tau32 >= m_tau32_lower );
}
