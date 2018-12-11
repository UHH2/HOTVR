#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include <vector>

/** \brief Search all generator level hadronically decaying top
 * quarks.
 *
 * The found top quarks are written to the event as GenTop. The GenTop
 * class stores the four-vectors of the top quark and it's decay
 * products.
 */

class GenTop {
 public:
  GenTop() {};
  LorentzVector get_top() {return m_top;}
  LorentzVector get_b() {return m_b;}
  LorentzVector get_W() {return m_W;}
  LorentzVector get_q1() {return m_q1;}
  LorentzVector get_q2() {return m_q2;}
  
  void set_top(LorentzVector v4) {m_top = v4;}
  void set_b(LorentzVector v4) {m_b = v4;}
  void set_W(LorentzVector v4) {m_W = v4;}
  void set_q1(LorentzVector v4) {m_q1 = v4;}
  void set_q2(LorentzVector v4) {m_q2 = v4;}

 private:
  LorentzVector m_top, m_b, m_W, m_q1, m_q2;
};

class HadronicTop: public uhh2::AnalysisModule {
 public:
  explicit HadronicTop(uhh2::Context &ctx,
		       const std::string & handle_name = "HadronicTop");
  virtual bool process(uhh2::Event &event) override;

  virtual ~HadronicTop();

 private:
  uhh2::Event::Handle<std::vector<GenTop> > h_tophad;
};
