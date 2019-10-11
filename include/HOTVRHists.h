#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Muon.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/common/include/JetIds.h"

#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HOTVRGenIds.h"

#include <vector>

namespace uhh2 {
  /** \brief Hist class for HOTVR jets
   *
   * Fills histogramms for all important variables concerning HOTVR.
   *
   */
  class HOTVRHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<TopJetId> &topjetid = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRHists();
  protected:
    // HOTVR
    TH1F *N_HotvrTopjets, *Pt_HotvrTopjets, *Eta_HotvrTopjets, *M_HotvrTopjets, *A_HotvrTopjets, *Msd_HotvrTopjets;
    TH1F *NSub_HotvrTopjets, *Fpt_HotvrTopjets, *Mpair_HotvrTopjets, *Tau32_HotvrTopjets;
    TH1F *DeltaPhi_MET_HotvrTopjets, *DeltaR_Muo_HotvrTopjets, *DeltaPhi_Muo_HotvrTopjets, *DeltaR_Ele_HotvrTopjets, *DeltaPhi_Ele_HotvrTopjets;
    TH1F *Pt_HotvrTopjets_Sub1, *A_HotvrTopjets_Sub1, *Pt_HotvrTopjets_Sub2, *A_HotvrTopjets_Sub2, *Pt_HotvrTopjets_Sub3, *A_HotvrTopjets_Sub3;
    TH1F *N_HotvrTopjet1, *Pt_HotvrTopjet1,  *Eta_HotvrTopjet1, *M_HotvrTopjet1, *A_HotvrTopjet1;
    TH1F *NSub_HotvrTopjet1, *Fpt_HotvrTopjet1, *Mpair_HotvrTopjet1, *Tau32_HotvrTopjet1;
    TH1F *DeltaPhi_MET_HotvrTopjet1, *DeltaR_Muo_HotvrTopjet1, *DeltaPhi_Muo_HotvrTopjet1, *DeltaR_Ele_HotvrTopjet1, *DeltaPhi_Ele_HotvrTopjet1;
    TH1F *Pt_HotvrTopjet1_Sub1, *A_HotvrTopjet1_Sub1, *Pt_HotvrTopjet1_Sub2, *A_HotvrTopjet1_Sub2, *Pt_HotvrTopjet1_Sub3, *A_HotvrTopjet1_Sub3;
    TH2F *Pt_vs_Eta_HotvrTopjets, *Pt_vs_Eta_HotvrTopjet1, *Pt_vs_Eta_HotvrRebin;
    TH1F *Pt_rebin_HotvrTopjets, *EtaAbs_HotvrTopjets;
    TH1F *NLeadingTopjet;

    boost::optional<TopJetId> m_topjetid;
  };

 /** \brief Hist class for HOTVR genjets
   *
   * Fills histogramms for all important variables concerning HOTVR.
   *
   */
  class HOTVRGenHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRGenHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<GenTopJetId> &gentopjetid = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRGenHists();
  protected:
    // HOTVR
    TH1F *N_HotvrTopjets, *Pt_HotvrTopjets, *Eta_HotvrTopjets, *M_HotvrTopjets, *R_HotvrTopjets;
    TH1F *NSub_HotvrTopjets, *Fpt_HotvrTopjets, *Mpair_HotvrTopjets, *Tau32_HotvrTopjets;
    TH1F *DeltaR_L_HotvrTopjets, *DeltaPhi_L_HotvrTopjets;
    TH1F *Pt_HotvrTopjets_Sub1, *Pt_HotvrTopjets_Sub2, *Pt_HotvrTopjets_Sub3;
    TH1F *N_HotvrTopjet1, *Pt_HotvrTopjet1, *Eta_HotvrTopjet1, *M_HotvrTopjet1, *R_HotvrTopjet1;
    TH1F *NSub_HotvrTopjet1, *Fpt_HotvrTopjet1, *Mpair_HotvrTopjet1, *Tau32_HotvrTopjet1;
    TH1F *DeltaR_L_HotvrTopjet1, *DeltaPhi_L_HotvrTopjet1;
    TH1F *Pt_HotvrTopjet1_Sub1, *Pt_HotvrTopjet1_Sub2, *Pt_HotvrTopjet1_Sub3;
 
    boost::optional<GenTopJetId> m_gentopjetid;
  };

 /** \brief Hist class for pile-up correction studies
   *
   * Fills number of events (n) and number of events surviving p_T >
   * 200 GeV cut before (u) and after (c) L1 correction.
   *
   */
  class HOTVRPileUpHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRPileUpHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~HOTVRPileUpHists();
  protected:

    TH1F *n,*u,*c;

  };

  class HOTVRPerformanceHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRPerformanceHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~HOTVRPerformanceHists();
  protected:
    TH1F *hist_delta_pt_gen_reco;

  };

   
}
