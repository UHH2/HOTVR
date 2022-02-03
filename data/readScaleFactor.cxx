// includes only needed for nicely printing scale factors to terminal:
#include <iostream>
#include <sstream>
#include <iomanip>


//______________________________
// CONSTANTS AND TYPEDEFS

enum class EYear {
  isUL17,
  isUL18,
};

const std::map<EYear, TString> kRootFiles = {
  {EYear::isUL17, "TopTaggingScaleFactors_RunIISummer19UL17_PUPPIv15.root"},
  {EYear::isUL18, "TopTaggingScaleFactors_RunIISummer19UL18_PUPPIv15.root"},
};

enum class ETagger {
  isAK8,
  isAK8_subjetbtag,
  isHOTVR,
};

enum class EWorkingPoint {
  isWP0p38, // for AK8 taggers
  isWP0p47, // %
  isWP0p52, // %
  isWP0p61, // %
  isWP0p69, // %
  isDefault, // for HOTVR only
};

enum class EMergeCategory {
  isFullyMerged,
  isNotMerged,
};

const std::map<ETagger, TString> kTaggerToString = {
  {ETagger::isAK8, "AK8_PUPPI"},
  {ETagger::isAK8_subjetbtag, "AK8_PUPPI_subjetbtag"},
  {ETagger::isHOTVR, "HOTVR_PUPPI"},
};

const std::map<EWorkingPoint, TString> kWorkingPointToString = {
  {EWorkingPoint::isWP0p38, "wp0p38_vt"},
  {EWorkingPoint::isWP0p47, "wp0p47_t"},
  {EWorkingPoint::isWP0p52, "wp0p52_m"},
  {EWorkingPoint::isWP0p61, "wp0p61_l"},
  {EWorkingPoint::isWP0p69, "wp0p69_vl"},
  {EWorkingPoint::isDefault, "default"},
};

const std::map<EWorkingPoint, TString> kWorkingPointToMisTagEffString = {
  {EWorkingPoint::isWP0p38, "mis0p001"},
  {EWorkingPoint::isWP0p47, "mis0p005"},
  {EWorkingPoint::isWP0p52, "mis0p010"},
  {EWorkingPoint::isWP0p61, "mis0p025"},
  {EWorkingPoint::isWP0p69, "mis0p050"},
};

const std::map<EMergeCategory, TString> kMergeCategoryToString = {
  {EMergeCategory::isFullyMerged, "FullyMerged"},
  {EMergeCategory::isNotMerged, "NotMerged"},
};

typedef struct {
  Double_t fCentral = 1; // central value
  Double_t fErrTotUp = 0; // total uncertainty
  Double_t fErrTotDown = 0; // %
  Double_t fErrStatUp = 0; // statistical component
  Double_t fErrStatDown = 0; // %
  Double_t fErrSystUp = 0; // systematic component
  Double_t fErrSystDown = 0; // %
} ScaleFactor_t;


//______________________________
// FUNCTIONS

ScaleFactor_t loadScaleFactor(const Double_t jetPt, const ETagger tagger, const EWorkingPoint wp, const EMergeCategory mergeCat, const EYear year) {
  if((wp == EWorkingPoint::isDefault && tagger != ETagger::isHOTVR) || (tagger == ETagger::isHOTVR && wp != EWorkingPoint::isDefault)) {
    throw std::invalid_argument("Given combination of tagger/WP not valid");
  }
  TFile *rootFile = TFile::Open(kRootFiles.at(year).Data(), "READ");
  const TString graph_basename = kTaggerToString.at(tagger)+"_"+kWorkingPointToString.at(wp)+(tagger == ETagger::isAK8 ? "_"+kWorkingPointToMisTagEffString.at(wp) : "")+"/"+kMergeCategoryToString.at(mergeCat)+"_";
  const TGraphAsymmErrors *graph_tot = (TGraphAsymmErrors*)rootFile->Get((graph_basename+"tot").Data());
  const TGraphAsymmErrors *graph_stat = (TGraphAsymmErrors*)rootFile->Get((graph_basename+"stat").Data());
  const TGraphAsymmErrors *graph_syst = (TGraphAsymmErrors*)rootFile->Get((graph_basename+"syst").Data());

  UInt_t ibin = 0;
  Double_t x;
  Double_t y;
  graph_tot->GetPoint(ibin, x, y);

  ScaleFactor_t sf;

  const Double_t minimum_jetPt = x - graph_tot->GetErrorXlow(ibin);
  if(jetPt < minimum_jetPt) {
    std::cerr << "No scale factor for such low jet pT available. Will set scale factor to 1 with zero uncertainty" << std::endl;
    return sf;
  }

  while(!(jetPt >= x - graph_tot->GetErrorXlow(ibin) && jetPt < x + graph_tot->GetErrorXhigh(ibin))) { // check if jetPt within pT interval of bin with bin number "ibin"
    ++ibin;
    graph_tot->GetPoint(ibin, x, y);
    if(ibin + 1 >= graph_tot->GetN()) break; // use scale factor of last bin if given jetPt exceeds highest bin edge of the graph
  }

  sf = ScaleFactor_t{
    .fCentral = y,
    .fErrTotUp = graph_tot->GetErrorYhigh(ibin),
    .fErrTotDown = graph_tot->GetErrorYlow(ibin),
    .fErrStatUp = graph_stat->GetErrorYhigh(ibin),
    .fErrStatDown = graph_stat->GetErrorYlow(ibin),
    .fErrSystUp = graph_syst->GetErrorYhigh(ibin),
    .fErrSystDown = graph_syst->GetErrorYlow(ibin),
  };

  rootFile->Close();

  return sf;
}

void printScaleFactor(const ScaleFactor_t & sf, const Int_t precision = 3) {
  std::stringstream ss;
  ss << "SF: ";
  ss << std::fixed << std::setprecision(precision) << sf.fCentral;
  ss << " +" << std::fixed << std::setprecision(precision) << sf.fErrTotUp;
  ss << "/";
  ss << "-" << std::fixed << std::setprecision(precision) << sf.fErrTotDown;
  ss << " (tot)";
  ss << " [";
  ss << " +" << std::fixed << std::setprecision(precision) << sf.fErrStatUp;
  ss << "/";
  ss << "-" << std::fixed << std::setprecision(precision) << sf.fErrStatDown;
  ss << " (stat) |";
  ss << " +" << std::fixed << std::setprecision(precision) << sf.fErrSystUp;
  ss << "/";
  ss << "-" << std::fixed << std::setprecision(precision) << sf.fErrSystDown;
  ss << " (syst) ]";
  std::cout << ss.str() << std::endl;
}


//______________________________
// MACRO

void readScaleFactor() {
  // Exemplary configuration:
  const Double_t jetPt = 401.816; // in GeV
  const ETagger tagger = ETagger::isHOTVR;
  const EWorkingPoint wp = EWorkingPoint::isDefault;
  //const EMergeCategory mergeCat = EMergeCategory::isFullyMerged; // in practice, this needs to be determined for each individual jet; see TWiki how to do it
  const EMergeCategory mergeCat = EMergeCategory::isNotMerged; 
  const EYear year = EYear::isUL18;

  const ScaleFactor_t sf = loadScaleFactor(jetPt, tagger, wp, mergeCat, year);
  printScaleFactor(sf);
}
