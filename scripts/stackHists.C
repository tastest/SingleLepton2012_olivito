#include <vector>
#include <iostream>
#include <string>

#include "TFile.h"
#include "plotUtilities.C"
#include "stackHists.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"

bool alreadyInitialized_ = false;
std::string currentPath_;

void initialize( char* path, bool doData ) {

  if (currentPath_.size()) {
    if (currentPath_ == std::string(path)) {
      std::cout << "Using current dir" << std::endl;
      return;
    }

    std::cout << "Re-intializing inputs" << std::endl;

    if (doData && data_) data_->Close();
    wjets_->Close();
    zjets_->Close();
    tt_->Close();
    vv_->Close();
    t_->Close();
    ttV_->Close();
    vvv_->Close();
    mcfiles_.clear();
    mclabels_.clear();

    // should delete old pointers here..
  }

  std::cout << "Reading files from path: " << path << std::endl;

  if (doData) data_ = new TFile(Form("%s/data_histos.root",path));
  wjets_ = new TFile(Form("%s/wjets_histos.root",path));
  zjets_ = new TFile(Form("%s/zjets_histos.root",path));
  tt_ = new TFile(Form("%s/ttbar_histos.root",path));
  vv_ = new TFile(Form("%s/VV_histos.root",path));
  t_ = new TFile(Form("%s/single_top_histos.root",path));
  ttV_ = new TFile(Form("%s/ttV_histos.root",path));
  vvv_ = new TFile(Form("%s/VVV_histos.root",path));

  mcfiles_.push_back(wjets_); mclabels_.push_back("wjets");
  mcfiles_.push_back(zjets_); mclabels_.push_back("zjets");
  mcfiles_.push_back(tt_);    mclabels_.push_back("ttbar");
  mcfiles_.push_back(vv_);    mclabels_.push_back("VV");
  mcfiles_.push_back(t_);     mclabels_.push_back("single_top");
  mcfiles_.push_back(ttV_);   mclabels_.push_back("ttV");
  mcfiles_.push_back(vvv_);   mclabels_.push_back("VVV");

  currentPath_ = std::string(path);
  return;
}

void stackHist( char* path, char* hist, char* flavor, char* dir, int nbins, float xmin, float xmax, char* xtitle ) {
  initialize(path);
  TCanvas* c = new TCanvas();
  compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir , 
		 nbins ,  xmin , xmax , xtitle ,
		 true , true , true , false );

  return;
}

//TH1F* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData, int rebin, bool normalize ) {
TGraphErrors* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData, int rebin, bool normalize, float mcnorm ) {
//TCanvas* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData, int rebin, bool normalize, float mcnorm ) {
  initialize(path, doData);
  histStyle style = getHistStyle(hist);
  TCanvas* c = new TCanvas();
  return compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
  			doData , doData , true , style.log, normalize , false , mcnorm );
  // compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  // 		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
  // 			doData , doData , true , style.log, normalize , false , mcnorm );

  //  return c;
}

histStyle getHistStyle( char* hist ) {

  histStyle style;
  std::string histString = std::string(hist);

  if (histString == "h_ht30") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (histString == "h_ht40") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (histString == "h_ht40up") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{up} [GeV]";
    style.log = true;
  } else if (histString == "h_ht40dn") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{dn} [GeV]";
    style.log = true;
  } else if (histString == "h_dilmass") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 200.;
    style.xtitle = "m(ll) [GeV]";
    style.log = true;
  } else if (histString == "h_dilpt") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    // style.nbins = 160;
    // style.xmin = 0.;
    // style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (histString == "h_dilpt_ee") {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (histString == "h_dilpt_mm") {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (histString == "h_dileta") {
    style.nbins = 30;
    style.xmin = -6.0;
    style.xmax = 6.0;
    style.xtitle = "#eta(ll)";
    style.log = false;
  } else if (histString == "h_dilrapidity") {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (histString == "h_dilrapidity_lowdilpt") {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (histString == "h_dilrapidity_highdilpt") {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (histString == "h_dildeltar") {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(ll)";
    style.log = false;
  } else if (histString == "h_dildeltaphi") {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(ll)";
    style.log = false;
  } else if (histString == "h_lep1pt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep1) [GeV]";
    style.log = true;
  } else if (histString == "h_lep2pt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep2) [GeV]";
    style.log = true;
  } else if (histString == "h_lep3pt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep3) [GeV]";
    style.log = true;
  } else if (histString == "h_lep1mt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep1) [GeV]";
    style.log = true;
  } else if (histString == "h_lep2mt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep2) [GeV]";
    style.log = true;
  } else if (histString == "h_lep3mt") {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep3) [GeV]";
    style.log = true;
  } else if (histString == "h_lep1phi") {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep1)";
    style.log = false;
  } else if (histString == "h_lep2phi") {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep2)";
    style.log = false;
  } else if (histString == "h_jet1eta") {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (histString == "h_jet1rapidity") {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (histString == "h_jet1eta_lowdilpt") {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (histString == "h_jet1rapidity_lowdilpt") {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (histString == "h_jet1eta_highdilpt") {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (histString == "h_jet1rapidity_highdilpt") {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (histString == "h_pfmet") {
    style.nbins = 400;
    style.xmin = 0.;
    style.xmax = 200.;
    style.xtitle = "pf MET [GeV]";
    style.log = true;
  } else if (histString == "h_njets") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets)";
    style.log = true;
  } else if (histString == "h_njetsup") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES up";
    style.log = true;
  } else if (histString == "h_njetsdn") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES down";
    style.log = true;
  } else if (histString == "h_njetsupdiff") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets) up - N(jets)";
    style.log = true;
  } else if (histString == "h_nnonbjets") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets)";
    style.log = true;
  } else if (histString == "h_nnonbjetsup") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES up";
    style.log = true;
  } else if (histString == "h_nnonbjetsdn") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES down";
    style.log = true;
  } else if (histString == "h_njets40") {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets)";
    style.log = true; 
  } else if (histString == "h_nbcsvm") {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Medium)";
    style.log = true;
  } else if (histString == "h_nbcsvl") {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Loose)";
    style.log = true;
  } else if (histString == "h_nmistags") {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(mistags, csv Medium)";
    style.log = true;
  } else if (histString == "h_ngenb") {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(reco jets matched to gen b quarks)";
    style.log = true;
  } else if (histString == "h_nlep") {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(lep)";
    style.log = true;
  } else if (histString == "h_sumpt") {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (histString == "h_sumpt_ee") {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (histString == "h_sumpt_mm") {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (histString == "h_sumptup") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (histString == "h_sumptdn") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (histString == "h_sumptup10") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up x 10 [GeV]";
    style.log = true;
  } else if (histString == "h_sumptdn10") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down x 10 [GeV]";
    style.log = true;
  } else if (histString == "h_ptpar") {
    style.nbins = 100;
    style.xmin = -100.;
    style.xmax = 400.;
    style.xtitle = "jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (histString == "h_ptperp") {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "jet vec sum p_{T}^{#perp} [GeV]";
    style.log = true;
  } else if (histString == "h_nonbsumpt") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (histString == "h_nonbsumptup") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (histString == "h_nonbsumptdn") {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (histString == "h_dilpt_minus_ptpar_b2b") {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "dilep p_{T} - jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (histString == "h_deltaphi_dil_sum") {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(dilep, jet vec sum)";
    style.log = false;
  } else if (histString == "h_deltaphi_z_sum") {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(gen Z, jet system)";
    style.log = false;
  } else if (histString == "h_deltaphi_ttbar_nonbsum") {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(true ttbar, non-btagged jet vec sum)";
    style.log = false;
  } else if (histString == "h_nvtx") {
    style.nbins = 40;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "N(vtx)";
    style.log = false;

  } else {
    std::cout << "ERROR: hist name " << histString << " not recognized" << std::endl;
  }

  return style;
}
