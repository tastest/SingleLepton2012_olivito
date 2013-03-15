#include <vector>
#include <iostream>
#include <string>

#include "TFile.h"
#include "plotUtilities.C"
#include "stackHists.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

bool alreadyInitialized_ = false;
std::string currentPath_;

void initialize( const char* path, bool doData, int doSig ) {

  if (currentPath_.size()) {
    if (currentPath_ == std::string(path)) {
      std::cout << "Using current dir" << std::endl;
      return;
    }

    std::cout << "Re-intializing inputs" << std::endl;

    if (doData && data_) data_->Close();
    wjets_->Close();
    //    wbb_->Close();
    zjets_->Close();
    tt_->Close();
    //    tt0l_->Close();
    tt1l_->Close();
    tt2l_->Close();
    vv_->Close();
    t_->Close();
    ttV_->Close();
    vvv_->Close();
    rare_->Close();
    whbb_->Close();

    if (doSig > 0) {
      TChiwh_250_1_->Close();
    }
    if (doSig > 1) {
      TChiwh_200_1_->Close();
    }
    if (doSig > 2) {
      TChiwh_300_1_->Close();
    }
    if (doSig > 3) {
      TChiwh_150_1_->Close();
      TChiwh_250_25_->Close();
      TChiwh_250_50_->Close();
      TChiwh_350_1_->Close();
    }

    mcfiles_.clear();
    mclabels_.clear();

    // should delete old pointers here..
  }

  std::cout << "Reading files from path: " << path << std::endl;

  if (doData) data_ = new TFile(Form("%s/data_histos.root",path));
  //  wjets_ = new TFile(Form("%s/wjets_histos.root",path));
  wjets_ = new TFile(Form("%s/wjets_comb_histos.root",path));
  //  wbb_ = new TFile(Form("%s/wbb_histos.root",path));
  zjets_ = new TFile(Form("%s/zjets_histos.root",path));
  tt_ = new TFile(Form("%s/ttbar_histos.root",path));
  //  tt0l_ = new TFile(Form("%s/ttbar_0l_histos.root",path));
  // tt1l_ = new TFile(Form("%s/ttbar_1l_histos.root",path));
  // tt2l_ = new TFile(Form("%s/ttbar_2l_histos.root",path));
  tt1l_ = new TFile(Form("%s/ttbar_mg_1l_histos.root",path));
  tt2l_ = new TFile(Form("%s/ttbar_mg_2l_histos.root",path));
  vv_ = new TFile(Form("%s/VV_histos.root",path));
  t_ = new TFile(Form("%s/single_top_histos.root",path));
  ttV_ = new TFile(Form("%s/ttV_histos.root",path));
  vvv_ = new TFile(Form("%s/VVV_histos.root",path));
  rare_ = new TFile(Form("%s/rare_histos.root",path));
  whbb_ = new TFile(Form("%s/whbb_histos.root",path));

  if (doSig > 0) {
    TChiwh_250_1_ = new TFile(Form("%s/TChiwh_250_1_histos.root",path));
  }
  if (doSig > 1) {
    TChiwh_200_1_ = new TFile(Form("%s/TChiwh_200_1_histos.root",path));
  }
  if (doSig > 2) {
    TChiwh_300_1_ = new TFile(Form("%s/TChiwh_300_1_histos.root",path));
  }
  if (doSig > 3) {
    TChiwh_150_1_ = new TFile(Form("%s/TChiwh_150_1_histos.root",path));
    TChiwh_250_25_ = new TFile(Form("%s/TChiwh_250_25_histos.root",path));
    TChiwh_250_50_ = new TFile(Form("%s/TChiwh_250_50_histos.root",path));
    TChiwh_350_1_ = new TFile(Form("%s/TChiwh_350_1_histos.root",path));
  }

  //  mcfiles_.push_back(tt_);    mclabels_.push_back("ttbar");
  mcfiles_.push_back(tt2l_);    mclabels_.push_back("ttbar 2l");
  mcfiles_.push_back(tt1l_);    mclabels_.push_back("ttbar 1l");
  //  mcfiles_.push_back(tt0l_);    mclabels_.push_back("ttbar 0l");
  mcfiles_.push_back(wjets_); mclabels_.push_back("wjets");
  //  mcfiles_.push_back(wbb_); mclabels_.push_back("wbb");
  mcfiles_.push_back(t_);     mclabels_.push_back("single_top");
  //  mcfiles_.push_back(zjets_); mclabels_.push_back("zjets");
  //  mcfiles_.push_back(vv_);    mclabels_.push_back("VV");
  //  mcfiles_.push_back(ttV_);   mclabels_.push_back("ttV");
  //  mcfiles_.push_back(vvv_);   mclabels_.push_back("VVV");
  mcfiles_.push_back(rare_);     mclabels_.push_back("rare");
  mcfiles_.push_back(whbb_);     mclabels_.push_back("whbb");

  if (doSig == 1) {
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  }
  else if (doSig == 2) {
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  }
  else if (doSig == 3) {
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
    mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("TChiwh_300_1");
  }
  else if (doSig > 3) {
    mcfiles_.push_back(TChiwh_150_1_);     mclabels_.push_back("TChiwh_150_1");
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
    mcfiles_.push_back(TChiwh_250_25_);     mclabels_.push_back("TChiwh_250_25");
    mcfiles_.push_back(TChiwh_250_50_);     mclabels_.push_back("TChiwh_250_50");
    mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("TChiwh_300_1");
    mcfiles_.push_back(TChiwh_350_1_);     mclabels_.push_back("TChiwh_350_1");
  }

  currentPath_ = std::string(path);
  return;
}

void stackHist( const char* path, const char* hist, const char* flavor, const char* dir, int nbins, float xmin, float xmax, const char* xtitle ) {
  initialize(path);
  TCanvas* c = new TCanvas();
  compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir , 
		 nbins ,  xmin , xmax , xtitle ,
		 true , true , true , false );

  return;
}

//TH1F* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize ) {
// TGraphErrors* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize, float mcnorm ) {
TCanvas* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize, float mcnorm, int doSig ) {
  initialize(path, doData, doSig);
  histStyle style = getHistStyle(hist);
  TCanvas* c = new TCanvas();
  compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
  			doData , doData , true , style.log, normalize , false , mcnorm );
  // compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  // 		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
  // 			doData , doData , true , style.log, normalize , false , mcnorm );

  return c;
}

void printYieldsDir( const char* path, const char* dir, bool doData, int latex, int doSig ) {
  initialize(path, doData, doSig);
  printYields(mcfiles_, mclabels_, data_, dir, doData, latex);
}

void saveAllHists( const char* path, const char* dir, const char* flavor, const char* outpath, bool doData, int rebin, bool normalize, float mcnorm ) {

  std::vector<std::string> histnames;
  histnames.push_back("h_bbmass");
  histnames.push_back("h_pfmet");
  histnames.push_back("h_lep1mt");
  histnames.push_back("h_mt2bl");
  histnames.push_back("h_mt2b");
  histnames.push_back("h_mt2w");
  histnames.push_back("h_bbpt");
  histnames.push_back("h_wpt");
  histnames.push_back("h_lep1pt");
  histnames.push_back("h_lep1metdphi");
  histnames.push_back("h_njets");
  histnames.push_back("h_njetsalleta");
  histnames.push_back("h_nbjets");
  histnames.push_back("h_jet1pt");
  histnames.push_back("h_jet2pt");
  histnames.push_back("h_bjet1pt");
  histnames.push_back("h_bjet2pt");
  histnames.push_back("h_jet1csv");
  histnames.push_back("h_jet2csv");
  //  if (TString(dir).Contains("cr3") || TString(dir).Contains("cr6")) {
  if (TString(dir).Contains("cr3")) {
    //    histnames.clear();
    histnames.push_back("h_pseudomet_lep");
    histnames.push_back("h_pseudomt_lep");
    histnames.push_back("h_pseudomt2bl");
    histnames.push_back("h_leppt");
    histnames.push_back("h_dphi_pseudomet_lep");
    histnames.push_back("h_dildphi");
  }

  for (unsigned int i = 0; i < histnames.size(); ++i) {
    std::string label_flavor = std::string(flavor);
    if (std::string(flavor).length() > 0) histnames[i] += std::string("_") + std::string(flavor);
    else if (TString(dir).Contains("cr3") || TString(dir).Contains("cr6")) label_flavor = "all";
    TCanvas* c = stackHistAuto(path, histnames[i].c_str(), label_flavor.c_str(), dir, doData, rebin, normalize, mcnorm);
    c->SaveAs(Form("%s/%s_%s.eps",outpath,dir,TString(histnames[i]).ReplaceAll("h_","").Data()));
  }

}

histStyle getHistStyle( const char* hist ) {

  histStyle style;
  std::string histString = std::string(hist);

  if (matchHistName(histString,"h_ht30")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40up")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{up} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40dn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{dn} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilmass")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 200.;
    style.xtitle = "m(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    // style.nbins = 160;
    // style.xmin = 0.;
    // style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_ee")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_mm")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dileta")) {
    style.nbins = 30;
    style.xmin = -6.0;
    style.xmax = 6.0;
    style.xtitle = "#eta(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity_lowdilpt")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity_highdilpt")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dildeltar")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dildeltaphi")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep2pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep3pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep1mt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "m_{T}(lep1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep2mt")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep3mt")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep1phi")) {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep1)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep2phi")) {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep2)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep1,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_mt2b")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{T2}^{b} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_mt2bl")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{T2}^{bl} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_mt2w")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{T2}^{W} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet1eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1eta_lowdilpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity_lowdilpt")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1eta_highdilpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity_highdilpt")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(jet1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet2pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(jet2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet1pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(bjet1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet2pt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(bjet2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet1eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet2eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_bbpt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(b#bar{b}) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bbmass")) {
    style.nbins = 50;
    style.xmin = 0.;
    style.xmax = 1000.;
    style.xtitle = "M(b#bar{b}) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bbdr")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdphi")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_wpt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(W) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pfmet")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pf MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_njets")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, |#eta| < 2.4)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsalleta")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, |#eta| < 4.7)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsup")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES up";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsdn")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES down";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsupdiff")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets) up - N(jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nbjets")) {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(b jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjets")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjetsup")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES up";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjetsdn")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES down";
    style.log = true;
  } else if (matchHistName(histString,"h_njets40")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets)";
    style.log = true; 
  } else if (matchHistName(histString,"h_nbcsvm")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Medium)";
    style.log = true;
  } else if (matchHistName(histString,"h_nbcsvl")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Loose)";
    style.log = true;
  } else if (matchHistName(histString,"h_nmistags")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(mistags, csv Medium)";
    style.log = true;
  } else if (matchHistName(histString,"h_ngenb")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(reco jets matched to gen b quarks)";
    style.log = true;
  } else if (matchHistName(histString,"h_nlep")) {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(lep)";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt_ee")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt_mm")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptup")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptdn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptup10")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up x 10 [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptdn10")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down x 10 [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ptpar")) {
    style.nbins = 100;
    style.xmin = -100.;
    style.xmax = 400.;
    style.xtitle = "jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ptperp")) {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "jet vec sum p_{T}^{#perp} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumpt")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumptup")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumptdn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_minus_ptpar_b2b")) {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "dilep p_{T} - jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_deltaphi_dil_sum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(dilep, jet vec sum)";
    style.log = false;
  } else if (matchHistName(histString,"h_deltaphi_z_sum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(gen Z, jet system)";
    style.log = false;
  } else if (matchHistName(histString,"h_deltaphi_ttbar_nonbsum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(true ttbar, non-btagged jet vec sum)";
    style.log = false;
  } else if (matchHistName(histString,"h_nvtx")) {
    style.nbins = 40;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "N(vtx)";
    style.log = false;
  } else if (matchHistName(histString,"h_vtx")) {
    style.nbins = 40;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "N(vtx)";
    style.log = false;
  } else if (matchHistName(histString,"h_pseudomet_lep")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pseudomt_lep")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo m_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pseudomt2bl")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo M_{T2}^{bl} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_leppt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dildphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep1,lep2)";
    style.log = false;
  } else if (matchHistName(histString,"h_dphi_pseudomet_lep")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep,pseudo MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1csv")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 1.0;
    style.xtitle = "CSV(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet2csv")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 1.0;
    style.xtitle = "CSV(jet2)";
    style.log = false;

  } else {
    std::cout << "ERROR: hist name " << histString << " not recognized" << std::endl;
  }

  return style;
}

bool matchHistName(const TString& histname, const TString& matchname) {

  // first check for a simple match
  if (histname.EqualTo(matchname)) return true;

  // otherwise check for a name match with a channel suffix
  if (histname.Contains(matchname)) {
    if (histname.EndsWith("_e")) return true;
    if (histname.EndsWith("_m")) return true;
    if (histname.EndsWith("_ee")) return true;
    if (histname.EndsWith("_mm")) return true;
    if (histname.EndsWith("_em")) return true;
  }

  return false;
}

bool matchHistName(const TString& histname, const char* matchname) {
  return matchHistName(histname,TString(matchname));
}
