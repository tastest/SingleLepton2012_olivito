#ifndef STACKHISTS
#define STACKHISTS

#include "TFile.h"
#include <vector>
#include <string>
#include "plotUtilities.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"

struct histStyle {
  int nbins;
  float xmin;
  float xmax;
  std::string xtitle;
  bool log;
};

TFile *data_;
TFile *tt_;
TFile *tt0l_;
TFile *tt1l_;
TFile *tt2l_;
TFile *wjets_;
TFile *zjets_;
TFile *vv_;
TFile *t_;
TFile *ttV_;
TFile *vvv_;
TFile *rare_;

TFile *TChiwh_150_1_;
TFile *TChiwh_200_1_;
TFile *TChiwh_250_1_;
TFile *TChiwh_300_1_;
TFile *TChiwh_350_1_;

vector<TFile*>  mcfiles_;
vector<char*>   mclabels_;

histStyle getHistStyle(char* hist);
void initialize( char* path, bool doData = true, bool doSig = false );
//TH1F* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData = true, int rebin = 1, bool normalize = false );
TGraphErrors* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData = true, 
                             int rebin = 1, bool normalize = false, float mcnorm = -1. );
//TCanvas* stackHistAuto( char* path, char* hist, char* flavor, char* dir, bool doData = true, 
//                             int rebin = 1, bool normalize = false, float mcnorm = -1. );
void printYieldsDir( char* path, char* dir, bool doData, bool latex );

#endif
