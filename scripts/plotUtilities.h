#ifndef PLOTUTILITIES
#define PLOTUTILITIES

#include "TFile.h"
#include "TChain.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <string>

struct dataMCHists{
  TPad* comppad;
  TH1F* datahist;
  TH1F* mctothist;
  TH1F* ratio;
};

//TH1F* compareDataMC( std::vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
TGraphErrors* compareDataMC( std::vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
		    int nbins ,  float xmin , float xmax , const char* xtitle , 
		    bool overlayData = true , bool residual = false, bool drawLegend = true , 
			     bool log = false , bool normalize = false , bool fit = false, float mcnorm = -1. );

void initSymbols(bool);
void  printLine(bool);
void printYields( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* dir , bool doData, bool latex = false );
void fillYieldHist( TH1F* hin, TH1F* hyield, int bin );
//void deleteHistos();
TLegend *getLegend( vector<char*> labels , bool overlayData, 
		    float x1 = 0.7, float y1 = 0.45 , float x2 = 0.87 , float y2 = 0.94 );

TGraphAsymmErrors* makeBand(TH1F* centhist, TH1F* uphist, TH1F* dnhist);
TCanvas* compareNormalized(std::string histname, TFile* file1, TFile* file2, TFile* file3, int rebin = 1, bool norm = true);

char* pm;         
char* delim;      
char* delimstart; 
char* delimend;   
char* ee;         
char* mm;         
char* em;         
char* e;         
char* m;         

int   width1;
int   width2;
int   linelength;

#endif
