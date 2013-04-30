#include "plotUtilities.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
#include "TObject.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>
#include "TF1.h"

float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

void printLine( int latex ){

  if( latex == 1 ){
    cout << "\\hline" << endl;
  }
  else if( latex == 2 ){
    cout << "*\\hline" << endl;
  }
  else{
    for( int i = 0 ; i < linelength ; ++i ) cout << "-";
    cout << endl;
  }
}

void printHeader(bool dilep){

  if (dilep) {
    cout << delimstart << setw(width1) << "Sample"    << setw(width2)
	 << delim      << setw(width1) << ee          << setw(width2)
	 << delim      << setw(width1) << mm          << setw(width2)
	 << delim      << setw(width1) << em          << setw(width2)
	 << delim      << setw(width1) << "total"     << setw(width2) 
	 << delimend   << endl;
  } else {
    cout << delimstart << setw(width1) << "Sample"    << setw(width2)
	 << delim      << setw(width1) << e          << setw(width2)
	 << delim      << setw(width1) << m          << setw(width2)
	 << delim      << setw(width1) << "total"     << setw(width2) 
	 << delimend   << endl;
  }
}


void print( TH1F* h , string label , bool dilep, bool correlatedError ){

  stringstream se;
  stringstream sm;
  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if (dilep) {
    if( label == "data" ){
      see  << Form( "%.0f" , h->GetBinContent(1) );
      smm  << Form( "%.0f" , h->GetBinContent(2) );
      sem  << Form( "%.0f" , h->GetBinContent(3) );
      stot << Form( "%.0f" , h->GetBinContent(4)       );
    } else if( label == "data/MC" ){
      see  << Form( "%.2f" , h->GetBinContent(1) );
      smm  << Form( "%.2f" , h->GetBinContent(2) );
      sem  << Form( "%.2f" , h->GetBinContent(3) );
      stot << Form( "%.2f" , h->GetBinContent(4)       );
    }else{
      see  << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
      smm  << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
      sem  << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );
    
      float error = 0;
      if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
      else                  error = h->GetBinError(4);
    
      stot << Form( "%.1f" , h->GetBinContent(4)       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(width1) << label      << setw(width2)
	 << delim      << setw(width1) << see.str()  << setw(width2)
	 << delim      << setw(width1) << smm.str()  << setw(width2)
	 << delim      << setw(width1) << sem.str()  << setw(width2)
	 << delim      << setw(width1) << stot.str() << setw(width2)
	 << delimend   << endl;
  }

  else {
    if( label == "data" ){
      se  << Form( "%.0f" , h->GetBinContent(1) );
      sm  << Form( "%.0f" , h->GetBinContent(2) );
      stot << Form( "%.0f" , h->GetBinContent(3)       );
    } else if( label == "data/MC" ){
      se  << Form( "%.2f" , h->GetBinContent(1) );
      sm  << Form( "%.2f" , h->GetBinContent(2) );
      stot << Form( "%.2f" , h->GetBinContent(3)       );
    }else{
      se  << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
      sm  << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
    
      float error = 0;
      if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2);
      else                  error = h->GetBinError(3);
    
      stot << Form( "%.1f" , h->GetBinContent(3)       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(width1) << label      << setw(width2)
	 << delim      << setw(width1) << se.str()  << setw(width2)
	 << delim      << setw(width1) << sm.str()  << setw(width2)
	 << delim      << setw(width1) << stot.str() << setw(width2)
	 << delimend   << endl;
  }
  
}

void printYields( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* dir , bool doData, int latex ){

  initSymbols( latex );
  bool dilep = false;
  if (TString(dir).Contains("cr3") || TString(dir).Contains("cr4")) dilep = true;

  // h_events:
  //  bin 1: unweighted event yield
  //  bin 2: weighted event yield
  TString histname = "h_events";
  if (!TString(dir).Length() == 0) histname = Form("%s/",dir) + histname;
  TString histname_e = histname + "_e";
  TString histname_m = histname + "_m";
  TString histname_ee = histname + "_ee";
  TString histname_mm = histname + "_mm";
  TString histname_em = histname + "_em";

  printLine(latex);
  printHeader(dilep);
  printLine(latex);

  TH1F* hyield = 0;
  if (dilep) hyield = new TH1F("hyield","yield",4,0,4);
  else hyield = new TH1F("hyield","yield",3,0,3);
  hyield->Sumw2();
  TH1F* hmctot = 0;

  //----------------------
  // print SM MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < mcfiles.size() ; imc++){

    bool correlatedError = false;

    if( TString(labels[imc]).Contains("TChi") ) continue;

    if (dilep) {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_ee = (TH1F*)mcfiles[imc]->Get(histname_ee.Data());
      TH1F* mchist_mm = (TH1F*)mcfiles[imc]->Get(histname_mm.Data());
      TH1F* mchist_em = (TH1F*)mcfiles[imc]->Get(histname_em.Data());

      fillYieldHist(mchist,hyield,4);
      fillYieldHist(mchist_ee,hyield,1);
      fillYieldHist(mchist_mm,hyield,2);
      fillYieldHist(mchist_em,hyield,3);
    }
    else {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_e = (TH1F*)mcfiles[imc]->Get(histname_e.Data());
      TH1F* mchist_m = (TH1F*)mcfiles[imc]->Get(histname_m.Data());

      fillYieldHist(mchist,hyield,3);
      fillYieldHist(mchist_e,hyield,1);
      fillYieldHist(mchist_m,hyield,2);
    }

    if( imc == 0 ) hmctot = (TH1F*) hyield->Clone();
    else           hmctot->Add(hyield);
    
    print( hyield , labels[imc] , dilep, correlatedError );

    //hyield->Reset();
  }

  printLine(latex);

  //-------------------------------
  // print sum of SM MC samples
  //-------------------------------

  print( hmctot , "total SM MC", dilep );

  printLine(latex);
 
  if (doData) {

    if (dilep) {
      TH1F* datahist = (TH1F*)datafile->Get(histname.Data());
      TH1F* datahist_ee = (TH1F*)datafile->Get(histname_ee.Data());
      TH1F* datahist_mm = (TH1F*)datafile->Get(histname_mm.Data());
      TH1F* datahist_em = (TH1F*)datafile->Get(histname_em.Data());

      fillYieldHist(datahist,hyield,4);
      fillYieldHist(datahist_ee,hyield,1);
      fillYieldHist(datahist_mm,hyield,2);
      fillYieldHist(datahist_em,hyield,3);
    }
    else {
      TH1F* datahist = (TH1F*)datafile->Get(histname.Data());
      TH1F* datahist_e = (TH1F*)datafile->Get(histname_e.Data());
      TH1F* datahist_m = (TH1F*)datafile->Get(histname_m.Data());

      fillYieldHist(datahist,hyield,3);
      fillYieldHist(datahist_e,hyield,1);
      fillYieldHist(datahist_m,hyield,2);
    }

    print( hyield , "data", dilep );
    
    printLine(latex);

    TH1F* hratio = (TH1F*)hyield->Clone("hratio");
    hratio->Divide(hmctot);

    print( hratio , "data/MC", dilep );
    
    printLine(latex);
  }

  //----------------------
  // print SUSY MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < mcfiles.size() ; imc++){

    if( !TString(labels[imc]).Contains("TChi") )   continue;

    if (dilep) {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_ee = (TH1F*)mcfiles[imc]->Get(histname_ee.Data());
      TH1F* mchist_mm = (TH1F*)mcfiles[imc]->Get(histname_mm.Data());
      TH1F* mchist_em = (TH1F*)mcfiles[imc]->Get(histname_em.Data());

      fillYieldHist(mchist,hyield,4);
      fillYieldHist(mchist_ee,hyield,1);
      fillYieldHist(mchist_mm,hyield,2);
      fillYieldHist(mchist_em,hyield,3);

    }
    else {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_e = (TH1F*)mcfiles[imc]->Get(histname_e.Data());
      TH1F* mchist_m = (TH1F*)mcfiles[imc]->Get(histname_m.Data());

      fillYieldHist(mchist,hyield,3);
      fillYieldHist(mchist_e,hyield,1);
      fillYieldHist(mchist_m,hyield,2);
    }

    print( hyield , labels[imc], dilep );

  }

  printLine(latex);

  delete hyield;
  delete hmctot;
}

void fillYieldHist( TH1F* hin, TH1F* hyield, int bin ) {

  Double_t err = 0.;
  if (hin) {
    // float nmc = hin->IntegralAndError(0,-1,err);
    // hyield->SetBinContent(bin,nmc);
    // hyield->SetBinError(bin,err);
    // these assume h_events with bin 2 having the weighted event yields
    hyield->SetBinContent(bin,hin->GetBinContent(2));
    hyield->SetBinError(bin,hin->GetBinError(2));
  } else {
    hyield->SetBinContent(bin,0.);
    hyield->SetBinError(bin,0.);
  }

}



#include <TList.h>
#include <TIterator.h>

void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();
   while ((obj=iter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
}

void initSymbols( int latex ){

  //-------------------------------------------------------
  // table format
  //-------------------------------------------------------

  width1      = 22;
  width2      = 5;
  linelength  = (width1+width2)*5+1;

  //-------------------------------------------------------
  // symbols
  //-------------------------------------------------------
  
  if( latex == 1 ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "ee";
    mm         = "$\\mu\\mu$";
    em         = "e$\\mu$";
    e         = "e";
    m         = "$\\mu$";
  } else if (latex == 2) {
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "*";
    delimend   = "\\\\";
    ee         = "ee";
    mm         = "$\\mu\\mu$";
    em         = "e$\\mu$";
    e         = "e";
    m         = "$\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    ee         = "ee";
    mm         = "mm";
    em         = "em";
    e         = "e";
    m         = "m";
  }

}


TLegend *getLegend( vector<char*> labels , bool overlayData, float x1, float y1, float x2, float y2){

  //int colors[]={6,2,7,4,5,8,9,15,12};
  //  int colors[]={4,7,2,5,8,9,15,12,1};
  int colors[]={7,4,2,5,8,9,15,12,1};

  int sigcolors[]={2,4,3};  

  TLegend *leg = new TLegend(x1,y1,x2,y2);

  TH1F*    datahist = new TH1F("datahist","datahist",1,0,1);
  datahist->Sumw2();

  if( overlayData ) leg->AddEntry(datahist,"data");

  const unsigned int nmc = labels.size();
  TH1F*    mchist[nmc];

  //-----------------
  // SM samples
  //-----------------

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( TString(t).Contains("TChiwh") )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    // formatting
    mchist[imc]->SetFillColor( colors[imc] );

    if( strcmp("tt",t)      == 0 ) t = "t#bar{t}";
    if( strcmp("ttbar",t)   == 0 ) t = "t#bar{t}";
    if( strcmp("ttbar 1l",t)== 0 ) t = "t#bar{t} #rightarrow 1l";
    if( strcmp("ttbar 2l",t)== 0 ) t = "t#bar{t} #rightarrow 2l";
    if( strcmp("ttll",t)    == 0 ) t = "t#bar{t} #rightarrow ll";
    if( strcmp("tttau",t)   == 0 ) t = "t#bar{t} #rightarrow l#tau/#tau#tau";
    if( strcmp("ttfake",t)  == 0 ) t = "t#bar{t} #rightarrow fake";
    if( strcmp("t",t)       == 0 ) t = "single top";
    if( strcmp("single_top",t)       == 0 ) t = "single top";
    if( strcmp("wjets",t)   == 0 ) t = "W+jets";
    if( strcmp("zjets",t)   == 0 ) t = "Z+jets";
    if( strcmp("ww",t)      == 0 ) t = "WW";
    if( strcmp("wz",t)      == 0 ) t = "WZ";
    if( strcmp("zz",t)      == 0 ) t = "ZZ";
    if( strcmp("whbb",t)      == 0 ) t = "WH #rightarrow l#nubb";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  //-----------------
  // sig samples
  //-----------------

  int nsigmc = 0;

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //  for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( !TString(t).Contains("TChiwh") )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    // formatting
    mchist[imc]->SetFillColor( 0 );
    mchist[imc]->SetLineColor(sigcolors[nsigmc]);
    mchist[imc]->SetLineWidth(2);
    //    mchist[imc]->SetLineStyle(2);
    ++nsigmc;

    if( strcmp("ttall",t) == 0 ) t = "t#bar{t}";
    if( strcmp("t",t)     == 0 ) t = "single top";
    if( strcmp("wjets",t) == 0 ) t = "W+jets";
    if( strcmp("WW",t)    == 0 ) t = "W^{+}W^{-}";
    if( strcmp("WZ",t)    == 0 ) t = "W^{#pm}Z^{0}";
    if( strcmp("ZZ",t)    == 0 ) t = "Z^{0}Z^{0}";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  
  return leg;

}

//TH1F* compareDataMC( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
TGraphErrors* compareDataMC( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , 
                    const char* flavor , const char* dir ,
		    int nbins ,  float xmin , float xmax ,  
		    const char* xtitle , bool overlayData , bool residual , bool drawLegend , bool log , 
			     bool normalize , bool fit, float mcnorm ){

  dataMCHists datamchists;
  TH1F* ratio = 0;
  TGraphErrors* gratio = 0;

  TString fullhistname = histname;
  if (!TString(dir).Length() == 0) fullhistname = Form("%s/",dir) + fullhistname;

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->Draw();
    plotpad->cd();
    if( log ) plotpad->SetLogy();
  }
  else{
    if( log ) gPad->SetLogy();
  }

  cout << "Plotting " << fullhistname << endl;

  //int colors[]={6,2,7,4,5,8,9,15,12};
  //int colors[]={kRed+2,5,7,5,5,8,9,15,12};
  // zmet colors
  //  int colors[]={4,7,2,5,8,9,15,12,1};

  //  int colors[]={4,7,2,5,8,9,15,12,1};
  int colors[]={7,4,2,5,8,9,15,12,1};

  int sigcolors[]={2,4,3};

  assert( mcfiles.size() == labels.size() );
  const unsigned int nmc = mcfiles.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = 0;
  vector<TH1F*> mchist;
  mchist.resize(nmc);
  TH1F* datahist = 0;
  if (overlayData) {
    TH1F* datahist_tmp = (TH1F*)datafile->Get(fullhistname.Data());
    if (datahist_tmp) {
      datahist = (TH1F*)datahist_tmp->Clone(fullhistname+TString("_data"));
    } else {
      std::cout << "WARNING: no hist with name " << fullhistname << " found for data" << std::endl;
      datahist = 0;
    }
  }
  float nmctot = 0.0;
  float errmctot = 0.0;

  for( int imc = nmc-1 ; imc > -1 ; imc-- ){
    TH1F* mchist_tmp = (TH1F*)mcfiles[imc]->Get(fullhistname.Data());
    if (mchist_tmp) {
      mchist[imc] = (TH1F*)mchist_tmp->Clone(fullhistname+Form("_%s",labels[imc]));
      Double_t err;
      float nmc = mchist[imc]->IntegralAndError(0,-1,err);
      if (TString(labels[imc]).Contains("TChiwh")) {
	// signal: don't add to bg total
	cout << "MC signal yield " << labels[imc] << " " << Form("%.2f",nmc) << " +/- " << Form("%.2f",err) << endl;
      } else {
	nmctot += mchist[imc]->Integral();
	Double_t err;
	float nmc = mchist[imc]->IntegralAndError(0,-1,err);
	errmctot = sqrt(pow(errmctot,2) + pow(err,2));
	cout << "MC yield " << labels[imc] << " " << Form("%.2f",nmc) << " +/- " << Form("%.2f",err) << endl;
      }
    } else {
      std::cout << "WARNING: no hist with name " << fullhistname << " found for sample " << labels[imc] << std::endl;
      mchist[imc] = 0;
    }
  }

  cout << "MC bg total yield " << Form("%.2f",nmctot) << " +/- " << Form("%.2f",errmctot) << endl;

  float ndata = 0.;

  if (overlayData) {
    ndata = datahist->Integral();
    cout << "data yield " << ndata << endl;
  }

  float SF = ndata/nmctot;
  if( normalize ){
    cout << "Data, MC, SF " << ndata << ", " << nmctot << ", " << SF << endl;
    if (mcnorm < 0.) cout << "Scaling MC by " << SF << "!!!!!!" << endl;
    else cout << "Scaling MC by input factor of " << mcnorm << "!!!!!!" << endl;
  }

  int nmcsig = 0;
  std::vector<TH1F*> mcsighist;

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc-1 ; imc > -1 ; imc-- ){

    // skip samples without hists
    if (!mchist[imc]) continue;

    mchist[imc]->Sumw2();
    float xmin_original = mchist[imc]->GetXaxis()->GetXmin();
    float xmax_original = mchist[imc]->GetXaxis()->GetXmax();
    int nbins_original = mchist[imc]->GetNbinsX();
    float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    float xbinsize = (xmax - xmin)/float(nbins);
    int rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) mchist[imc]->Rebin(rebinFactor);
    mchist[imc]->SetAxisRange(xmin,xmax-0.1*xbinsize,"X");
    Double_t err;
    float lastbin_and_overflow =  mchist[imc]->IntegralAndError(nbins,-1,err);
    mchist[imc]->SetBinContent(nbins,lastbin_and_overflow);
    mchist[imc]->SetBinError(nbins,err);

    if( normalize ) {
      if (mcnorm < 0.) mchist[imc]->Scale(SF);
      else mchist[imc]->Scale(mcnorm);
    }

    if( TString( labels.at(imc) ).Contains("TChiwh") ){
      // signal
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineWidth(2);
      //      mchist[imc]->SetLineStyle(2);
      ++nmcsig;
      mcsighist.push_back(mchist[imc]);
    }else{
      // bg
      mchist[imc]->SetLineColor(kBlack);
      mchist[imc]->SetFillColor( colors[imc] );
      mcstack->Add( mchist[imc] );

      //if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
      //    if( imc == nmc-1 ) mctothist = (TH1F*) mchist[imc]->Clone();
      if( !mctothist ) mctothist = (TH1F*) mchist[imc]->Clone();
      else               mctothist->Add(mchist[imc]);

    }

  }

  mctothist->SetLineColor(kBlack);

  if( overlayData ){

    datahist->Sumw2();
    float xmin_original = datahist->GetXaxis()->GetXmin();
    float xmax_original = datahist->GetXaxis()->GetXmax();
    int nbins_original = datahist->GetNbinsX();
    float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    float xbinsize = (xmax - xmin)/float(nbins);
    int rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) datahist->Rebin(rebinFactor);
    datahist->SetAxisRange(xmin,xmax-0.1*xbinsize,"X");
    Double_t err;
    float lastbin_and_overflow =  datahist->IntegralAndError(nbins,-1,err);
    datahist->SetBinContent(nbins,lastbin_and_overflow);
    datahist->SetBinError(nbins,err);

    float max = datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin());
    if( mctothist->GetMaximum() > max ) max = mctothist->GetMaximum();
    if( log ) datahist->SetMaximum( 90 * max );
    else      datahist->SetMaximum( 1.4 * max );

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetXaxis()->SetTitleSize(0.05);
    datahist->GetXaxis()->SetLabelSize(0.04);
    datahist->GetXaxis()->SetTitleOffset(1.2);
    datahist->SetLineColor(kBlack);
    datahist->SetMarkerColor(kBlack);
    datahist->SetMarkerStyle(20);
    datahist->Draw("E1");
    mcstack->Draw("samehist");
    if (nmcsig > 0) {
      for (int i = (int) mcsighist.size() - 1; i > -1 ; i--) {
        mcsighist[i]->SetLineColor( sigcolors[mcsighist.size()-1-i] );
	mcsighist[i]->Draw("samehist");
      }
    }
    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");
    
    if(!log) {
      datahist->GetYaxis()->SetRangeUser(0.,1.4*max);
      if (TString(histname).Contains("dphi")) datahist->GetYaxis()->SetRangeUser(0.,3.0*max);
    }

    datamchists.datahist = datahist;
    
  }
  else{

    float max = mctothist->GetMaximum() + mctothist->GetBinError(mctothist->GetMaximumBin());
    if( log ) mctothist->SetMaximum( 90 * max );
    else      mctothist->SetMaximum( 1.4 * max );

    mctothist->GetXaxis()->SetTitle(xtitle);
    mctothist->Draw("hist");
    mcstack->Draw("hist same");
    mctothist->Draw("hist same axis");
    if (nmcsig > 0) {
      for (int i = (int) mcsighist.size() - 1; i > -1 ; i--) {
        mcsighist[i]->SetLineColor( sigcolors[mcsighist.size()-1-i] );
	mcsighist[i]->Draw("samehist");
      }
    }

    if(!log) mctothist->GetYaxis()->SetRangeUser(0.,1.4*max);
  }

  datamchists.mctothist = mctothist;

  if( drawLegend ){
    TLegend* myleg = getLegend( labels , overlayData );
    myleg->Draw();
  }

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  text->DrawLatex(0.2,0.88,"CMS Preliminary");
  //text->DrawLatex(0.2,0.83,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
  text->DrawLatex(0.2,0.83,"#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  if     ( TString(flavor).Contains("ee")  ) text->DrawLatex(0.2,0.78,"Events with ee");
  else if( TString(flavor).Contains("mm")  ) text->DrawLatex(0.2,0.78,"Events with #mu#mu");
  else if( TString(flavor).Contains("em")  ) text->DrawLatex(0.2,0.78,"Events with e#mu");
  else if( TString(flavor).Contains("all") ) text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu/e#mu");
  else if( TString(flavor).Contains("sf")  ) text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu");
  else if( TString(flavor).Contains("3l")  ) text->DrawLatex(0.2,0.78,"Events with eee/ee#mu/#mu#mue/#mu#mu#mu");
  else if( TString(flavor).Contains("e")  ) text->DrawLatex(0.2,0.78,"Events with e");
  else if( TString(flavor).Contains("m")  ) text->DrawLatex(0.2,0.78,"Events with #mu");
  else if( TString(flavor).Contains("sl")  ) text->DrawLatex(0.2,0.78,"Events with e/#mu");
  else                                       text->DrawLatex(0.2,0.78,"Events with e/#mu");
  //else                                       text->DrawLatex(0.2,0.78,"Events with e#mu");
  //else                                       text->DrawLatex(0.2,0.78,"Events with #mu#mu");

  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,0.98);
    respad->Draw();
    respad->cd();
    respad->SetTopMargin(0.05);

    gPad->SetGridy();

    //    if (false) {
    if ((TString(histname).Contains("sumpt") || TString(histname).Contains("dilpt")) 
        && !TString(histname).Contains("rapidity") ) {
      // int newnbins = 18;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,520};

      // int newnbins = 13;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,350,500};

      // binning for ttbar high stats sumpt plots
      int newnbins = 10;
      Double_t xbins[] = {0,20,40,60,80,100,150,200,250,350,500};

      // used for most plots of sumpt and dilpt
      //int newnbins = 14;
      //Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,350,500,800};

      // new z+jets binning for flatter systs
      // int newnbins = 15;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,400,500,800};

      // for monojet dilpt plot
      // int newnbins = 7;
      // Double_t xbins[] = {0,150,200,250,300,350,500,800};

      // int newnbins = 6;
      // Double_t xbins[] = {0,50,100,150,200,300,520};

      // int newnbins = 8;
      // Double_t xbins[] = {0,20,40,60,80,100,150,250,520};

      // int newnbins = 3;
      // Double_t xbins[] = {0,100,200,520};

      // int newnbins = 2;
      // Double_t xbins[] = {0,150,500};

      ratio = (TH1F*) datahist->Rebin(newnbins,Form("%s_ratio",datahist->GetName()),xbins);
      TH1F* mcdenom = (TH1F*) mctothist->Rebin(newnbins,Form("%s_denom",mctothist->GetName()),xbins);
      ratio->Divide(mcdenom);
    } else {
      ratio = (TH1F*) datahist->Clone(Form("%s_ratio",datahist->GetName()));
      ratio->Divide(mctothist);
    }

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.15);
    //ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetTitle("data/MC  ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(0.7);

    TF1* fpol1 = new TF1("fpol1","pol1",datahist->GetXaxis()->GetXmin(),datahist->GetXaxis()->GetXmax());
    if( fit ) ratio->Fit(fpol1,"R");

    ratio->Draw();

    gratio = new TGraphErrors(ratio);
    for (int i=0; i<ratio->GetNbinsX(); ++i) {
      gratio->SetPointError(i,ratio->GetBinWidth(i+1)/2.,gratio->GetErrorY(i));
    }
    gratio->SetLineColor(kBlack);
    gratio->SetMarkerColor(kBlack);
    gratio->SetMarkerStyle(20);
    gratio->Draw("psame");

    TLine line;
    line.SetLineWidth(1);
    //    line.DrawLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
    line.DrawLine(xmin,1,xmax,1);

    datamchists.ratio = ratio;
  }

  datamchists.comppad = fullpad;
  //  return datamchists;
  //  return ratio;
  return gratio;
}

//____________________________________________________________________________
TGraphAsymmErrors* makeBand(TH1F* centhist, TH1F* uphist, TH1F* dnhist) {

  int nbins = centhist->GetNbinsX();

  // makes a graph with the correct central values
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nbins);

  // need to set correct upper and lower errors to create band

  // assume all bins are same width..
  float binwidth = centhist->GetBinWidth(1);

  for (int ibin=1; ibin < nbins+1; ++ibin) {
    float valx = centhist->GetBinCenter(ibin);
    float valy = centhist->GetBinContent(ibin);
    graph->SetPoint(ibin-1,valx,valy);
    float uperr = valy - uphist->GetBinContent(ibin);
    float dnerr = dnhist->GetBinContent(ibin) - valy;
    float err = max(uperr,dnerr);
    if (uperr < 0.) uperr = dnerr;
    else if (dnerr < 0.) dnerr = uperr;
    //    graph->SetPointError(ibin-1,binwidth/2.,binwidth/2.,err,err);
    graph->SetPointError(ibin-1,binwidth/2.,binwidth/2.,uperr,dnerr);
  }

  return graph;
}

//____________________________________________________________________________
TGraphAsymmErrors* makeBand(TGraphErrors* centgraph, TGraphErrors* upgraph, TGraphErrors* dngraph) {

  int nbins = centgraph->GetN();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nbins);

  for (int ibin=0; ibin < nbins; ++ibin) {
    Double_t valx, valy;
    centgraph->GetPoint(ibin,valx,valy);
    graph->SetPoint(ibin,valx,valy);
    Double_t errx = centgraph->GetErrorX(ibin);
    Double_t upy, dny;
    upgraph->GetPoint(ibin,valx,upy);
    dngraph->GetPoint(ibin,valx,dny);
    float uperr = valy - upy;
    float dnerr = dny - valy;
    float err = max(uperr,dnerr);
    if (uperr < 0.) uperr = dnerr;
    else if (dnerr < 0.) dnerr = uperr;
    //    graph->SetPointError(ibin,errx,errx,err,err);
    graph->SetPointError(ibin,errx,errx,uperr,dnerr);
  }

  return graph;
}

// //____________________________________________________________________________
// TCanvas* compareNormalized(std::string histname, TFile* fsig, TFile* fsig2, TFile* fbg, int rebin, bool norm) {

//   TCanvas* c = new TCanvas(Form("c_%s",histname.c_str()),Form("c_%s",histname.c_str()));

//   TH1F* hsig = (TH1F*)((TH1F*)fsig->Get(histname.c_str()))->Clone(Form("%s_fsig",histname.c_str()));
//   TH1F* hsig2 = 0;
//   if (fsig2) hsig2 = (TH1F*)((TH1F*)fsig2->Get(histname.c_str()))->Clone(Form("%s_fsig2",histname.c_str()));
//   TH1F* hbg = (TH1F*)((TH1F*)fbg->Get(histname.c_str()))->Clone(Form("%s_fbg",histname.c_str()));

//   hsig->SetMarkerColor(kBlue);
//   hsig->SetLineColor(kBlue);
//   hsig->SetLineWidth(2);
//   if (hsig2) {
//     hsig2->SetMarkerColor(kGreen);
//     hsig2->SetLineColor(kGreen);
//     hsig2->SetLineWidth(2);
//   }
//   hbg->SetMarkerColor(kRed);
//   hbg->SetLineColor(kRed);
//   hbg->SetLineWidth(2);

//   if (rebin > 1) {
//     hsig->Rebin(rebin);
//     if (hsig2) hsig2->Rebin(rebin);
//     hbg->Rebin(rebin);
//   }

//   TH1F* hsig_norm = 0;
//   TH1F* hsig2_norm = 0;
//   TH1F* hbg_norm = 0;

//   if (norm) {
//     hsig_norm = (TH1F*)hsig->DrawNormalized("histe");
//     if (hsig2) hsig2_norm = (TH1F*)hsig2->DrawNormalized("histe same");
//     hbg_norm = (TH1F*)hbg->DrawNormalized("histe same");
//   } else {
//     hsig_norm = hsig;
//     if (hsig2) hsig2_norm = hsig2;
//     hbg_norm = hbg;
//     hsig_norm->Draw("histe");
//     if (hsig2) hsig2_norm->Draw("histe same");
//     hbg_norm->Draw("histe same");
//   }

//   if (hbg_norm->GetMaximum() > hsig_norm->GetMaximum()) {
//     hsig_norm->GetYaxis()->SetRangeUser(1E-4,1.1*hbg_norm->GetMaximum());
//   }

//   if (hsig2_norm && (hsig2_norm->GetMaximum() > hsig_norm->GetMaximum())) {
//     hsig_norm->GetYaxis()->SetRangeUser(1E-4,1.1*hsig2_norm->GetMaximum());
//   }

//   // std::string label_sig = "Wbb+jets";
//   // std::string label_sig2 = "TChiwh_150_1";
//   // std::string label_bg = "W+jets";
//   std::string label_sig = "TChiwh_250_1";
//   std::string label_sig2 = "TChiwh_150_1";
//   std::string label_bg = "all bkg";

//   TLegend *leg = new TLegend(0.66,0.77,0.93,0.9);
//   leg->SetFillColor(0);
//   leg->AddEntry(hsig_norm,label_sig.c_str(),"l");
//   if (hsig2_norm) leg->AddEntry(hsig2_norm,label_sig2.c_str(),"l");
//   leg->AddEntry(hbg_norm,label_bg.c_str(),"l");
//   leg->Draw("same");

//   gPad->Modified();

//   return c;
// }

//____________________________________________________________________________
TCanvas* compareNormalized(std::string histname, TFile* f1, std::string label1, TFile* f2, std::string label2, int rebin, bool norm, TFile* f3, std::string label3) {

  TH1F* h1 = (TH1F*)f1->Get(histname.c_str());
  TH1F* h2 = (TH1F*)f2->Get(histname.c_str());
  TH1F* h3 = 0;
  if (f3) h3 = (TH1F*)f3->Get(histname.c_str());

  TCanvas* c = compareNormalized(h1,label1,h2,label2,rebin,norm,h3,label3);

  return c;
}

//____________________________________________________________________________
TCanvas* compareNormalized(TH1F* h1, std::string label1, TH1F* h2, std::string label2, int rebin, bool norm, TH1F* h3, std::string label3) {

  TCanvas* c = new TCanvas(Form("c_%s",h1->GetName()),Form("c_%s",h1->GetName()));

  TH1F* h1_clone = (TH1F*)h1->Clone(Form("%s_clone",h1->GetName()));
  TH1F* h2_clone = (TH1F*)h2->Clone(Form("%s_clone",h2->GetName()));
  TH1F* h3_clone = 0;
  if (h3) h3_clone = (TH1F*)h3->Clone(Form("%s_clone",h3->GetName()));

  h1_clone->SetMarkerColor(kBlue);
  h1_clone->SetLineColor(kBlue);
  h1_clone->SetLineWidth(2);
  h2_clone->SetMarkerColor(kRed);
  h2_clone->SetLineColor(kRed);
  h2_clone->SetLineWidth(2);
  if (h3_clone) {
    h3_clone->SetMarkerColor(7);
    h3_clone->SetLineColor(7);
    h3_clone->SetLineWidth(2);
  }

  if (rebin > 1) {
    h1_clone->Rebin(rebin);
    if (h3_clone) h3_clone->Rebin(rebin);
    h2_clone->Rebin(rebin);
  }

  TH1F* h1_norm = 0;
  TH1F* h2_norm = 0;
  TH1F* h3_norm = 0;

  if (norm) {
    h1_norm = (TH1F*)h1_clone->DrawNormalized("histe");
    h2_norm = (TH1F*)h2_clone->DrawNormalized("histe same");
    if (h3_clone) h3_norm = (TH1F*)h3_clone->DrawNormalized("histe same");
  } else {
    h1_norm = h1_clone;
    h2_norm = h2_clone;
    if (h3_clone) h3_norm = h3_clone;
    h1_norm->Draw("histe");
    h2_norm->Draw("histe same");
    if (h3_clone) h3_norm->Draw("histe same");
  }

  if (h2_norm->GetMaximum() > h1_norm->GetMaximum()) {
    h1_norm->GetYaxis()->SetRangeUser(1E-4,1.1*h2_norm->GetMaximum());
  }

  if (h3_norm && (h3_norm->GetMaximum() > h1_norm->GetMaximum())) {
    h1_norm->GetYaxis()->SetRangeUser(1E-4,1.1*h3_norm->GetMaximum());
  }

  TLegend *leg = new TLegend(0.66,0.77,0.93,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(h1_norm,label1.c_str(),"l");
  leg->AddEntry(h2_norm,label2.c_str(),"l");
  if (h3_norm) leg->AddEntry(h3_norm,label3.c_str(),"l");
  leg->Draw("same");

  gPad->Modified();

  return c;
}

// -----------------------------------------------

TH1F* cumulate (TH1F* in, bool increasing) {
  TH1F* h_out = new TH1F(in->GetName() + TString("tmp"), in->GetTitle(), in->GetNbinsX(),
			 in->GetBinLowEdge(1), in->GetBinLowEdge(in->GetNbinsX() + 1));
  h_out->Sumw2();
  h_out->SetFillColor(in->GetFillColor());
  h_out->SetFillStyle(in->GetFillStyle());
  h_out->SetLineStyle(in->GetLineStyle());
  h_out->SetLineColor(in->GetLineColor());
  h_out->SetLineWidth(in->GetLineWidth());
  double sum = 0;
  double err2 = 0;
  if (increasing) {
    for (int j = 0; j <= in->GetNbinsX() + 1; ++j) {
      sum += in->GetBinContent(j);
      err2 += in->GetBinError(j)*in->GetBinError(j);
      h_out->SetBinContent(j, sum);
      h_out->SetBinError(j, sqrt(err2));
    }
  } else {
    for (int j = in->GetNbinsX() + 1; j >= 0; --j) {
      sum += in->GetBinContent(j);
      err2 += in->GetBinError(j)*in->GetBinError(j);
      h_out->SetBinContent(j, sum);
      h_out->SetBinError(j, sqrt(err2));
    }
  }
  return h_out;
}

// ---------------------------------------------------

TGraphErrors* eff_rej (TH1F* signal, TH1F* background, bool normalize, bool increasing, bool print) {

  TH1F* sig = (TH1F*)signal->Clone("h_tmp_s");
  if (normalize) sig->Scale(1 / sig->Integral(0, sig->GetNbinsX() + 1));
  TH1F* bg = (TH1F*)background->Clone("h_tmp_bg");
  if (normalize) bg->Scale(1 / bg->Integral(0, bg->GetNbinsX() + 1));
  TH1F* sig_cum = cumulate(sig, increasing);
  TH1F* bg_cum = cumulate(bg, increasing);
  TGraphErrors* ret = new TGraphErrors(signal->GetNbinsX());
  for (int i = 1; i <= signal->GetNbinsX(); ++i) {
    const double x = sig_cum->GetBinCenter(i);
    const double s = sig_cum->GetBinContent(i);
    if (!normalize && (s < 3.0)) continue;
    const double b = bg_cum->GetBinContent(bg_cum->FindBin(x));
    const double errs = sig_cum->GetBinError(i);
    const double errb = bg_cum->GetBinError(bg_cum->FindBin(x));
    ret->SetPoint(i - 1, b, s); // gotta love offsets
    ret->SetPointError(i - 1, errb, errs);
    if (print) {
      std::cout << "x value: " << x << ", s: " << s << ", b: " << b;
      if (b > 0.) std::cout << ", s/b: " << s/b << ", s/sqrt(b): " << s/sqrt(b)
			    << ", sig: " << sqrt(2*( (s+b) * log(1+(s/b)) - s) );
      std::cout << std::endl;
    }
  }
  return ret;
}

// ---------------------------------------------------
// function to make graph of S/sqrt(B)

TGraph* s_over_rootb (TH1F* signal, TH1F* background, bool increasing, bool do_s_over_b, bool print) {

  TH1F* sig = (TH1F*)signal->Clone("h_tmp_s");
  TH1F* bg = (TH1F*)background->Clone("h_tmp_bg");
  TH1F* sig_cum = cumulate(sig, increasing);
  TH1F* bg_cum = cumulate(bg, increasing);
  TGraph* ret = new TGraphErrors(signal->GetNbinsX());
  for (int i = 1; i <= signal->GetNbinsX(); ++i) {
    const double x = sig_cum->GetBinCenter(i);
    const double s = sig_cum->GetBinContent(i);
    if (s < 3.0) continue;
    const double b = bg_cum->GetBinContent(bg_cum->FindBin(x));
    //    const double errs = sig_cum->GetBinError(i);
    //    const double errb = bg_cum->GetBinError(bg_cum->FindBin(x));
    double s_over_rootb = 0.;
    if (b > 0.) {
      s_over_rootb = s/sqrt(b);
      if (do_s_over_b) s_over_rootb = s/b;
    }
    ret->SetPoint(i - 1, x, s_over_rootb); // gotta love offsets
    //    ret->SetPointError(i - 1, errb, errs);
    if (print) {
      std::cout << "x value: " << x << ", s: " << s << ", b: " << b;
      if (b > 0.) std::cout << ", s/b: " << s/b << ", s/sqrt(b): " << s/sqrt(b)
			    << ", sig: " << sqrt(2*( (s+b) * log(1+(s/b)) - s) );
      std::cout << std::endl;
    }
  }
  return ret;
}

//____________________________________________________________________
TLegend* init_legend(float x1, float y1, float x2, float y2) {
  // intializes a legend with white background, no border
  TLegend* l = new TLegend(x1,y1,x2,y2);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->SetShadowColor(0);
  return l;
}

//______________________________________________________________________________
TLegend* legendize(TCanvas* c, const TString& opt, const TString& label1, const TString& label2, const TString& label3 ) {

  TString labels[3] = {label1,label2,label3};
  // makes a basic legend by finding all TH1s and graphs in a canvas
  TList* clist = c->GetListOfPrimitives();
  TLegend* l = init_legend();
  // iterate over list of primitives
  int nhists = 0;
  TIter next(clist);
  TObject* obj = 0;
  while ((obj = next())) {
    // check for TH1 or TGraph
    if (TString(obj->ClassName()).Contains("TH1")) {
      TH1* h = (TH1*)obj;
      if (nhists <= 2) {
	l->AddEntry(h,labels[nhists],opt);
      } else {
	l->AddEntry(h,"",opt);
      }
      nhists++;
    } else if (TString(obj->ClassName()).Contains("TGraph")) {
      TGraph* h = (TGraph*)obj;
      if (nhists <= 2) {
	l->AddEntry(h,labels[nhists],opt);
      } else {
	l->AddEntry(h,"",opt);
      }
      nhists++;
    }
  }
  l->Draw("same");

  return l;
}
