#include "plotUtilities.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
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

void printLine( bool latex ){

  if( latex ){
    cout << "\\hline" << endl;
  }
  else{
    for( int i = 0 ; i < linelength ; ++i ) cout << "-";
    cout << endl;
  }
}

void printHeader(){

  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee          << setw(width2)
       << delim      << setw(width1) << mm          << setw(width2)
       << delim      << setw(width1) << em          << setw(width2)
       << delim      << setw(width1) << "total"     << setw(width2) 
       << delimend   << endl;

}


void print( TH1F* h , string label , bool correlatedError = false ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if( label == "data" ){
    see  << Form( "%.0f" , h->GetBinContent(1) );
    smm  << Form( "%.0f" , h->GetBinContent(2) );
    sem  << Form( "%.0f" , h->GetBinContent(3) );
    stot << Form( "%.0f" , h->Integral()       );
  }else{
    //see  << Form( "%.1f" , h->GetBinContent(1) );
    //smm  << Form( "%.1f" , h->GetBinContent(2) );
    //sem  << Form( "%.1f" , h->GetBinContent(3) );
    //stot << Form( "%.1f" , h->Integral()       );
    
    see  << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
    smm  << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
    sem  << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );
    
    float error = 0;
    if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
    else                  error = histError(h,1,4);
    
    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << see.str()  << setw(width2)
       << delim      << setw(width1) << smm.str()  << setw(width2)
       << delim      << setw(width1) << sem.str()  << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
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

void initSymbols( bool latex ){

  //-------------------------------------------------------
  // table format
  //-------------------------------------------------------

  width1      = 15;
  width2      = 4;
  linelength  = (width1+width2)*5+1;

  //-------------------------------------------------------
  // symbols
  //-------------------------------------------------------
  
  if( latex ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "ee";
    mm         = "$\\mu\\mu$";
    em         = "e$\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    ee         = "ee";
    mm         = "mm";
    em         = "em";
  }

}


TLegend *getLegend( vector<char*> labels , bool overlayData, float x1, float y1, float x2, float y2){

  //int colors[]={6,2,7,4,5,8,9,15,12};
  int colors[]={4,7,2,5,8,9,15,12,1};
  
  TLegend *leg = new TLegend(x1,y1,x2,y2);

  TH1F*    datahist = new TH1F("datahist","datahist",1,0,1);
  datahist->Sumw2();

  if( overlayData ) leg->AddEntry(datahist,"data");

  const int nmc = labels.size();
  TH1F*    mchist[nmc];

  //-----------------
  // SM samples
  //-----------------

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( TString(t).Contains("LM") )   continue;
    if( TString(t).Contains("T2tt") ) continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("LM") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    if( strcmp("tt",t)      == 0 ) t = "t#bar{t}";
    if( strcmp("ttbar",t)   == 0 ) t = "t#bar{t}";
    if( strcmp("ttll",t)    == 0 ) t = "t#bar{t} #rightarrow ll";
    if( strcmp("tttau",t)   == 0 ) t = "t#bar{t} #rightarrow l#tau/#tau#tau";
    if( strcmp("ttfake",t)  == 0 ) t = "t#bar{t} #rightarrow fake";
    if( strcmp("t",t)       == 0 ) t = "single top";
    if( strcmp("wjets",t)   == 0 ) t = "W+jets";
    if( strcmp("zjets",t)   == 0 ) t = "Z+jets";
    if( strcmp("ww",t)      == 0 ) t = "WW";
    if( strcmp("wz",t)      == 0 ) t = "WZ";
    if( strcmp("zz",t)      == 0 ) t = "ZZ";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  //-----------------
  // LM samples
  //-----------------

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( ! (TString(t).Contains("LM")||TString(t).Contains("T2tt")) )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("T2tt") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

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
TGraphErrors* compareDataMC( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
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

  int colors[]={4,7,2,5,8,9,15,12,1};

  assert( mcfiles.size() == labels.size() );
  const unsigned int nmc = mcfiles.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = new TH1F();
  vector<TH1F*> mchist;
  mchist.resize(nmc);
  TH1F* datahist = 0;
  if (overlayData) {
    datahist = (TH1F*)((TH1F*)datafile->Get(fullhistname.Data()))->Clone(fullhistname+TString("_data"));
  }
  float nmctot = 0.0;

  for( int imc = nmc-1 ; imc > -1 ; imc-- ){
    mchist[imc] = (TH1F*)((TH1F*)mcfiles[imc]->Get(fullhistname.Data()))->Clone(fullhistname+Form("_%s",labels[imc]));
    nmctot += mchist[imc]->Integral();
    cout << "MC yield " << labels[imc] << " " << Form("%.2f",mchist[imc]->Integral()) << endl;
  }

  cout << "MC total yield " << Form("%.2f",nmctot) << endl;

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

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc-1 ; imc > -1 ; imc-- ){

    mchist[imc]->Sumw2();
    float xmin_original = mchist[imc]->GetXaxis()->GetXmin();
    float xmax_original = mchist[imc]->GetXaxis()->GetXmax();
    int nbins_original = mchist[imc]->GetNbinsX();
    float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    float xbinsize = (xmax - xmin)/float(nbins);
    int rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) mchist[imc]->Rebin(rebinFactor);
    mchist[imc]->SetAxisRange(xmin,xmax,"X");

    if( normalize ) {
      if (mcnorm < 0.) mchist[imc]->Scale(SF);
      else mchist[imc]->Scale(mcnorm);
    }

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("T2tt") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    mcstack->Add( mchist[imc] );

    //if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
    if( imc == nmc-1 ) mctothist = (TH1F*) mchist[imc]->Clone();
    else               mctothist->Add(mchist[imc]);

  }

  if( overlayData ){

    datahist->Sumw2();
    float xmin_original = datahist->GetXaxis()->GetXmin();
    float xmax_original = datahist->GetXaxis()->GetXmax();
    int nbins_original = datahist->GetNbinsX();
    float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    float xbinsize = (xmax - xmin)/float(nbins);
    int rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) datahist->Rebin(rebinFactor);
    datahist->SetAxisRange(xmin,xmax,"X");

    float max = datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin());
    if( mctothist->GetMaximum() > max ) max = mctothist->GetMaximum();
    if( log ) datahist->SetMaximum( 90 * max );
    else      datahist->SetMaximum( 1.4 * max );

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetXaxis()->SetTitleSize(0.05);
    datahist->GetXaxis()->SetLabelSize(0.04);
    datahist->GetXaxis()->SetTitleOffset(1.2);
    datahist->Draw("E1");
    mcstack->Draw("samehist");
    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");
    
    if(!log) datahist->GetYaxis()->SetRangeUser(0.,1.4*max);

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
  else                                       text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu");
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
