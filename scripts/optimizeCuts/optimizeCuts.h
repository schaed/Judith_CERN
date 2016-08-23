#ifndef OPTIMIZECUTS_HH
#define OPTIMIZECUTS_HH

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TError.h>
#include <TPad.h>

#define MAX_HITS 3

#define T0NBINS 100
#define T0MIN 100.
#define T0MAX 300.
#define T0FITRANGELOW 180.
#define T0FITRANGEHIGH 240.
#define T0NPAR 6
#define T0PARINITPLATEAU 200.
#define T0PARINITMULEFT 200.
#define T0PARINITSIGMALEFT 2.
#define T0PARINITMURIGHT 225.
#define T0PARINITSIGMARIGHT 2.
#define T0PARINITOFFSET 10.
#define T0PARLIMITLOWPLATEAU 100.
#define T0PARLIMITLOWMULEFT 190.
#define T0PARLIMITLOWSIGMALEFT 1.
#define T0PARLIMITLOWMURIGHT 215.
#define T0PARLIMITLOWSIGMARIGHT 1.
#define T0PARLIMITLOWOFFSET 0.
#define T0PARLIMITHIGHPLATEAU 300.
#define T0PARLIMITHIGHMULEFT 210.
#define T0PARLIMITHIGHSIGMALEFT 3.
#define T0PARLIMITHIGHMURIGHT 235.
#define T0PARLIMITHIGHSIGMARIGHT 3.
#define T0PARLIMITHIGHOFFSET 100.
#define T0LEGENDXLOW 0.58
#define T0LEGENDYLOW 0.78
#define T0LEGENDXHIGH 0.89
#define T0LEGENDYHIGH 0.89
#define T0STATSXLOW  0.13
#define T0STATSYLOW 0.56
#define T0STATSXHIGH 0.42
#define T0STATSYHIGH 0.88
#define T0RESULTSXLOW  0.13
#define T0RESULTSYLOW 0.44
#define T0RESULTSXHIGH 0.42
#define T0RESULTSYHIGH 0.54
#define T0CUTNSIGMA 3.

#define CHARGENBINS 500
#define CHARGEMIN 0.
#define CHARGEMAX 0.5
#define CHARGEFITRANGELOW 0.015
#define CHARGEFITRANGEHIGH 0.25
#define CHARGENPAR 5
#define CHARGEPARINITSCALE 1.5
#define CHARGEPARINITSIGMA 0.04
#define CHARGEPARINITMU 0.1
#define CHARGEPARINITBGSCALE 0.2
#define CHARGEPARINITBGASYN 0.005
#define CHARGEPARLIMITLOWSCALE 0.
#define CHARGEPARLIMITLOWSIGMA 0.
#define CHARGEPARLIMITLOWMU 0. 
#define CHARGEPARLIMITLOWBGSCALE 0.
#define CHARGEPARLIMITLOWBGASYN 0.
#define CHARGEPARLIMITHIGHSCALE 10.
#define CHARGEPARLIMITHIGHSIGMA 10.
#define CHARGEPARLIMITHIGHMU 10.
#define CHARGEPARLIMITHIGHBGSCALE 10.
#define CHARGEPARLIMITHIGHBGASYN 10.
#define CHARGEQCROSS 0.06
#define CHARGEFITLEGENDXLOW 0.40
#define CHARGEFITLEGENDYLOW 0.74
#define CHARGEFITLEGENDXHIGH 0.89
#define CHARGEFITLEGENDYHIGH 0.89
#define CHARGEFITSTATSXLOW  0.50
#define CHARGEFITSTATSYLOW 0.43
#define CHARGEFITSTATSXHIGH 0.89
#define CHARGEFITSTATSYHIGH 0.73
#define CHARGEFITLEGENDEXCLUSIONXLOW 0.58
#define CHARGEFITLEGENDEXCLUSIONYLOW 0.36
#define CHARGEFITLEGENDEXCLUSIONXHIGH 0.89
#define CHARGEFITLEGENDEXCLUSIONYHIGH 0.39

#define TIMINGNBINS 100
#define TIMINGMIN 0.
#define TIMINGMAX 50.
#define TIMINGLEGENDXLOW 0.58
#define TIMINGLEGENDYLOW 0.78
#define TIMINGLEGENDXHIGH 0.89
#define TIMINGLEGENDYHIGH 0.89

#define CANVASSIZE 1000.
#define LEGENDLINECOLOR 0
#define LEGENDFILLCOLOR 0

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TH1F **allocateH1Array(const unsigned int nH1,
		       const char *name,
		       const char *title,
		       const unsigned int nBinsX,
		       const double xMin,
		       const double xMax,
		       const char *xTitle,
		       const unsigned int color){

  cout << " - allocating H1 array " << endl;
  
  TH1F **h1 = new TH1F*[nH1];
  for(unsigned int iH1=0; iH1<nH1; iH1++){
    char nameStr[100];
    sprintf(nameStr, name, iH1);
    char titleStr[100];
    sprintf(titleStr, title, iH1);
    h1[iH1] = new TH1F(nameStr, titleStr,
		       nBinsX, xMin, xMax);
    h1[iH1] -> SetLineColor(1);
    h1[iH1] -> GetXaxis() -> SetTitle(xTitle);
    h1[iH1] -> SetFillColor(color);
    h1[iH1] -> SetDirectory(0);
  }

  return h1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TH2F **allocateH2Array(const unsigned int nH2,
		       const char *name,
		       const char *title,
		       const unsigned int nBinsX,
		       const double xMin,
		       const double xMax,
		       const char *xTitle,
		       const unsigned int nBinsY,
		       const double yMin,
		       const double yMax,
		       const char *yTitle){

  cout << " - allocating H2 array " << endl;
  
  TH2F **h2 = new TH2F*[nH2];
  for(unsigned int iH2=0; iH2<nH2; iH2++){
    char nameStr[100];
    sprintf(nameStr, name, iH2);
    char titleStr[100];
    sprintf(titleStr, title, iH2);
    h2[iH2] = new TH2F(nameStr, titleStr,
		       nBinsX, xMin, xMax,
		       nBinsY, yMin, yMax);
    h2[iH2] -> GetXaxis() -> SetTitle(xTitle);
    h2[iH2] -> GetYaxis() -> SetTitle(yTitle);
    h2[iH2] -> SetDirectory(0);
  }

  return h2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TF1 **allocateFArray(const unsigned int nF,
		     const char *expr,
		     const double rangeLow,
		     const double rangeHigh){

  cout << " - allocating F array " << endl;

  TF1 **f = new TF1*[nF];
  for(unsigned int iF=0; iF<nF; iF++){
    char name[100];
    sprintf(name, "channel %d", iF);
    f[iF] = new TF1(name, expr, rangeLow, rangeHigh);
  }

  return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TCanvas **allocateCanvasArray(const unsigned int nCC, 
			      const char *flag,
			      const double posX,
			      const double posY,
			      const double widthX,
			      const double widthY,
			      const bool logX,
			      const bool logY,
			      const bool logZ){

  cout << " - allocating Canvas array " << endl;

  TCanvas **cc = new TCanvas*[nCC];
  for(unsigned int iCC=0; iCC<nCC; iCC++){
    char name[100];
    sprintf(name, flag, iCC);
    cc[iCC] = new TCanvas(name, name,
			  posX, posY, widthX, widthY);
    cc[iCC] -> SetLogx(logX);
    cc[iCC] -> SetLogy(logY);
    cc[iCC] -> SetLogz(logZ);
  }

  return cc;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void deleteFArray(const unsigned int nF, 
		  TF1 **f){

  cout << " - deleting F array " << endl;

  for(unsigned int iF=0; iF<nF; iF++){
    delete f[iF];
  }
  delete[] f;
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void deleteH1Array(const unsigned int nH1, 
		   TH1F **h1){

  cout << " - deleting H1 array " << endl;

  for(unsigned int iH1=0; iH1<nH1; iH1++){
    delete h1[iH1];
  }
  delete[] h1;
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void deleteH2Array(const unsigned int nH2, 
		   TH2F **h2){

  cout << " - deleting H2 array " << endl;

  for(unsigned int iH2=0; iH2<nH2; iH2++){
    delete h2[iH2];
  }
  delete[] h2;
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void deleteCanvasArray(const unsigned int nCC, 
		       TCanvas **cc){

  cout << " - deleting Canvas array " << endl;

  for(unsigned int iCC=0; iCC<nCC; iCC++){
    delete cc[iCC];
  }
  delete[] cc;
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void setTimingCuts(double *cut_Timing){ 

  cout << " - setting timing cuts " << endl;

  cut_Timing[0] = 15.;
  cut_Timing[1] = 10.;
  cut_Timing[2] = 15.;

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void setFunctionParameters(const unsigned int nHits, 
			   const unsigned int nPar, 
			   TF1 **f,
			   double *par, 
			   string *parName, 
			   double *parLimitLow, 
			   double *parLimitHigh){

  cout << " - setting function parameters " << endl;

  for(unsigned int iHit=0; iHit<nHits; iHit++){
    for(unsigned int iPar=0; iPar<nPar; iPar++){
      f[iHit] -> SetParameter(iPar, par[iPar]);
      f[iHit] -> SetParName(iPar, parName[iPar].c_str());
      f[iHit] -> SetParLimits(iPar, parLimitLow[iPar], parLimitHigh[iPar]);
    }
  }

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TF1 **optimizeCutT0(const unsigned int nHits,
		    TH1F **h1,
		    double *cutLow,
		    double *cutHigh,
		    double *nEventsTotal,
		    double *nEventsTotalErr,
		    double *nEventsBg,
		    double *nEventsBgErr,
		    double *nEventsSignal,
		    double *nEventsSignalErr,
		    double *bgRejection,
		    double *bgRejectionErr){

  cout << " - optimizing T0 cut " << endl;

  TCanvas *cTmp = new TCanvas();
  
  TF1 **f = allocateFArray(nHits,
			   "([0] / 2.) * ( (1. + TMath::Erf( (x-[1]) / (sqrt(2.) * [2]) )) * (1. + TMath::Erf( -(x-[3]) / (sqrt(2.) * [4]) )) ) + [5]",
			   T0FITRANGELOW,
			   T0FITRANGEHIGH);

  double par[T0NPAR];
  par[0] = T0PARINITPLATEAU;
  par[1] = T0PARINITMULEFT;
  par[2] = T0PARINITSIGMALEFT;
  par[3] = T0PARINITMURIGHT;
  par[4] = T0PARINITSIGMARIGHT;
  par[5] = T0PARINITOFFSET;
  string parName[T0NPAR];
  parName[0] = "plateau #left[ns^{-1}#right]";
  parName[1] = "#mu_{L} [ns]";
  parName[2] = "#sigma_{L} [ns]";
  parName[3] = "#mu_{R} [ns]";
  parName[4] = "#sigma_{R} [ns]";
  parName[5] = "offset #left[ns^{-1}#right]";
  double parLimitLow[T0NPAR];
  parLimitLow[0] = T0PARLIMITLOWPLATEAU;
  parLimitLow[1] = T0PARLIMITLOWMULEFT;
  parLimitLow[2] = T0PARLIMITLOWSIGMALEFT;
  parLimitLow[3] = T0PARLIMITLOWMURIGHT;
  parLimitLow[4] = T0PARLIMITLOWSIGMARIGHT;
  parLimitLow[5] = T0PARLIMITLOWOFFSET;
  double parLimitHigh[T0NPAR];
  parLimitHigh[0] = T0PARLIMITHIGHPLATEAU;
  parLimitHigh[1] = T0PARLIMITHIGHMULEFT;
  parLimitHigh[2] = T0PARLIMITHIGHSIGMALEFT;
  parLimitHigh[3] = T0PARLIMITHIGHMURIGHT;
  parLimitHigh[4] = T0PARLIMITHIGHSIGMARIGHT;
  parLimitHigh[5] = T0PARLIMITHIGHOFFSET;

  setFunctionParameters(nHits, T0NPAR, f,
			par, parName, parLimitLow, parLimitHigh);

  for(unsigned int iHit=0; iHit<nHits; iHit++){
    h1[iHit] -> Fit(f[iHit], "R && Q");
    cutLow[iHit] = f[iHit] -> GetParameter(1) - T0CUTNSIGMA * f[iHit] -> GetParameter(2);
    cutHigh[iHit] = f[iHit] -> GetParameter(3) + T0CUTNSIGMA * f[iHit] -> GetParameter(4);
    nEventsTotal[iHit] = f[iHit] -> Integral(cutLow[iHit], cutHigh[iHit]);
    nEventsTotalErr[iHit] = f[iHit] -> IntegralError(cutLow[iHit], cutHigh[iHit]);
    double rangeSize = cutHigh[iHit] - cutLow[iHit];
    nEventsBg[iHit] = rangeSize * f[iHit] -> GetParameter(5);
    nEventsBgErr[iHit] = rangeSize * f[iHit] -> GetParError(5);
    nEventsSignal[iHit] = nEventsTotal[iHit] - nEventsBg[iHit];
    nEventsSignalErr[iHit] = sqrt(nEventsTotalErr[iHit] * nEventsTotalErr[iHit] - nEventsBgErr[iHit] * nEventsBgErr[iHit]);
    bgRejection[iHit] = 1. - nEventsBg[iHit] / nEventsSignal[iHit];
    bgRejectionErr[iHit] = sqrt(nEventsSignal[iHit] * nEventsSignal[iHit] * nEventsBgErr[iHit] * nEventsBgErr[iHit] + nEventsBg[iHit] * nEventsBg[iHit] * nEventsSignalErr[iHit] * nEventsSignalErr[iHit]) / (nEventsSignal[iHit] * nEventsSignal[iHit]);
  }

  delete cTmp;
  return f;  
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void fitCharge(TH1F *h1, 
	       TF1 *fTotal,
	       double &cutChargeSharing_Charge,
	       double &cutNEventsBg_Charge,
	       double &cutNEventsBgErr_Charge,
	       double &cutNEventsSignal_Charge,
	       double &cutNEventsSignalErr_Charge,
	       double &cutCSRejection_Charge,
	       double &cutCSRejectionErr_Charge){
  
  cout << " - fitting charge distribution " << endl;

  // temporary fit of the BG only
  TF1 *fBg = new TF1("fBg", "[0]/(x-[1])", CHARGEFITRANGELOW, CHARGEQCROSS);
  h1 -> Fit(fBg, "R && Q");
  double p0 = fBg -> GetParameter(0);
  double p1 = fBg -> GetParameter(1);
  delete fBg;

  //  temporary fit of the signal only
  TF1 *fSignal = new TF1("fSignal", "( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) )", CHARGEQCROSS, CHARGEFITRANGEHIGH);
  fSignal -> SetParameter(0, 2.);
  fSignal -> SetParameter(1, 0.05);
  fSignal -> SetParameter(2, 0.1);
  h1 -> Fit(fSignal, "R && Q");

  // final fit BG+signal
  fTotal -> SetParameter(0, fSignal -> GetParameter(0));
  fTotal -> SetParameter(1, fSignal -> GetParameter(1));
  fTotal -> SetParameter(2, fSignal -> GetParameter(2));
  fTotal -> SetParameter(3, p0);
  fTotal -> SetParameter(4, p1);
  h1 -> Fit(fTotal, "R && Q");

  // finding BG-signal intersection point
  TF1 *fIntersect = new TF1("fIntersect",
			    "TMath::Abs(( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) ) - ([3] / (x-[4])) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[1]))) / 2. )",
			    0., fTotal -> GetParameter(2));
  fIntersect -> SetParameter(0, fTotal -> GetParameter(0));
  fIntersect -> SetParameter(1, fTotal -> GetParameter(1));
  fIntersect -> SetParameter(2, fTotal -> GetParameter(2));
  fIntersect -> SetParameter(3, fTotal -> GetParameter(3));
  fIntersect -> SetParameter(4, fTotal -> GetParameter(4));
  cutChargeSharing_Charge = fIntersect -> GetMinimumX();
  delete fIntersect;

  // this part does not work yet
  /* // computing rejection power */
  /* fSignal -> SetParameter(0, fTotal -> GetParameter(0)); */
  /* fSignal -> SetParameter(1, fTotal -> GetParameter(1)); */
  /* fSignal -> SetParameter(2, fTotal -> GetParameter(2)); */
  /* fSignal -> SetParError(0, fTotal -> GetParError(0)); */
  /* fSignal -> SetParError(1, fTotal -> GetParError(1)); */
  /* fSignal -> SetParError(2, fTotal -> GetParError(2)); */
  /* fSignal -> SetRange(CHARGEMIN, CHARGEMAX); */
  /* cutNEventsSignal_Charge = fSignal -> Integral(cutChargeSharing_Charge, CHARGEMAX); */
  /* //  cutNEventsSignalErr_Charge = fSignal -> IntegralError(cutChargeSharing_Charge, 1000.); */
  /* fBg = new TF1("fBg", "[0]/(x-[1]) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[3]))) / 2. ", CHARGEMIN, CHARGEMAX); */
  /* fBg -> SetParameter(0, fTotal -> GetParameter(3)); */
  /* fBg -> SetParameter(1, fTotal -> GetParameter(4)); */
  /* fBg -> SetParameter(2, fTotal -> GetParameter(2)); */
  /* fBg -> SetParameter(3, fTotal -> GetParameter(1)); */
  /* fBg -> SetParError(0, fTotal -> GetParError(3)); */
  /* fBg -> SetParError(1, fTotal -> GetParError(4)); */
  /* fBg -> SetParError(2, fTotal -> GetParError(2)); */
  /* fBg -> SetParError(3, fTotal -> GetParError(1)); */
  /* cutNEventsBg_Charge = fBg -> Integral(cutChargeSharing_Charge, CHARGEMAX); */
  /* //  cutNEventsBgErr_Charge = fBg -> IntegralError(cutChargeSharing_Charge, 1000.); */

  /* delete fBg; */

  delete fSignal;

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TF1 **optimizeChargeCut(const unsigned int nHits,
			TH1F **h1_bg,
			TH2F **h2,
			TH2F **h2_bgSubtracted,
			TH1F **h1_bgSubtracted,
			double *cutLow_T0,
			double *cutHigh_T0,
			double *cutChargeSharing_Charge,
			double *cutNEventsBg_Charge,
			double *cutNEventsBgErr_Charge,
			double *cutNEventsSignal_Charge,
			double *cutNEventsSignalErr_Charge,
			double *cutCSRejection_Charge,
			double *cutCSRejectionErr_Charge){

  cout << " - optimizing charge cut " << endl;

  TCanvas *cTmp = new TCanvas();

  TF1 **f = allocateFArray(nHits,
			   "( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) ) + ([3] / (x-[4])) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[1]))) / 2. ",
			   CHARGEFITRANGELOW,
			   CHARGEFITRANGEHIGH);

  double par[CHARGENPAR];
  par[0] = CHARGEPARINITSCALE;
  par[1] = CHARGEPARINITSIGMA;
  par[2] = CHARGEPARINITMU;
  par[3] = CHARGEPARINITBGSCALE;
  par[4] = CHARGEPARINITBGASYN;
  string parName[CHARGENPAR];
  parName[0] = "scale #left[mV^{-1}#right]";
  parName[1] = "#sigma [mV]";
  parName[2] = "#mu [mV]";
  parName[3] = "CS scale [mV]";
  parName[4] = "CS asyn. #left[mV^{-1}#right]";
  double parLimitLow[CHARGENPAR];
  parLimitLow[0] = CHARGEPARLIMITLOWSCALE;
  parLimitLow[1] = CHARGEPARLIMITLOWSIGMA;
  parLimitLow[2] = CHARGEPARLIMITLOWMU;
  parLimitLow[3] = CHARGEPARLIMITLOWBGSCALE;
  parLimitLow[4] = CHARGEPARLIMITLOWBGASYN;
  double parLimitHigh[CHARGENPAR];
  parLimitHigh[0] = CHARGEPARLIMITHIGHSCALE;
  parLimitHigh[1] = CHARGEPARLIMITHIGHSIGMA;
  parLimitHigh[2] = CHARGEPARLIMITHIGHMU;
  parLimitHigh[3] = CHARGEPARLIMITHIGHBGSCALE;
  parLimitHigh[4] = CHARGEPARLIMITHIGHBGASYN;

  setFunctionParameters(nHits, CHARGENPAR, f,
			par, parName, parLimitLow, parLimitHigh);

  for(unsigned int iHit=0; iHit<nHits; iHit++){

    // calculating number of bins outside CHARGE cut region
    double nBinsT0 = T0NBINS * (cutLow_T0[iHit] - T0MIN + T0MAX - cutHigh_T0[iHit]) / (T0MAX - T0MIN);
    cout << "- channel " << iHit << ": "
	 << "nBinsT0 = " << nBinsT0
	 << endl;

    // subtracting background from 2D histogram
    TH2F *h2tmp = new TH2F("h2", "h2", T0NBINS, T0MIN, T0MAX, CHARGENBINS, CHARGEMIN, CHARGEMAX);
    for(unsigned int iY=0; iY<CHARGENBINS; iY++){
      double charge = h1_bg[iHit] -> GetBinContent(iY+1) / nBinsT0;
      for(unsigned int iX=0; iX<T0NBINS; iX++){
	h2tmp -> SetBinContent(iX+1, iY+1, charge);
	h2_bgSubtracted[iHit] -> SetBinContent(iX+1, iY+1, h2[iHit] -> GetBinContent(iX+1, iY+1));
      }
    }
    h2tmp -> Sumw2();
    h2_bgSubtracted[iHit] -> Sumw2();
    h2_bgSubtracted[iHit] -> Add(h2tmp, -1.);
    delete h2tmp;

    // calculating 1D signal histogram
    TH1D *h1Tmp = h2_bgSubtracted[iHit] -> ProjectionY("tmp", 
						       h2_bgSubtracted[iHit] -> GetXaxis() -> FindBin(cutLow_T0[iHit]),
						       h2_bgSubtracted[iHit] -> GetXaxis() -> FindBin(cutHigh_T0[iHit]),
						       "e");
    for(unsigned int iY=0; iY<CHARGENBINS; iY++){
      h1_bgSubtracted[iHit] -> SetBinContent(iY, h1Tmp -> GetBinContent(iY));
    }
    delete h1Tmp;

    // fitting signal
    fitCharge(h1_bgSubtracted[iHit], f[iHit],
	      cutChargeSharing_Charge[iHit],
	      cutNEventsBg_Charge[iHit],
	      cutNEventsBgErr_Charge[iHit],
	      cutNEventsSignal_Charge[iHit],
	      cutNEventsSignalErr_Charge[iHit],
	      cutCSRejection_Charge[iHit],
	      cutCSRejectionErr_Charge[iHit]);

  }

  delete cTmp;
  return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void fillExclusionPlot(TH1F *h1,
		       const double rangeLow,
		       const double rangeHigh,
		       const double val){

  cout << " - filling exclusion plot " << endl;

  for(int iBin=1; iBin<h1 -> GetNbinsX()+2; iBin++){
    if(h1 -> GetXaxis() -> GetBinLowEdge(iBin) >= rangeLow && h1 -> GetXaxis() -> GetBinLowEdge(iBin) <= rangeHigh) h1 -> SetBinContent(iBin+1, val);
  }
  h1 -> SetLineColor(1);
  h1 -> SetFillColor(1);
  h1 -> SetFillStyle(3004);

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotT0(TH1F *h1,
	    TH1F *h1_cut,
	    TF1 *f,
	    TLegend *leg, 
	    TPaveStats *ptstats,
	    TPaveStats *ptresults,
	    TH1F *h1_exclusionLeft,
	    TH1F *h1_exclusionRight,
	    const double cutTiming,
	    const double cutLow,
	    const double cutHigh,
	    const double cutBgRejection,
	    const double cutBgRejectionErr,
	    const unsigned int nPar){

  cout << " - plotting T0 " << endl;

  double yScaleFactor = 1.1;

  fillExclusionPlot(h1_exclusionLeft,
		    h1 -> GetXaxis() -> GetBinLowEdge(1), cutLow,
		    yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));
  fillExclusionPlot(h1_exclusionRight,
		    cutHigh, h1 -> GetBinLowEdge(h1 -> GetNbinsX()),
		    yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));

  h1 -> Draw();
  h1_cut -> Draw("same");
  h1_exclusionLeft -> Draw("same");
  h1_exclusionRight -> Draw("same");
  h1 -> GetYaxis() -> SetRangeUser(0, yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));

  leg = new TLegend(T0LEGENDXLOW, T0LEGENDYLOW, T0LEGENDXHIGH, T0LEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(1);
  char text[200];
  leg -> AddEntry(h1, "no cut", "f");
  sprintf(text, "timing > %.1lf ns", cutTiming);
  leg -> AddEntry(h1_cut, text, "f");
  leg -> AddEntry(h1_exclusionLeft, "cut region", "f");
  leg -> Draw();

  ptstats = new TPaveStats(T0STATSXLOW, T0STATSYLOW, T0STATSXHIGH, T0STATSYHIGH, "brNDC");
  ptstats -> SetName("stats");
  ptstats -> SetBorderSize(1);
  ptstats -> SetFillColor(LEGENDFILLCOLOR);
  ptstats -> SetLineColor(1);
  ptstats -> SetTextAlign(12);
  ptstats -> SetTextFont(42);
  sprintf(text, "#chi^{2} / ndf = %.2lf / %d", f -> GetChisquare(), f -> GetNDF());
  ptstats -> AddText(text);
  for(unsigned int iPar=0; iPar<nPar; iPar++){
    sprintf(text, "%s  = (%.1lf #pm %.1lf)", f -> GetParName(iPar), f -> GetParameter(iPar), f -> GetParError(iPar));
    ptstats -> AddText(text);
  }
  ptstats -> SetOptStat(0);
  ptstats -> SetOptFit(111);
  ptstats -> Draw();

  ptresults = new TPaveStats(T0RESULTSXLOW, T0RESULTSYLOW, T0RESULTSXHIGH, T0RESULTSYHIGH, "brNDC");
  ptresults -> SetTextSize(0.025);
  ptresults -> SetName("results");
  ptresults -> SetBorderSize(1);
  ptresults -> SetFillColor(LEGENDFILLCOLOR);
  ptresults -> SetLineColor(1);
  ptresults -> SetTextAlign(12);
  ptresults -> SetTextFont(42);
  sprintf(text, "T_{0} cut #in [%.1lf #pm %.1lf] ns", cutLow, cutHigh);
  ptresults -> AddText(text);
  sprintf(text, "bg rej. = (%.1lf #pm %.1lf) %%", 100. * cutBgRejection, 100. * cutBgRejectionErr);
  ptresults -> AddText(text);
  ptresults -> SetOptStat(0);
  ptresults -> SetOptFit(111);
  ptresults -> Draw();
  
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotCharge(TH1F *h1, 
		TH1F *h1_cut,
		TLegend *leg,
		const double cutLowT0,
		const double cutHighT0){

  cout << " - plotting charge " << endl;

  h1 -> Draw();
  h1_cut -> Draw("same");

  leg = new TLegend(T0LEGENDXLOW, T0LEGENDYLOW, T0LEGENDXHIGH, T0LEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(LEGENDLINECOLOR);
  char text[200];
  leg -> AddEntry(h1, "no cut", "f");
  sprintf(text, "T_{0} #notin [%.1lf:%.1lf] ns", cutLowT0, cutHighT0);
  leg -> AddEntry(h1_cut, text, "f");
  leg -> Draw();

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotChargeFit(TH1F *h1, 
		   TF1 *f,
		   TF1 *fBg,
		   TF1 *fSignal,
		   TH1F *h1_exclusion,
		   const double cut,
		   TLegend *leg, 
		   TPaveStats *ptstats,
		   TLegend *legExclusion,
		   const unsigned int nPar){

  cout << " - plotting charge fit " << endl;

  const double yScaleFactor = 1.1;

  gStyle -> SetOptFit(0);

  h1 -> GetYaxis() -> SetRangeUser(0., yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));
  h1 -> Draw("e");

  fillExclusionPlot(h1_exclusion,
                    h1 -> GetXaxis() -> GetBinLowEdge(1), cut,
                    yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));
  h1_exclusion -> Draw("same");

  fSignal -> SetParameter(0, f -> GetParameter(0));
  fSignal -> SetParameter(1, f -> GetParameter(1));
  fSignal -> SetParameter(2, f -> GetParameter(2));
  fSignal -> SetLineStyle(7);
  fSignal -> Draw("same");

  fBg -> SetParameter(0, f -> GetParameter(3));
  fBg -> SetParameter(1, f -> GetParameter(4));
  fBg -> SetParameter(2, f -> GetParameter(2));
  fBg -> SetParameter(3, f -> GetParameter(1));
  fBg -> SetLineStyle(8);
  fBg -> Draw("same");

  leg = new TLegend(CHARGEFITLEGENDXLOW, CHARGEFITLEGENDYLOW, CHARGEFITLEGENDXHIGH, CHARGEFITLEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(LEGENDLINECOLOR);
  leg -> AddEntry(h1, "data, bg subtracted", "f");
  leg -> AddEntry(f, "fit to data");
  leg -> AddEntry(fBg, "charge sharing component");
  leg -> AddEntry(fSignal, "signal component");
  leg -> Draw();

  leg = new TLegend(CHARGEFITLEGENDEXCLUSIONXLOW, CHARGEFITLEGENDEXCLUSIONYLOW, CHARGEFITLEGENDEXCLUSIONXHIGH, CHARGEFITLEGENDEXCLUSIONYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(LEGENDLINECOLOR);
  char text[200];
  sprintf(text, "Charge < %.3lf mV", cut);
  leg -> AddEntry(h1_exclusion, text, "f");
  leg -> Draw();

  ptstats = new TPaveStats(CHARGEFITSTATSXLOW, CHARGEFITSTATSYLOW, CHARGEFITSTATSXHIGH, CHARGEFITSTATSYHIGH, "brNDC");
  ptstats -> SetName("stats");
  ptstats -> SetBorderSize(1);
  ptstats -> SetFillColor(LEGENDFILLCOLOR);
  ptstats -> SetLineColor(LEGENDLINECOLOR);
  ptstats -> SetTextAlign(12);
  ptstats -> SetTextFont(42);
  sprintf(text, "#chi^{2} / ndf = %.2lf / %d", f -> GetChisquare(), f -> GetNDF());
  ptstats -> AddText(text);
  for(unsigned int iPar=0; iPar<nPar; iPar++){
    sprintf(text, "%s  = (%lf #pm %lf)", f -> GetParName(iPar), f -> GetParameter(iPar), f -> GetParError(iPar));
    ptstats -> AddText(text);
  }
  ptstats -> SetOptStat(0);
  ptstats -> SetOptFit(111);
  ptstats -> Draw();

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotTiming(TH1F *h1, 
		TH1F *h1Cut,
		TLegend *leg){

  cout << " - plotting timing" << endl;

  h1 -> Draw();
  h1Cut -> Draw("same");

  leg = new TLegend(TIMINGLEGENDXLOW, TIMINGLEGENDYLOW, TIMINGLEGENDXHIGH, TIMINGLEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(LEGENDLINECOLOR);
  char text[200];
  leg -> AddEntry(h1, "no cut", "f");
  leg -> AddEntry(h1Cut, "T_{0} and charge cuts", "f");
  leg -> Draw();

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void draw(const unsigned int nH2,
	  TH2F **h2_Timing_vs_T0,
	  TH2F **h2_Charge_vs_T0,
	  TH2F **h2_Charge_vs_Timing,
	  TH1F **h1_T0,
	  TH1F **h1_T0_TimingCut,
	  TH1F **h1_exclusionLeft_T0,
	  TH1F **h1_exclusionRight_T0,
	  double *cut_Timing,
	  const unsigned int T0Npar,
	  TF1 **f_T0,
	  double *cutLow_T0,
	  double *cutHigh_T0,
	  double *cutBgRejection_T0,
	  double *cutBgRejectionErr_T0,
	  TH1F **h1_Charge,
	  TH1F **h1_Charge_T0Cut,
	  TH2F **h2_Charge_vs_T0_bgSubtracted,
	  TH1F **h1_Charge_bgSubtracted,
	  TF1 **f_Charge,
	  double *cutChargeSharing_Charge,
	  TH1F **h1_exclusion_Charge,
	  TH1F **h1_Timing,
	  TH1F **h1_Timing_T0Cut_ChargeCut,
	  const string outputFolder,
	  const string flag,
	  const bool holdPlots){

  cout << " - drawing " << endl;

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);

  // allocating canvases
  TCanvas **cc_Timing_vs_T0 = allocateCanvasArray(nH2, "cc_Timing_vs_T0_%d",
  						  0, 0, CANVASSIZE, CANVASSIZE,
  						  false, true, true);

  TCanvas **cc_Charge_vs_T0 = allocateCanvasArray(nH2, "cc_Charge_vs_T0_%d",
						  CANVASSIZE, 0, CANVASSIZE, CANVASSIZE,
						  false, true, true);

  TCanvas **cc_Charge_vs_Timing = allocateCanvasArray(nH2, "cc_Charge_vs_Timing_%d",
  						      2.*CANVASSIZE, 0, CANVASSIZE, CANVASSIZE,
  						      false, true, true);

  TCanvas **cc_T0 = allocateCanvasArray(nH2, "cc_T0_%d",
  					0, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					false, false, false);

  TCanvas **cc_Charge = allocateCanvasArray(nH2, "cc_Charge_%d",
  					    CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					    true, true, false);

  TCanvas **cc_Charge_vs_T0_bgSubtracted = allocateCanvasArray(nH2, "cc_Charge_vs_T0_bgSubtracted_%d",
							       2.*CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
							       false, true, false);

  TCanvas **cc_Charge_bgSubtracted = allocateCanvasArray(nH2, "cc_Charge_bgSubtracted_%d",
							 2.*CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
							 false, false, false);
  TCanvas **cc_Timing = allocateCanvasArray(nH2, "cc_Timing_%d",
  					    CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					    false, false, false);

  // allocating functions
  TF1 **f_ChargeBg = allocateFArray(nH2,
				    "([0] / (x-[1])) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[3]))) / 2.",
				    CHARGEFITRANGELOW,
				    CHARGEFITRANGEHIGH);

  TF1 **f_ChargeSignal = allocateFArray(nH2,
					"( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) )",				       
					CHARGEFITRANGELOW,
					CHARGEFITRANGEHIGH);
  
  // allocating TLegends and TPaves
  TLegend **leg_T0 = new TLegend*[nH2];
  TPaveStats **ptstats_T0 = new TPaveStats*[nH2];
  TPaveStats **ptresults_T0 = new TPaveStats*[nH2];
  TLegend **leg_Charge = new TLegend*[nH2];
  TLegend **leg_ChargeFit = new TLegend*[nH2];
  TPaveStats **ptstats_ChargeFit = new TPaveStats*[nH2];
  TLegend **legExclusion_ChargeFit = new TLegend*[nH2];
  TLegend **leg_Timing = new TLegend*[nH2];

  // drawing
  for(unsigned int iH=0; iH<nH2; iH++){

    stringstream channel;
    channel << iH;

    cc_Timing_vs_T0[iH] -> cd();
    h2_Timing_vs_T0[iH] -> Draw("colz");
    cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.png").c_str());
    cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.pdf").c_str());
    cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.eps").c_str());
    cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.root").c_str());

    cc_Charge_vs_T0[iH] -> cd();
    h2_Charge_vs_T0[iH] -> Draw("colz");
    cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.png").c_str());
    cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.pdf").c_str());
    cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.eps").c_str());
    cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.root").c_str());

    cc_Charge_vs_Timing[iH] -> cd();
    h2_Charge_vs_Timing[iH] -> Draw("colz");
    cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.png").c_str());
    cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.pdf").c_str());
    cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.eps").c_str());
    cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.root").c_str());

    cc_T0[iH] -> cd();
    plotT0(h1_T0[iH], h1_T0_TimingCut[iH], f_T0[iH],
    	   leg_T0[iH], ptstats_T0[iH], ptresults_T0[iH],
    	   h1_exclusionLeft_T0[iH], h1_exclusionRight_T0[iH],
    	   cut_Timing[iH], cutLow_T0[iH], cutHigh_T0[iH],
    	   cutBgRejection_T0[iH], cutBgRejectionErr_T0[iH],
    	   T0Npar);
    cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.png").c_str());
    cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.pdf").c_str());
    cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.eps").c_str());
    cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.root").c_str());

    cc_Charge[iH] -> cd();
    plotCharge(h1_Charge[iH], h1_Charge_T0Cut[iH],
    	       leg_Charge[iH],
    	       cutLow_T0[iH], cutHigh_T0[iH]);
    cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.png").c_str());
    cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.pdf").c_str());
    cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.eps").c_str());
    cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.root").c_str());

    cc_Charge_vs_T0_bgSubtracted[iH] -> cd();
    h2_Charge_vs_T0_bgSubtracted[iH] -> SetMinimum(0);
    h2_Charge_vs_T0_bgSubtracted[iH] -> Draw("colz");
    cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.png").c_str());
    cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.pdf").c_str());
    cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.eps").c_str());
    cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.root").c_str());

    cc_Charge_bgSubtracted[iH] -> cd();
    plotChargeFit(h1_Charge_bgSubtracted[iH], 
		  f_Charge[iH], f_ChargeBg[iH], f_ChargeSignal[iH],
		  h1_exclusion_Charge[iH], cutChargeSharing_Charge[iH],
		  leg_ChargeFit[iH], ptstats_ChargeFit[iH], legExclusion_ChargeFit[iH],
		  CHARGENPAR);
    cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.png").c_str());
    cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.pdf").c_str());
    cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.eps").c_str());
    cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.root").c_str());

    cc_Timing[iH] -> cd();
    plotTiming(h1_Timing[iH], h1_Timing_T0Cut_ChargeCut[iH],
	       leg_Timing[iH]);
    cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.png").c_str());
    cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.pdf").c_str());
    cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.eps").c_str());
    cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.root").c_str());
  }

  // cleaning memory
  if(holdPlots) return ;
  deleteFArray(MAX_HITS, f_ChargeBg);
  deleteFArray(MAX_HITS, f_ChargeSignal);
  deleteCanvasArray(nH2, cc_Timing_vs_T0);
  deleteCanvasArray(nH2, cc_Charge_vs_T0);
  deleteCanvasArray(nH2, cc_Charge_vs_Timing);
  deleteCanvasArray(nH2, cc_T0);
  deleteCanvasArray(nH2, cc_Charge);
  deleteCanvasArray(nH2, cc_Charge_vs_T0_bgSubtracted);
  deleteCanvasArray(nH2, cc_Timing);
  return ;
}

#endif
