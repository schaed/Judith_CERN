#ifndef OPTIMIZECUTS_HH
#define OPTIMIZECUTS_HH

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
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
#include <TGraph.h>
#include <TVirtualFFT.h>

#define BC 25. // bunch crossing time, BC = 25 ns
#define MAX_HITS 3 // number of channels stored in the raw data file
#define CHI2LIMITLOW 0.1 // reduced chi2 lower limit for considering a fit as good
#define CHI2LIMITHIGH 10 // reduced chi2 higher limit for considering a fit as good
#define CHARGEFFTFREQCONVERSION 0.008 // constant parameter used to determine the optimal frequency cut for FFT filtering the charge distributions
#define HIGHSTATISTICSLIMIT 10000 // minimum number of events for considering the run as a high statistics one

#define HARM1MODNBINS 200
#define HARM1MODMIN 0
#define HARM1MODMAX 10
#define HARM1PHNBINS 200
#define HARM1PHMIN -3.14159
#define HARM1PHMAX 3.14159

#define T0NBINS 405
#define T0MIN -9
#define T0MAX 800.
#define T0REBINFACTOR 1 // rebinning factor for the clone histogram from which the T0 fit range is estimated (from the position of the maximum of the distribution)
#define T0FITRANGELOW 180.
#define T0FITRANGEHIGH 240.
#define T0NPAR 6
#define T0PARINITPLATEAU 200.
#define T0PARINITMULEFT 200.
#define T0PARINITSIGMALEFT 2.
#define T0PARINITDELTAMU 25.
#define T0PARINITSIGMARIGHT 2.
#define T0PARINITOFFSET 10.
#define T0PARLIMITLOWPLATEAU 100.
#define T0PARLIMITLOWMULEFT 190.
#define T0PARLIMITLOWSIGMALEFT 1.
#define T0PARLIMITLOWDELTAMU 20
#define T0PARLIMITLOWSIGMARIGHT 1.
#define T0PARLIMITLOWOFFSET 0.
#define T0PARLIMITHIGHPLATEAU 300.
#define T0PARLIMITHIGHMULEFT 210.
#define T0PARLIMITHIGHSIGMALEFT 3.
#define T0PARLIMITHIGHDELTAMU 30.
#define T0PARLIMITHIGHSIGMARIGHT 3.
#define T0PARLIMITHIGHOFFSET 100.
#define T0LEGENDXLOW 0.58
#define T0LEGENDYLOW 0.78
#define T0LEGENDXHIGH 0.89
#define T0LEGENDYHIGH 0.89
#define T0STATSXLOW  0.60
#define T0STATSYLOW 0.44
#define T0STATSXHIGH 0.89
#define T0STATSYHIGH 0.76
#define T0RESULTSXLOW  0.60
#define T0RESULTSYLOW 0.32
#define T0RESULTSXHIGH 0.89
#define T0RESULTSYHIGH 0.42
#define T0CUTNSIGMA 3. // number of sigmas for T0 cut

#define CHARGENBINS 500
#define CHARGEMIN 0.
#define CHARGEMAX 0.5
#define CHARGECUTTYPE 0 // 0: from crossing point between charge sharing curve and signal distribution curve, other types not defined yet
#define CHARGEFITRANGEFINDINGALGORITHM 1 // 0: from rebinned histogram, 1: from FFT filtered histogram
#define CHARGELEGENDXLOW 0.52
#define CHARGELEGENDYLOW 0.76
#define CHARGELEGENDXHIGH 0.89
#define CHARGELEGENDYHIGH 0.89
#define CHARGERESULTSXLOW  0.52
#define CHARGERESULTSYLOW 0.52
#define CHARGERESULTSXHIGH 0.89
#define CHARGERESULTSYHIGH 0.62
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

#define TIMINGNBINS 50
#define TIMINGMIN 0.
#define TIMINGMAX 200.
#define TIMINGMAXPRELIMINARY 50.
#define TIMINGPRELIMINARYCUT 5.
#define TIMINGLEGENDXLOW 0.52
#define TIMINGLEGENDYLOW 0.78
#define TIMINGLEGENDXHIGH 0.89
#define TIMINGLEGENDYHIGH 0.89
#define TIMINGNQUANTILES 1000
#define TIMINGCUTFRACTION 0.01 // fractional area for Timing cut
#define TIMINGRESULTSXLOW  0.52
#define TIMINGRESULTSYLOW 0.70
#define TIMINGRESULTSXHIGH 0.89
#define TIMINGRESULTSYHIGH 0.75
#define TIMINGCUTNSIGMA 3. // number of sigmas for T0 cut

#define CANVASSIZE 1000.
#define LEGENDLINECOLOR 0
#define LEGENDFILLCOLOR 0

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

struct plot_struct {
  TH1F **h1_Harm1Mod;
  TH1F **h1_Harm1Ph;
  TH1F **h1_Harm1Re;
  TH1F **h1_Harm1Im;
  TH2F **h2_Harm1Polar;
  TH2F **h2_Harm1Polar_bgSubtracted;
  TH2F **h2_Harm1Cartesian;
  TH2F **h2_Harm1Cartesian_bgSubtracted;
  TH1F **h1_Harm1Re_bgSubtracted;
  TH1F **h1_Harm1Im_bgSubtracted;
  TH1F **h1_Harm1Mod_bgSubtracted;
  TH1F **h1_Harm1Ph_bgSubtracted;
  TH2F **h2_Harm1Mod_vs_T0;
  TH2F **h2_Timing_vs_T0;
  TH2F **h2_Charge_vs_T0;
  TH2F **h2_Charge_vs_Timing;
  TH1F **h1_T0;
  TH1F **h1_T0_TimingCut;
  TH1F **h1_exclusionLeft_T0;
  TH1F **h1_exclusionRight_T0;
  TH1F **h1_Charge;
  TH1F **h1_Charge_T0Cut_TimingCut;
  TH1F **h1_Charge_T0Cut_bg;
  TH2F **h2_Charge_vs_T0_bgSubtracted;
  TH1F **h1_Charge_bgSubtracted;
  TH1F **h1_Charge_bgSubtracted_filtered;
  TH1F **h1_exclusion_Charge;
  TH1F **h1_Timing;
  TH1F **h1_Timing_T0Cut_ChargeCut;
  TH1F **h1_exclusionLeft_Timing;
  TH1F **h1_exclusionRight_Timing;
} ;

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

const string extractFlag(const char *fileName){
  stringstream name_ss;
  name_ss << fileName;
  string name = name_ss.str();
  int pos1 = name.find("run");
  int pos2 = name.find(".root");
  return (name.substr(pos1, pos2-pos1) + "-cut");
}

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
		       const char *yTitle,
		       const unsigned int color,
		       ofstream &logfile){

  cout << " - allocating H1 array " << endl;
  logfile << "\t- in function allocateH1Array(): " << name << endl;
  
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
    h1[iH1] -> GetYaxis() -> SetTitle(yTitle);
    h1[iH1] -> GetYaxis() -> SetTitleOffset(1.4);
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
		       const char *yTitle,
		       ofstream &logfile){

  cout << " - allocating H2 array " << endl;
  logfile << "\t- in function allocateH2Array(): " << name << endl;
  
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
		     const double rangeHigh,
		     ofstream &logfile){

  cout << " - allocating F array " << endl;
  logfile << "\t- in function allocateFArray(): " << expr << endl;

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
			      const bool logZ,
			      ofstream &logfile){

  cout << " - allocating Canvas array " << endl;
  logfile << "\t- in function allocateCanvasArray(): " << flag << endl;

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
		  TF1 **f,
		  ofstream &logfile){

  cout << " - deleting F array " << endl;
  logfile << "\t- in function deleteFArray()" << endl;

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
		   TH1F **h1,
		   ofstream &logfile){

  cout << " - deleting H1 array " << endl;
  logfile << "\t- in function deleteH1Array()" << endl;

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
		   TH2F **h2,
		   ofstream &logfile){

  cout << " - deleting H2 array " << endl;
  logfile << "\t- in function deleteH2Array()" << endl;

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
		       TCanvas **cc,
		       ofstream &logfile){

  cout << " - deleting Canvas array " << endl;
  logfile << "\t- in function deleteCanvasArray()" << endl;

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

void subtractBgHarm1(const unsigned int nH2,
		     TH2F **h2Polar, 
		     TH2F **h2Polar_bgSubtracted,
		     TH2F **h2Cartesian, 
		     TH2F **h2Cartesian_bgSubtracted,
		     TH1F **h1Re,
		     TH1F **h1Im,
		     TH1F **h1Mod,
		     TH1F **h1Ph,
		     ofstream &logfile){

  cout << " - subtracting background in h1 distributions" << endl;
  logfile << "\t- in function subtractBgHarm1()" << endl;

  for(unsigned int iH2=0; iH2<nH2; iH2++){

    logfile << "\t\t- hit " << iH2 << ": polar coordinates" << endl;
    unsigned int halfRangeArg = h2Polar[iH2] -> GetNbinsX()/2;
    for(unsigned int ix=1; ix<=halfRangeArg; ix++){
      for(int iy=1; iy<=h2Polar[iH2] -> GetNbinsY(); iy++){
	h2Polar_bgSubtracted[iH2] -> SetBinContent(halfRangeArg+ix, iy, h2Polar[iH2] -> GetBinContent(halfRangeArg+ix, iy) - h2Polar[iH2] -> GetBinContent(ix, iy));
      }
    }

    logfile << "\t\t- hit " << iH2 << ": projecting on the phase axis" << endl;
    TH1D *h1TmpPolarX = h2Polar_bgSubtracted[iH2] -> ProjectionX(); 
    for(int ix=1; ix<=h1Ph[iH2] -> GetNbinsX(); ix++){
      h1Ph[iH2] -> SetBinContent(ix, h1TmpPolarX -> GetBinContent(ix));
    }

    logfile << "\t\t- hit " << iH2 << ": projecting on the module axis" << endl;
    TH1D *h1TmpPolarY = h2Polar_bgSubtracted[iH2] -> ProjectionY(); 
    for(int ix=1; ix<=h1Mod[iH2] -> GetNbinsX(); ix++){
      h1Mod[iH2] -> SetBinContent(ix, h1TmpPolarY -> GetBinContent(ix));
    }

    logfile << "\t\t- hit " << iH2 << ": cartesian coordinates" << endl;
    unsigned int halfRangeIm = h2Cartesian[iH2] -> GetNbinsY()/2;
    for(unsigned int iy=1; iy<=halfRangeIm; iy++){
      for(int ix=1; ix<=h2Cartesian[iH2] -> GetNbinsX(); ix++){
	h2Cartesian_bgSubtracted[iH2] -> SetBinContent(ix, halfRangeIm+iy, h2Cartesian[iH2] -> GetBinContent(ix, halfRangeIm+iy) - h2Cartesian[iH2] -> GetBinContent(ix, halfRangeIm-iy));
      }
    }

    logfile << "\t\t- hit " << iH2 << ": projecting on the real axis" << endl;
    TH1D *h1TmpCartesianX = h2Cartesian_bgSubtracted[iH2] -> ProjectionX(); 
    for(int ix=1; ix<=h1Re[iH2] -> GetNbinsX(); ix++){
      h1Re[iH2] -> SetBinContent(ix, h1TmpCartesianX -> GetBinContent(ix));
    }

    logfile << "\t\t- hit " << iH2 << ": projecting on the imaginary axis" << endl;
    TH1D *h1TmpCartesianY = h2Cartesian_bgSubtracted[iH2] -> ProjectionY(); 
    for(int iy=1; iy<=h1Im[iH2] -> GetNbinsX(); iy++){
      h1Im[iH2] -> SetBinContent(iy, h1TmpCartesianY -> GetBinContent(iy));
    }

    /* logfile << "\t\t- hit " << iH2 << ": transforming onto polar coordinate distributions" << endl; */
    /* for(int ix=1; ix<=h2Cartesian_bgSubtracted[iH2] -> GetNbinsX(); ix++){ */
    /*   for(int iy=1; iy<=h2Cartesian_bgSubtracted[iH2] -> GetNbinsY(); iy++){ */
    /* 	double re = h2Cartesian_bgSubtracted[iH2] -> GetXaxis() -> GetBinLowEdge(ix); */
    /* 	double im = h2Cartesian_bgSubtracted[iH2] -> GetYaxis() -> GetBinLowEdge(iy); */
    /* 	double mod = sqrt(re * re + im * im); */
    /* 	double ph = atan(im / re); */
    /* 	h1Mod[iH2] -> Fill(mod, h2Cartesian_bgSubtracted[iH2] -> GetBinContent(ix, iy)); */
    /* 	h1Ph[iH2] -> Fill(ph, h2Cartesian_bgSubtracted[iH2] -> GetBinContent(ix, iy)); */
    /*   } */
    /* } */

  }

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void setTimingCuts(const unsigned int nHits,
		   double *cut_Timing,
		   ofstream &logfile){ 

  cout << " - setting timing cuts " << endl;
  logfile << "\t- in function setTimingCuts()" << endl;

  for(unsigned int i=0; i<nHits; i++){
    cut_Timing[i] = TIMINGPRELIMINARYCUT;
  }

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
			   double *parLimitHigh,
			   ofstream &logfile){

  cout << " - setting function parameters " << endl;
  logfile << "\t\t- in function setFunctionParameters()" << endl;

  for(unsigned int iHit=0; iHit<nHits; iHit++){
    logfile << "\t\t- hit " << iHit << ":" << endl;
    for(unsigned int iPar=0; iPar<nPar; iPar++){
      f[iHit] -> SetParameter(iPar, par[iPar]);
      f[iHit] -> SetParName(iPar, parName[iPar].c_str());
      f[iHit] -> SetParLimits(iPar, parLimitLow[iPar], parLimitHigh[iPar]);
      logfile << "\t\t\t- parameter " << iPar
	      << ", name = " << parName[iPar]
	      << ", value = " << par[iPar]
	      << ", limits = " << parLimitLow[iPar] << " - " << parLimitHigh[iPar]
	      << endl;
    }
  }

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void checkChi2(TF1 *f, 
	       const double limitLow, 
	       const double limitHigh,
	       const char *indentLevel,
	       ofstream &logfile){

  logfile << indentLevel << "- in function checkChi2()" << endl;

  double chi2Red = f -> GetChisquare() / f -> GetNDF();
  logfile << indentLevel << "- redChi2 = " << chi2Red << endl;
  if(chi2Red < CHI2LIMITLOW || chi2Red > CHI2LIMITHIGH){
    logfile << "WARNING - reduced chi2 from fit outside optimal range [ " << limitLow << " : " << limitHigh << " ]" << endl;
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
		    double *bgRejectionErr,
		    ofstream &logfile){

  logfile << "\t- in function optimiteCutT0()" << endl;
  cout << " - optimizing T0 cut " << endl;

  TCanvas *cTmp = new TCanvas();

  const char *f_expr = "([0] / 2.) * ( (1. + TMath::Erf( (x-[1]) / (sqrt(2.) * [2]) )) * (1. + TMath::Erf( -(x-([1]+[3])) / (sqrt(2.) * [4]) )) ) + [5]";
  double f_rangeLow = T0FITRANGELOW;
  double f_rangeHigh = T0FITRANGEHIGH;
  TF1 **f = allocateFArray(nHits,
			   f_expr,
			   f_rangeLow, f_rangeHigh,
			   logfile);
  logfile << "\t- fitting function = " << f_expr << endl;
  logfile << "\t- fitting range = [ " << f_rangeLow << " : " << f_rangeHigh << " ]" << endl;

  double par[T0NPAR];
  par[0] = T0PARINITPLATEAU;
  par[1] = T0PARINITMULEFT;
  par[2] = T0PARINITSIGMALEFT;
  par[3] = T0PARINITDELTAMU;
  par[4] = T0PARINITSIGMARIGHT;
  par[5] = T0PARINITOFFSET;
  string parName[T0NPAR];
  parName[0] = "plateau #left[ns^{-1}#right]";
  parName[1] = "#mu_{L} [ns]";
  parName[2] = "#sigma_{L} [ns]";
  parName[3] = "#Delta#mu [ns]";
  parName[4] = "#sigma_{R} [ns]";
  parName[5] = "offset #left[ns^{-1}#right]";
  double parLimitLow[T0NPAR];
  parLimitLow[0] = T0PARLIMITLOWPLATEAU;
  parLimitLow[1] = T0PARLIMITLOWMULEFT;
  parLimitLow[2] = T0PARLIMITLOWSIGMALEFT;
  parLimitLow[3] = T0PARLIMITLOWDELTAMU;
  parLimitLow[4] = T0PARLIMITLOWSIGMARIGHT;
  parLimitLow[5] = T0PARLIMITLOWOFFSET;
  double parLimitHigh[T0NPAR];
  parLimitHigh[0] = T0PARLIMITHIGHPLATEAU;
  parLimitHigh[1] = T0PARLIMITHIGHMULEFT;
  parLimitHigh[2] = T0PARLIMITHIGHSIGMALEFT;
  parLimitHigh[3] = T0PARLIMITHIGHDELTAMU;
  parLimitHigh[4] = T0PARLIMITHIGHSIGMARIGHT;
  parLimitHigh[5] = T0PARLIMITHIGHOFFSET;

  logfile << "\t- setting function parameters" << endl;
  setFunctionParameters(nHits, T0NPAR, f,
			par, parName, parLimitLow, parLimitHigh,
			logfile);

  logfile << "\t- fitting T0 distributions: fit range is going to be re-adjusted" << endl;
  for(unsigned int iHit=0; iHit<nHits; iHit++){

    // finding fit range
    logfile << "\t\t- finding fit range for hit " << iHit << ": " << endl;
    TH1F *hTmp = (TH1F *) h1[iHit] -> Clone();
    hTmp -> Rebin(T0REBINFACTOR);
    double mid = hTmp -> GetBinCenter(hTmp -> GetMaximumBin());
    double max = hTmp -> GetBinContent(hTmp -> GetMaximumBin()) / T0REBINFACTOR;
    delete hTmp;
    logfile << "\t\t\t- position of the maximum = " << mid << endl;
    logfile << "\t\t\t- value of the maximum = " << max << endl;
    double rangeLow = mid - 2. *BC;
    double rangeHigh =  mid + 2. *BC;
    logfile << "\t\t\t- new fit range = [ " << rangeLow << " - " << rangeHigh << " ]" << endl;
    f[iHit] -> SetRange(rangeLow, rangeHigh);
    double plateau = max;
    double plateauLimitLow = 0.;
    double plateauLimitHigh = max;
    logfile << "\t\t\t- new plateau initialization value = " << plateau << ", range = [ " << plateauLimitLow << " - " << plateauLimitHigh << " ]" << endl;
    f[iHit] -> SetParameter(0, plateau);
    f[iHit] -> SetParLimits(0, plateauLimitLow, plateauLimitHigh);
    double muLeft = mid - 0.5 * BC;
    double muLeftLimitLow = mid-BC;
    double muLeftLimitHigh = mid;
    logfile << "\t\t\t- new muLeft initialization value = " << muLeft << ", range = [ " << muLeftLimitLow << " - " << muLeftLimitHigh << " ]" << endl;
    f[iHit] -> SetParameter(1, muLeft);
    f[iHit] -> SetParLimits(1, muLeftLimitLow, muLeftLimitHigh);
    /* double muRight = mid + 0.5 * BC; */
    /* double muRightLimitLow = mid; */
    /* double muRightLimitHigh = mid+BC; */
    /* logfile << "\t\t\t- new muRight initialization value = " << muRight << ", range = [ " << muRightLimitLow << " - " << muRightLimitHigh << " ]" << endl; */
    /* f[iHit] -> SetParameter(3, muRight); */
    /* f[iHit] -> SetParLimits(3, muRightLimitLow, muRightLimitHigh); */

    // fitting
    logfile << "\t\t- fitting " << iHit << ":" << endl;
    h1[iHit] -> Fit(f[iHit], "R && Q");
    for(unsigned int iPar=0; iPar<T0NPAR; iPar++){
      logfile << "\t\t\t- par[" << iPar << "] (" << parName[iPar] << ") = ( " << f[iHit] -> GetParameter(iPar) << " +- " << f[iHit] -> GetParError(iPar) << " )" << endl;
    }
    logfile << "\t\t\t- chi2 = " << f[iHit] -> GetChisquare() << endl;
    logfile << "\t\t\t- ndf = " << f[iHit] -> GetNDF() << endl;
    checkChi2(f[iHit], CHI2LIMITLOW, CHI2LIMITHIGH, "\t\t\t\t", logfile);

    // calculating cuts and bg rejection
    logfile << "\t\t- calculating cuts and bg rejection" << endl;
    cutLow[iHit] = f[iHit] -> GetParameter(1) - T0CUTNSIGMA * f[iHit] -> GetParameter(2);
    cutHigh[iHit] = f[iHit] -> GetParameter(1) + f[iHit] -> GetParameter(3) + T0CUTNSIGMA * f[iHit] -> GetParameter(4);
    if(cutLow[iHit] < 0.){
      logfile << "WARNING: cutLow < 0 => forcing cutLow = 0" << endl;
      cutLow[iHit] = 0.;
    }
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

void filterHistogram(TH1F *h1, 
		     const double chargeQuantile,
		     ofstream &logfile){

  logfile << "\t\t\t\t\t- in function filterHistogram()" << endl;  

  Int_t n = h1 -> GetNbinsX();

  // finding optimal frequency cut
  double cut = CHARGEFFTFREQCONVERSION * n / chargeQuantile;
  logfile << "\t\t\t\t\t- frequency cut = " << cut << endl;  

  // calculating FFT
  logfile << "\t\t\t\t\t- calculating FFT" << endl;  
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = h1 -> FFT(hm, "MAG");

  // applying frequency cut
  logfile << "\t\t\t\t\t- applying frequency cut" << endl;  
  Double_t re_full[n];
  Double_t im_full[n];
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  fft -> GetPointsComplex(re_full, im_full);
  for(unsigned i=cut; i<n-cut; i++){
    re_full[i] = 0;
    im_full[i] = 0;
  }

  // calculating backward transform
  logfile << "\t\t\t\t\t- calculating reverse FFT" << endl;  
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back -> SetPointsComplex(re_full, im_full);
  fft_back -> Transform();
  TH1 *hb = 0;

  // filling histogram of backward transform
  logfile << "\t\t\t\t\t- filling filtered histogram" << endl;  
  hb = TH1::TransformHisto(fft_back, hb, "Re");
  for(int i=0; i<n; i++){
    h1 -> SetBinContent(i, hb -> GetBinContent(i) / n);
  }

  // cleaning memory
  logfile << "\t\t\t\t\t- cleaning memory" << endl;  
  delete hm;
  delete hb;
  delete fft_back;

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

int findChargeSharingRange(TH1F *h1,
			   TH1F *h1_filtered,
			   const unsigned int algorythm, // see #define CHARGEFITRANGEFINDINGALGORITHM
			   double &min,
			   double &cross,
			   const double chargeQuantile,
			   ofstream &logfile){

  cout << " - finding charge sharing range" << endl;  
  logfile << "\t\t\t\t- in function findChargeSharingRange()" << endl;

  // transforming histogram
  TH1F *hTmp = (TH1F *) h1 -> Clone();  
  if(algorythm == 0){
    logfile << "\t\t\t\t- algorythm type: rebinned histogram" << endl;
    hTmp -> Rebin(10);
  }
  else if(algorythm == 1){
    logfile << "\t\t\t\t- algorythm type: FFT filter of histogram" << endl;
    filterHistogram(hTmp, chargeQuantile, logfile);
  }
  else{
    cout << " - ERROR!!! invalid algorythm type: " << algorythm << endl;
    logfile << "ERROR: invalid algorythm type: " << algorythm << endl;
    return 1;
  }

  // finding first maximum and first minimum
  for(int i=1; i<h1 -> GetNbinsX()-1; i++){ // first maximum
    if(hTmp -> GetBinContent(i) > hTmp -> GetBinContent(i+1)){
      min = hTmp -> GetBinLowEdge(i);
      break;
    }
  }
  for(int i=h1 -> FindBin(min); i<h1 -> GetNbinsX()-1; i++){ // first minimum
    if(hTmp -> GetBinContent(i) <= 0.) continue;
    if(hTmp -> GetBinContent(i) < hTmp -> GetBinContent(i+1)){
      cross = hTmp -> GetBinLowEdge(i+1);
      break;
    }
  }

  // copying to h1_filtered
  for(int i=1; i<hTmp -> GetNbinsX(); i++){
    h1_filtered -> SetBinContent(i, hTmp -> GetBinContent(i));
  }
  delete hTmp;

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

int fitCharge(TH1F *h1, 
	      TH1F *h1_filtered, 
	      TF1 *fTotal,
	      double &cutChargeSharing_Charge,
	      double &cutNEventsBg_Charge,
	      double &cutNEventsBgErr_Charge,
	      double &cutNEventsSignal_Charge,
	      double &cutNEventsSignalErr_Charge,
	      double &cutCSRejection_Charge,
	      double &cutCSRejectionErr_Charge,
	      ofstream &logfile){
  
  cout << " - fitting charge distribution " << endl;
  logfile << "\t\t\t- in function fitCharge()" << endl;

  // eliminating negative counts
  logfile << "\t\t\t- eliminating negative counts" << endl;  
  for(int i=0; i<h1 -> GetNbinsX(); i++){
    if(h1 -> GetBinContent(i) < 0) h1 -> SetBinContent(i, 0);
  }

  // finding signal high limit
  logfile << "\t\t\t- finding signal high limit" << endl;
  unsigned int nQuantiles = 100;
  Double_t xq[nQuantiles];
  Double_t yq[nQuantiles];
  for (unsigned int iq=0; iq<nQuantiles; iq++){
    xq[iq] = Float_t(iq+1) / nQuantiles;
  }
  h1 -> GetQuantiles(nQuantiles, yq, xq);
  TGraph *gr = new TGraph(nQuantiles, xq, yq);
  double quantileCut = 0.95;
  double chargeQuantile = gr -> Eval(quantileCut);
  logfile << "\t\t\t\t- charge value containing " << quantileCut << "\% of events = " << chargeQuantile << endl;
  delete gr;

  // finding approximate BG range
  logfile << "\t\t\t- finding approximate charge sharing background range" << endl;
  double chargeMin = -1.;
  double chargeCross = -1;
  if(findChargeSharingRange(h1, h1_filtered, CHARGEFITRANGEFINDINGALGORITHM,
			    chargeMin, chargeCross,
			    chargeQuantile,
			    logfile)){
    cout << " - ERROR!!! cannot find charge sharing range" << endl;
    logfile << "ERROR: cannot find charge sharing range" << endl;
    return 1;
  }
  cout << " - chargeMin = " << chargeMin << endl;
  logfile << "\t\t\t\t- chargeMin = " << chargeMin << endl;
  cout << " - chargeCross = " << chargeCross << endl;
  logfile << "\t\t\t\t- chargeCross = " << chargeCross << endl;
  if(chargeMin >= chargeCross){
    cout << " - ERROR!!! chargeMin >= chargeCross" << endl;
    logfile << "ERROR: chargeMin >= chargeCross" << endl;
    return 1;
  }
  if(chargeMin < 0.){
    cout << " - WARNING!!! - chargeMin invalid value: " << chargeMin << endl;
    logfile << "WARNING: chargeMin invalid value: " << chargeMin << endl;
  }
  if(chargeCross < 0.){
    cout << " - WARNING!!! - chargeCross invalid value: " << chargeCross << endl;
    logfile << "WARNING: chargeCross invalid value: " << chargeCross << endl;
  }

  // temporary fit of the BG only
  logfile << "\t\t\t- temporary fit to charge sharing tail only" << endl;
  const char *fBg_expr = "[0]/(x-[1])";
  double fBg_rangeLow = chargeMin + 0.5*(chargeCross-chargeMin);
  double fBg_rangeHigh = chargeCross;
  TF1 *fBg = new TF1("fBg", 
		     fBg_expr,
		     fBg_rangeLow, fBg_rangeHigh);
  logfile << "\t\t\t\t- fitting function = " << fBg_expr << endl;
  logfile << "\t\t\t\t- fitting range = [ " << fBg_rangeLow << " : " << fBg_rangeHigh << " ]" << endl;
  logfile << "\t\t\t\t- uninitialized parameters (this may be optimized)" << endl;
  h1 -> Fit(fBg, "R && Q");
  double p0 = fBg -> GetParameter(0);
  double p1 = fBg -> GetParameter(1);
  logfile << "\t\t\t\t- result: " << endl;
  logfile << "\t\t\t\t\t- p0 = " << p0 << endl;
  logfile << "\t\t\t\t\t- p1 = " << p1 << endl;
  logfile << "\t\t\t\t\t- chi2 = " << fBg -> GetChisquare() << endl;
  logfile << "\t\t\t\t\t- ndf = " << fBg -> GetNDF() << endl;
  delete fBg;

  //  temporary fit of the signal only
  logfile << "\t\t\t- temporary fit to the signal distribution only" << endl;
  const char *fSignal_expr = "( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) )";
  double fSignal_rangeLow = chargeCross;
  double fSignal_rangeHigh = chargeQuantile;
  TF1 *fSignal = new TF1("fSignal", 
			 fSignal_expr,
			 fSignal_rangeLow, fSignal_rangeHigh);
  logfile << "\t\t\t\t- fitting function = " << fSignal_expr << endl;
  logfile << "\t\t\t\t- fitting range = [ " << fSignal_rangeLow << " : " << fSignal_rangeHigh << " ]" << endl;
  logfile << "\t\t\t\t- initializing parameters:" << endl;
  h1 -> GetXaxis() -> SetRangeUser(fSignal_rangeLow, fSignal_rangeHigh);
  double scaleInit = h1 -> Integral();
  double muInit = h1 -> GetBinCenter(h1 -> GetMaximumBin());
  double sigmaInit = muInit / 2.;
  logfile << "\t\t\t\t\t- scale = " << scaleInit << endl;
  logfile << "\t\t\t\t\t- mu = " << muInit << endl;
  logfile << "\t\t\t\t\t- sigma = " << sigmaInit << endl;
  fSignal -> SetParameter(0, scaleInit);
  fSignal -> SetParameter(1, muInit);
  fSignal -> SetParameter(2, sigmaInit);
  h1 -> Fit(fSignal, "R && Q");

  // final fit BG+signal
  logfile << "\t\t\t- final fit" << endl;
  double fTotal_rangeLow = fBg_rangeLow;
  double fTotal_rangeHigh = fSignal_rangeHigh;
  fTotal -> SetRange(fTotal_rangeLow, fTotal_rangeHigh);
  logfile << "\t\t\t\t- fitting range = [ " << fTotal_rangeLow << " : " << fTotal_rangeHigh << " ]" << endl;
  logfile << "\t\t\t\t- initializing parameters from previous temporary fits" << endl;
  fTotal -> SetParameter(0, fSignal -> GetParameter(0));
  fTotal -> SetParameter(1, fSignal -> GetParameter(1));
  fTotal -> SetParameter(2, fSignal -> GetParameter(2));
  fTotal -> SetParameter(3, p0);
  fTotal -> SetParameter(4, p1);
  h1 -> Fit(fTotal, "R && Q");
  checkChi2(fTotal, CHI2LIMITLOW, CHI2LIMITHIGH, "\t\t\t\t\t", logfile);

  // finding BG-signal intersection point
  logfile << "\t\t\t- finding intersection point between charge sharing and signal distributions" << endl;
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

  // computing rejection power
  logfile << "\t\t\t- computing rejection power (yet to be implemented)" << endl;
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

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

TF1 **optimizeCutCharge(const unsigned int nHits,
			TH1F **h1_bg,
			TH2F **h2,
			TH2F **h2_bgSubtracted,
			TH1F **h1_bgSubtracted,
			TH1F **h1_bgSubtracted_filtered,
			double *cutLow_T0,
			double *cutHigh_T0,
			double *cutChargeSharing_Charge,
			double *cutNEventsBg_Charge,
			double *cutNEventsBgErr_Charge,
			double *cutNEventsSignal_Charge,
			double *cutNEventsSignalErr_Charge,
			double *cutCSRejection_Charge,
			double *cutCSRejectionErr_Charge,
			ofstream &logfile){

  cout << " - optimizing charge cut " << endl;
  logfile << "\t- in function optimizeCutCharge()" << endl;

  TCanvas *cTmp = new TCanvas();

  const char *f_expr = "( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) ) + ([3] / (x-[4])) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[1]))) / 2. ";
  TF1 **f = allocateFArray(nHits,
			   f_expr,
			   CHARGEFITRANGELOW,
			   CHARGEFITRANGEHIGH,
			   logfile);
  logfile << "\t- fitting function = " << f_expr << endl;

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

  logfile << "\t- setting function parameters" << endl;
  setFunctionParameters(nHits, CHARGENPAR, f,
			par, parName, parLimitLow, parLimitHigh,
			logfile);

  logfile << "\t- building charge distribution and fitting" << endl;
  for(unsigned int iHit=0; iHit<nHits; iHit++){

    // calculating number of bins outside CHARGE cut region
    double nBinsT0 = T0NBINS * (cutLow_T0[iHit] - T0MIN + T0MAX - cutHigh_T0[iHit]) / (T0MAX - T0MIN);
    logfile << "\t\t- hit " << iHit << ": number of T0 bins outside the T0 cut region (background region) = " << nBinsT0 << endl;
    cout << " - channel " << iHit << ": "
	 << "nBinsT0 = " << nBinsT0
	 << endl;

    // subtracting background from 2D histogram
    logfile << "\t\t- hit " << iHit << ": subtracting background from 2D histogram" << endl;
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
    logfile << "\t\t- hit " << iHit << ": projecting signal region onto 1D histogram" << endl;
    TH1D *h1Tmp = h2_bgSubtracted[iHit] -> ProjectionY("tmp", 
						       h2_bgSubtracted[iHit] -> GetXaxis() -> FindBin(cutLow_T0[iHit]),
						       h2_bgSubtracted[iHit] -> GetXaxis() -> FindBin(cutHigh_T0[iHit]),
						       "e");
    for(unsigned int iY=0; iY<CHARGENBINS; iY++){
      h1_bgSubtracted[iHit] -> SetBinContent(iY, h1Tmp -> GetBinContent(iY));
    }
    delete h1Tmp;

    // fitting signal
    logfile << "\t\t- hit " << iHit << ": fitting bg-subtracted charge distribution" << endl;
    if(fitCharge(h1_bgSubtracted[iHit], h1_bgSubtracted_filtered[iHit], f[iHit],
		 cutChargeSharing_Charge[iHit],
		 cutNEventsBg_Charge[iHit],
		 cutNEventsBgErr_Charge[iHit],
		 cutNEventsSignal_Charge[iHit],
		 cutNEventsSignalErr_Charge[iHit],
		 cutCSRejection_Charge[iHit],
		 cutCSRejectionErr_Charge[iHit],
		 logfile)){
      cout << " - ERROR!!! cannot fit charge" << endl;
      return NULL;
    }

  }

  delete cTmp;
  return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

double getQuantile(TH1F *h1, 
		   const unsigned int nQuantiles, 
		   const double quantileCut){
  Double_t xq[nQuantiles];
  Double_t yq[nQuantiles];
  for (unsigned int iq=0; iq<nQuantiles; iq++){
    xq[iq] = Float_t(iq+1) / nQuantiles;
  }
  h1 -> GetQuantiles(nQuantiles, yq, xq);
  TGraph *gr = new TGraph(nQuantiles, xq, yq);
  double quantile = gr -> Eval(quantileCut);
  gr -> SaveAs("gr.root");
  delete gr;
  return quantile;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void optimizeCutTiming(const unsigned int nHits,
		       TH1F **h1,
		       TGraph **gr_quantiles,
		       double *cutLow,
		       double *cutHigh,
		       ofstream &logfile){

  cout << " - optimizing Timing cut" << endl;
  logfile << "\t- in function optimizingCutTiming()" << endl;

  for(unsigned int iHit=0; iHit<nHits; iHit++){

    // first fit: simple Gaussian, to estimate the peak position
    logfile << "\t\t- hit " << iHit << ": fitting peak with Gaussian:" << endl;
    h1[iHit] -> GetXaxis() -> SetRangeUser(TIMINGPRELIMINARYCUT, TIMINGMAX);
    double scaleMax = h1[iHit] -> GetBinContent(h1[iHit] -> GetMaximumBin());
    double mu = h1[iHit] -> GetBinCenter(h1[iHit] -> GetMaximumBin());
    h1[iHit] -> GetXaxis() -> SetRangeUser(TIMINGMIN, TIMINGMAX);
    double sigma = 6.;
    logfile << "\t\t\t- initializing parameters:" << endl;
    logfile << "\t\t\t\t- scaleMax = " << scaleMax << endl;
    logfile << "\t\t\t\t- mu = " << mu << endl;
    logfile << "\t\t\t\t- sigma = " << sigma << endl;
    TF1 *gauss = new TF1("gauss", "gaus(0)", mu-5*sigma, mu+5*sigma);
    gauss -> SetParLimits(0, 0., scaleMax);
    gauss -> SetParameter(0, scaleMax/2.);
    gauss -> SetParameter(1, mu);
    gauss -> SetParameter(2, sigma);
    gauss -> SetParLimits(2, 0, 100.);
    h1[iHit] -> Fit(gauss, "Q && R");
    if(gauss -> GetParameter(1) <= 0.){
      logfile << "\t\t\t- first fit attempt failed, trying to fit in a narrower range" << endl;
      gauss -> SetRange(mu-1.5*sigma, mu+1.5*sigma);
      h1[iHit] -> Fit(gauss, "Q && R");
    }
    double scale = gauss -> GetParameter(0);
    mu = gauss -> GetParameter(1);
    sigma = gauss -> GetParameter(2);
    delete gauss;

    // second fit: adding background
    gauss = new TF1("gauss", "gaus(0)+pol1(3)", mu-2*sigma, mu+4*sigma);
    gauss -> SetParameter(0, scale);
    gauss -> SetParameter(1, mu);
    gauss -> SetParameter(2, sigma);
    gauss -> SetParLimits(2, 0, 100.);
    gauss -> SetParLimits(3, 0, 1000000.);
    h1[iHit] -> Fit(gauss, "Q && R");
    scale = gauss -> GetParameter(0);
    mu = gauss -> GetParameter(1);
    sigma = gauss -> GetParameter(2);
    double scaleErr = gauss -> GetParError(0);
    double muErr = gauss -> GetParError(1);
    double sigmaErr = gauss -> GetParError(2);
    double chi2 = gauss -> GetChisquare();
    double ndf = gauss -> GetNDF();
    delete gauss;
    logfile << "\t\t\t- parameters values:" << endl;
    logfile << "\t\t\t\t- scale = ( " << scale << " +- " << scaleErr << " )" << endl;
    logfile << "\t\t\t\t- mu = ( " << mu << " +- " << muErr << " )" << endl;
    logfile << "\t\t\t\t- sigma = ( " << sigma << " +- " << sigmaErr << " )" << endl;
    logfile << "\t\t\t\t- chi2/ndf = " << chi2/ndf << endl;
    cutLow[iHit] = mu - TIMINGCUTNSIGMA*sigma;
    cutHigh[iHit] = mu + TIMINGCUTNSIGMA*sigma;
    if(cutLow[iHit] < 0.) cutLow[iHit] = 0.;

  }

  return ;
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
		    0.9*h1 -> GetXaxis() -> GetBinLowEdge(1), cutLow,
		    yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));
  fillExclusionPlot(h1_exclusionRight,
		    cutHigh, 1.1*h1 -> GetBinLowEdge(h1 -> GetNbinsX()),
		    yScaleFactor * h1 -> GetBinContent(h1 -> GetMaximumBin()));

  h1 -> Draw();
  double rangeMin = cutLow - 1.5*BC;
  double rangeMax = cutLow + 4.5*BC;
  if(rangeMin < 0) rangeMin = 0.;
  h1 -> GetXaxis() -> SetRangeUser(rangeMin, rangeMax);
  h1 -> GetYaxis() -> SetRangeUser(0, yScaleFactor * h1_cut -> GetBinContent(h1_cut -> GetMaximumBin()));
  h1_cut -> Draw("same");
  h1_exclusionLeft -> Draw("same");
  h1_exclusionRight -> Draw("same");

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

  gPad -> RedrawAxis();
  
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
  gPad -> RedrawAxis();

  leg = new TLegend(CHARGELEGENDXLOW, CHARGELEGENDYLOW, CHARGELEGENDXHIGH, CHARGELEGENDYHIGH);
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

  const double yScaleFactor = 2.5;

  gStyle -> SetOptFit(0);

  h1 -> GetXaxis() -> SetRangeUser(0., f -> GetParameter(2) + 5. * f -> GetParameter(1));
  h1 -> GetYaxis() -> SetRangeUser(0., yScaleFactor * f -> Eval(f -> GetParameter(2)));
  h1 -> Draw();

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

  gPad -> RedrawAxis();

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotTiming(TH1F *h1, 
		TH1F *h1Cut,
		const double cutLow,
		const  double cutHigh,
		TLegend *leg,
		TPaveStats *ptresults,
		TH1F *h1_exclusionLeft,
		TH1F *h1_exclusionRight){

  cout << " - plotting timing" << endl;

  const double yScaleFactor = 1.1;

  h1 -> GetYaxis() -> SetRangeUser(0., yScaleFactor * h1Cut -> GetBinContent(h1Cut -> GetMaximumBin()));
  h1Cut -> GetYaxis() -> SetRangeUser(0., yScaleFactor * h1Cut -> GetBinContent(h1Cut -> GetMaximumBin()));
  double rangeMin = cutLow - 1.5*BC;
  double rangeMax = cutLow + 4.5*BC;
  if(rangeMin < 0) rangeMin = 0.;
  //  h1 -> GetXaxis() -> SetRangeUser(rangeMin, rangeMax);
  //  h1Cut -> GetXaxis() -> SetRangeUser(rangeMin, rangeMax);

  h1 -> Draw();
  h1Cut -> Draw("same");
  gPad -> RedrawAxis();

  fillExclusionPlot(h1_exclusionLeft,
                    h1 -> GetXaxis() -> GetBinLowEdge(1), cutLow,
		    yScaleFactor * h1 -> GetBinContent(h1Cut -> GetMaximumBin()));
  h1_exclusionLeft -> Draw("same");

  fillExclusionPlot(h1_exclusionRight,
		    cutHigh, 1.1 * h1 -> GetBinLowEdge(h1 -> GetNbinsX()),
                    yScaleFactor * h1 -> GetBinContent(h1Cut -> GetMaximumBin()));
  h1_exclusionRight -> Draw("same");

  leg = new TLegend(TIMINGLEGENDXLOW, TIMINGLEGENDYLOW, TIMINGLEGENDXHIGH, TIMINGLEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(1);
  char text[200];
  leg -> AddEntry(h1, "no cut", "f");
  leg -> AddEntry(h1Cut, "T_{0} and charge cuts", "f");
  leg -> AddEntry(h1_exclusionLeft, "cut region", "f");
  leg -> Draw();

  ptresults = new TPaveStats(TIMINGRESULTSXLOW, TIMINGRESULTSYLOW, TIMINGRESULTSXHIGH, TIMINGRESULTSYHIGH, "brNDC");
  ptresults -> SetTextSize(0.025);
  ptresults -> SetName("results");
  ptresults -> SetBorderSize(1);
  ptresults -> SetFillColor(LEGENDFILLCOLOR);
  ptresults -> SetLineColor(1);
  ptresults -> SetTextAlign(12);
  ptresults -> SetTextFont(42);
  sprintf(text, "Timing cut #in [%.1lf #pm %.1lf] ns", cutLow, cutHigh);
  ptresults -> AddText(text);
  ptresults -> SetOptStat(0);
  ptresults -> SetOptFit(111);
  ptresults -> Draw();

  gPad -> RedrawAxis();

  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void plotChargeAfterCuts(TH1F *h1, 
			 TH1F *h1Cut,
			 TLegend *leg,
			 TPaveStats *ptstats,
			 const double cutLow_T0,
			 const double cutHigh_T0,
			 const double cutLow_Timing,
			 const double cutHigh_Timing){

  h1 -> Draw();
  h1Cut -> Draw("same");

  leg = new TLegend(CHARGELEGENDXLOW, CHARGELEGENDYLOW, CHARGELEGENDXHIGH, CHARGELEGENDYHIGH);
  leg -> SetFillColor(LEGENDFILLCOLOR);
  leg -> SetLineColor(LEGENDLINECOLOR);
  char text[200];
  leg -> AddEntry(h1, "no cut", "f");
  leg -> AddEntry(h1Cut, "T_{0} and timing cuts", "f");
  leg -> Draw();

  ptstats = new TPaveStats(CHARGERESULTSXLOW, CHARGERESULTSYLOW, CHARGERESULTSXHIGH, CHARGERESULTSYHIGH, "brNDC");
  ptstats -> SetTextSize(0.03);
  ptstats -> SetName("results");
  ptstats -> SetBorderSize(1);
  ptstats -> SetFillColor(LEGENDFILLCOLOR);
  ptstats -> SetLineColor(LEGENDLINECOLOR);
  ptstats -> SetTextAlign(12);
  ptstats -> SetTextFont(42);
  sprintf(text, "T_{0} cut #in [%.1lf #pm %.1lf] ns", cutLow_T0, cutHigh_T0);
  ptstats -> AddText(text);
  sprintf(text, "Timing cut #in [%.1lf #pm %.1lf] ns", cutLow_Timing, cutHigh_Timing);
  ptstats -> AddText(text);
  ptstats -> SetOptStat(0);
  ptstats -> SetOptFit(111);
  ptstats -> Draw();


  gPad -> RedrawAxis();
  
  return ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void savePlots(const unsigned int nH2,
	       plot_struct plots,
	       double *cut_Timing,
	       const unsigned int T0Npar,
	       TF1 **f_T0,
	       double *cutLow_T0,
	       double *cutHigh_T0,
	       double *cutBgRejection_T0,
	       double *cutBgRejectionErr_T0,
	       TF1 **f_Charge,
	       double *cutChargeSharing_Charge,
	       TGraph **gr_Timing_T0Cut_ChargeCut_quantiles,
	       double *cutLow_Timing,
	       double *cutHigh_Timing,
	       const string outputFolder,
	       const string flag,
	       const bool drawPlots,
	       const bool holdPlots,
	       ofstream &logfile){

  cout << " - drawing " << endl;

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(0);

  // allocating canvases
  TCanvas **cc_Harm1Mod = allocateCanvasArray(nH2, "cc_Harm1Mod_channel_%d",
  					      0, 0, CANVASSIZE, CANVASSIZE,
  					      false, false, false,
  					      logfile);
  TCanvas **cc_Harm1Ph = allocateCanvasArray(nH2, "cc_Harm1Ph_channel_%d",
  					      0, 0, CANVASSIZE, CANVASSIZE,
  					      false, false, false,
  					      logfile);
  TCanvas **cc_Harm1Re = allocateCanvasArray(nH2, "cc_Harm1Re_channel_%d",
  					      0, 0, CANVASSIZE, CANVASSIZE,
  					      false, false, false,
  					      logfile);
  TCanvas **cc_Harm1Im = allocateCanvasArray(nH2, "cc_Harm1Im_channel_%d",
  					      0, 0, CANVASSIZE, CANVASSIZE,
  					      false, false, false,
  					      logfile);
  TCanvas **cc_Harm1Polar = allocateCanvasArray(nH2, "cc_Harm1Polar_channel_%d",
						0, 0, CANVASSIZE, CANVASSIZE,
						false, false, false,
						logfile);
  TCanvas **cc_Harm1Polar_bgSubtracted = allocateCanvasArray(nH2, "cc_Harm1Polar_bgSubtracted_channel_%d",
							     0, 0, CANVASSIZE, CANVASSIZE,
							     false, false, false,
							     logfile);
  TCanvas **cc_Harm1Cartesian = allocateCanvasArray(nH2, "cc_Harm1Cartesian_channel_%d",
						    0, 0, CANVASSIZE, CANVASSIZE,
						    false, false, false,
						    logfile);
  TCanvas **cc_Harm1Cartesian_bgSubtracted = allocateCanvasArray(nH2, "cc_Harm1Cartesian_bgSubtracted_channel_%d",
								 0, 0, CANVASSIZE, CANVASSIZE,
								 false, false, false,
								 logfile);
  TCanvas **cc_Harm1Mod_vs_T0 = allocateCanvasArray(nH2, "cc_Harm1Mod_vs_T0_channel_%d",
  						  0, 0, CANVASSIZE, CANVASSIZE,
  						  false, false, true,
						  logfile);
  TCanvas **cc_Timing_vs_T0 = allocateCanvasArray(nH2, "cc_Timing_vs_T0_channel_%d",
  						  0, 0, CANVASSIZE, CANVASSIZE,
  						  false, true, true,
						  logfile);

  TCanvas **cc_Charge_vs_T0 = allocateCanvasArray(nH2, "cc_Charge_vs_T0_channel_%d",
  						  CANVASSIZE, 0, CANVASSIZE, CANVASSIZE,
  						  false, true, true,
						  logfile);

  TCanvas **cc_Charge_vs_Timing = allocateCanvasArray(nH2, "cc_Charge_vs_Timing_channel_%d",
  						      2.*CANVASSIZE, 0, CANVASSIZE, CANVASSIZE,
  						      false, true, true,
						      logfile);

  TCanvas **cc_T0 = allocateCanvasArray(nH2, "cc_T0_channel_%d",
  					0, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					false, false, false,
					logfile);

  TCanvas **cc_Charge = allocateCanvasArray(nH2, "cc_Charge_channel_%d",
  					    CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					    true, true, false,
					    logfile);

  TCanvas **cc_Charge_T0Cut_TimingCut = allocateCanvasArray(nH2, "cc_Charge_T0Cut_TimingCut_channel_%d",
							    CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
							    true, true, false,
							    logfile);

  TCanvas **cc_Charge_vs_T0_bgSubtracted = allocateCanvasArray(nH2, "cc_Charge_vs_T0_bgSubtracted_channel_%d",
  							       2.*CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  							       false, true, false,
							       logfile);

  TCanvas **cc_Charge_bgSubtracted = allocateCanvasArray(nH2, "cc_Charge_bgSubtracted_channel_%d",
  							 2.*CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  							 false, false, false,
							 logfile);

  TCanvas **cc_Charge_bgSubtracted_filtered = allocateCanvasArray(nH2, "cc_Charge_bgSubtracted_filtered_channel_%d",
								  2.*CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
								  false, false, false,
								  logfile);

  TCanvas **cc_Timing = allocateCanvasArray(nH2, "cc_Timing_channel_%d",
  					    CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  					    false, false, false,
					    logfile);

  TCanvas **cc_Timing_quantiles = allocateCanvasArray(nH2, "cc_Timing_quantiles_channel_%d",
  						      CANVASSIZE, CANVASSIZE, CANVASSIZE, CANVASSIZE,
  						      false, false, false,
						      logfile);

  // allocating functions
  TF1 **f_ChargeBg = allocateFArray(nH2,
				    "([0] / (x-[1])) * (1. + TMath::Erf((-(x-[2])) / (sqrt(2.)*[3]))) / 2.",
				    CHARGEFITRANGELOW,
				    CHARGEFITRANGEHIGH,
				    logfile);

  TF1 **f_ChargeSignal = allocateFArray(nH2,
					"( [0] / ( sqrt(2. * TMath::Pi()) * [1] ) ) * exp( - (x-[2]) * (x-[2]) / (2. * [1] * [1]) )",				       
					CHARGEFITRANGELOW,
					CHARGEFITRANGEHIGH,
					logfile);
  
  // allocating TLegends and TPaves
  TLegend **leg_T0 = new TLegend*[nH2];
  TPaveStats **ptstats_T0 = new TPaveStats*[nH2];
  TPaveStats **ptstats_Charge_T0Cut_TimingCut = new TPaveStats*[nH2];
  TPaveStats **ptresults_T0 = new TPaveStats*[nH2];
  TLegend **leg_Charge = new TLegend*[nH2];
  TLegend **leg_Charge_T0Cut_TimingCut = new TLegend*[nH2];
  TLegend **leg_ChargeFit = new TLegend*[nH2];
  TPaveStats **ptstats_ChargeFit = new TPaveStats*[nH2];
  TLegend **legExclusion_ChargeFit = new TLegend*[nH2];
  TLegend **leg_Timing = new TLegend*[nH2];
  TPaveStats **ptstats_Timing = new TPaveStats*[nH2];

  // creating output file
  TFile *file = TFile::Open((outputFolder + flag + ".root").c_str(), "RECREATE");

  // drawing
  for(unsigned int iH=0; iH<nH2; iH++){

    stringstream channel;
    channel << iH;

    cc_Harm1Mod[iH] -> cd();
    plots.h1_Harm1Mod[iH] -> SetMinimum(0);
    plots.h1_Harm1Mod[iH] -> Draw();
    plots.h1_Harm1Mod_bgSubtracted[iH] -> Draw("same");
    gPad -> RedrawAxis();
    cc_Harm1Mod[iH] -> Write();
    if(drawPlots){
      cc_Harm1Mod[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod.png").c_str());
      cc_Harm1Mod[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod.pdf").c_str());
      cc_Harm1Mod[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod.eps").c_str());
    }

    cc_Harm1Ph[iH] -> cd();
    plots.h1_Harm1Ph[iH] -> Draw();
    plots.h1_Harm1Ph[iH] -> SetMinimum(0);
    plots.h1_Harm1Ph_bgSubtracted[iH] -> Draw("same");
    gPad -> RedrawAxis();
    cc_Harm1Ph[iH] -> Write();
    if(drawPlots){
      cc_Harm1Ph[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Ph.png").c_str());
      cc_Harm1Ph[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Ph.pdf").c_str());
      cc_Harm1Ph[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Ph.eps").c_str());
    }

    cc_Harm1Re[iH] -> cd();
    plots.h1_Harm1Re[iH] -> Draw();
    plots.h1_Harm1Re_bgSubtracted[iH] -> Draw("same");
    gPad -> RedrawAxis();
    cc_Harm1Re[iH] -> Write();
    if(drawPlots){
      cc_Harm1Re[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Re.png").c_str());
      cc_Harm1Re[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Re.pdf").c_str());
      cc_Harm1Re[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Re.eps").c_str());
    }

    cc_Harm1Im[iH] -> cd();
    plots.h1_Harm1Im[iH] -> Draw();
    plots.h1_Harm1Im_bgSubtracted[iH] -> Draw("same");
    gPad -> RedrawAxis();
    cc_Harm1Im[iH] -> Write();
    if(drawPlots){
      cc_Harm1Im[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Im.png").c_str());
      cc_Harm1Im[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Im.pdf").c_str());
      cc_Harm1Im[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Im.eps").c_str());
    }

    cc_Harm1Polar[iH] -> cd();
    plots.h2_Harm1Polar[iH] -> Draw("colz");
    cc_Harm1Polar[iH] -> Write();
    if(drawPlots){
      cc_Harm1Polar[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar.png").c_str());
      cc_Harm1Polar[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar.pdf").c_str());
      cc_Harm1Polar[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar.eps").c_str());
    }

    cc_Harm1Polar_bgSubtracted[iH] -> cd();
    plots.h2_Harm1Polar_bgSubtracted[iH] -> SetMinimum(0);
    plots.h2_Harm1Polar_bgSubtracted[iH] -> SetMaximum(plots.h2_Harm1Polar[iH] -> GetBinContent(plots.h2_Harm1Polar[iH] -> GetMaximumBin()));
    plots.h2_Harm1Polar_bgSubtracted[iH] -> Draw("colz");
    cc_Harm1Polar_bgSubtracted[iH] -> Write();
    if(drawPlots){
      cc_Harm1Polar_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar_bgSubtracted.png").c_str());
      cc_Harm1Polar_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar_bgSubtracted.pdf").c_str());
      cc_Harm1Polar_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Polar_bgSubtracted.eps").c_str());
    }

    cc_Harm1Cartesian[iH] -> cd();
    plots.h2_Harm1Cartesian[iH] -> Draw("colz");
    cc_Harm1Cartesian[iH] -> Write();
    if(drawPlots){
      cc_Harm1Cartesian[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian.png").c_str());
      cc_Harm1Cartesian[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian.pdf").c_str());
      cc_Harm1Cartesian[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian.eps").c_str());
    }

    cc_Harm1Cartesian_bgSubtracted[iH] -> cd();
    plots.h2_Harm1Cartesian_bgSubtracted[iH] -> SetMinimum(0);
    plots.h2_Harm1Cartesian_bgSubtracted[iH] -> SetMaximum(plots.h2_Harm1Cartesian[iH] -> GetBinContent(plots.h2_Harm1Cartesian[iH] -> GetMaximumBin()));
    plots.h2_Harm1Cartesian_bgSubtracted[iH] -> Draw("colz");
    cc_Harm1Cartesian_bgSubtracted[iH] -> Write();
    if(drawPlots){
      cc_Harm1Cartesian_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian_bgSubtracted.png").c_str());
      cc_Harm1Cartesian_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian_bgSubtracted.pdf").c_str());
      cc_Harm1Cartesian_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Cartesian_bgSubtracted.eps").c_str());
    }

    cc_Harm1Mod_vs_T0[iH] -> cd();
    plots.h2_Harm1Mod_vs_T0[iH] -> Draw("colz");
    cc_Harm1Mod_vs_T0[iH] -> Write();
    if(drawPlots){
      cc_Harm1Mod_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod_vs_T0.png").c_str());
      cc_Harm1Mod_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod_vs_T0.pdf").c_str());
      cc_Harm1Mod_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Harm1Mod_vs_T0.eps").c_str());
    }

    cc_Timing_vs_T0[iH] -> cd();
    plots.h2_Timing_vs_T0[iH] -> Draw("colz");
    cc_Timing_vs_T0[iH] -> Write();
    if(drawPlots){
      cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.png").c_str());
      cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.pdf").c_str());
      cc_Timing_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_vs_T0.eps").c_str());
    }

    cc_Charge_vs_T0[iH] -> cd();
    plots.h2_Charge_vs_T0[iH] -> Draw("colz");
    cc_Charge_vs_T0[iH] -> Write();
    if(drawPlots){
      cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.png").c_str());
      cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.pdf").c_str());
      cc_Charge_vs_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0.eps").c_str());
    }

    cc_Charge_vs_Timing[iH] -> cd();
    plots.h2_Charge_vs_Timing[iH] -> Draw("colz");
    cc_Charge_vs_Timing[iH] -> Write();
    if(drawPlots){
      cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.png").c_str());
      cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.pdf").c_str());
      cc_Charge_vs_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_Timing.eps").c_str());
    }

    cc_T0[iH] -> cd();
    plotT0(plots.h1_T0[iH], plots.h1_T0_TimingCut[iH], f_T0[iH],
    	   leg_T0[iH], ptstats_T0[iH], ptresults_T0[iH],
    	   plots.h1_exclusionLeft_T0[iH], plots.h1_exclusionRight_T0[iH],
    	   cut_Timing[iH], cutLow_T0[iH], cutHigh_T0[iH],
    	   cutBgRejection_T0[iH], cutBgRejectionErr_T0[iH],
    	   T0Npar);
    cc_T0[iH] -> Write();
    // tmp
    cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.png").c_str());
    if(drawPlots){
      cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.png").c_str());
      cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.pdf").c_str());
      cc_T0[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_T0.eps").c_str());
    }

    cc_Charge[iH] -> cd();
    plotCharge(plots.h1_Charge[iH], plots.h1_Charge_T0Cut_bg[iH],
    	       leg_Charge[iH],
    	       cutLow_T0[iH], cutHigh_T0[iH]);
    cc_Charge[iH] -> Write();
    if(drawPlots){
      cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.png").c_str());
      cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.pdf").c_str());
      cc_Charge[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge.eps").c_str());
    }

    cc_Charge_T0Cut_TimingCut[iH] -> cd();
    plotChargeAfterCuts(plots.h1_Charge[iH], plots.h1_Charge_T0Cut_TimingCut[iH],
			leg_Charge_T0Cut_TimingCut[iH],
			ptstats_Charge_T0Cut_TimingCut[iH],
			cutLow_T0[iH], cutHigh_T0[iH],
			cutLow_Timing[iH], cutHigh_Timing[iH]);
    cc_Charge_T0Cut_TimingCut[iH] -> Write();
    // tmp
    cc_Charge_T0Cut_TimingCut[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_T0Cut_TimingCut.png").c_str());
    if(drawPlots){
      cc_Charge_T0Cut_TimingCut[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_T0Cut_TimingCut.png").c_str());
      cc_Charge_T0Cut_TimingCut[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_T0Cut_TimingCut.pdf").c_str());
      cc_Charge_T0Cut_TimingCut[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_T0Cut_TimingCut.eps").c_str());
    }

    cc_Charge_vs_T0_bgSubtracted[iH] -> cd();
    plots.h2_Charge_vs_T0_bgSubtracted[iH] -> SetMinimum(0);
    plots.h2_Charge_vs_T0_bgSubtracted[iH] -> Draw("colz");
    cc_Charge_vs_T0_bgSubtracted[iH] -> Write();
    if(drawPlots){
      cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.png").c_str());
      cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.pdf").c_str());
      cc_Charge_vs_T0_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_vs_T0_bgSubtracted.eps").c_str());
    }

    cc_Charge_bgSubtracted[iH] -> cd();
    plotChargeFit(plots.h1_Charge_bgSubtracted[iH], 
		  f_Charge[iH], f_ChargeBg[iH], f_ChargeSignal[iH],
		  plots.h1_exclusion_Charge[iH], cutChargeSharing_Charge[iH],
		  leg_ChargeFit[iH], ptstats_ChargeFit[iH], legExclusion_ChargeFit[iH],
		  CHARGENPAR);
    cc_Charge_bgSubtracted[iH] -> Write();
    if(drawPlots){
      cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.png").c_str());
      cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.pdf").c_str());
      cc_Charge_bgSubtracted[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted.eps").c_str());
    }

    cc_Charge_bgSubtracted_filtered[iH] -> cd();
    plots.h1_Charge_bgSubtracted_filtered[iH] -> Draw();
    cc_Charge_bgSubtracted_filtered[iH] -> Write();
    if(drawPlots){
      cc_Charge_bgSubtracted_filtered[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted_filtered.png").c_str());
      cc_Charge_bgSubtracted_filtered[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted_filtered.pdf").c_str());
      cc_Charge_bgSubtracted_filtered[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Charge_bgSubtracted_filtered.eps").c_str());
    }

    cc_Timing[iH] -> cd();
    plotTiming(plots.h1_Timing[iH], plots.h1_Timing_T0Cut_ChargeCut[iH],
    	       cutLow_Timing[iH], cutHigh_Timing[iH],
    	       leg_Timing[iH], ptstats_Timing[iH],
    	       plots.h1_exclusionLeft_Timing[iH], plots.h1_exclusionRight_Timing[iH]);
    cc_Timing[iH] -> Write();
    // tmp
    cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.png").c_str());
    if(drawPlots){
      cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.png").c_str());
      cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.pdf").c_str());
      cc_Timing[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing.eps").c_str());
    }

    cc_Timing_quantiles[iH] -> cd();
    gr_Timing_T0Cut_ChargeCut_quantiles[iH] -> GetXaxis() -> SetTitle("fractional area");
    gr_Timing_T0Cut_ChargeCut_quantiles[iH] -> GetYaxis() -> SetTitle("Timing [ns]");
    gr_Timing_T0Cut_ChargeCut_quantiles[iH] -> SetLineWidth(3);
    gr_Timing_T0Cut_ChargeCut_quantiles[iH] -> Draw("al");
    cc_Timing_quantiles[iH] -> Write();
    if(drawPlots){
      cc_Timing_quantiles[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_quantiles.png").c_str());
      cc_Timing_quantiles[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_quantiles.pdf").c_str());
      cc_Timing_quantiles[iH] -> SaveAs((outputFolder + flag + "_channel_" + channel.str() + "_Timing_quantiles.eps").c_str());
    }
  }

  file -> Write();
  file -> Close();
  delete file;

  // cleaning memory
  if(holdPlots) return ;
  deleteFArray(MAX_HITS, f_ChargeBg, logfile);
  deleteFArray(MAX_HITS, f_ChargeSignal, logfile);
  deleteCanvasArray(nH2, cc_Harm1Mod, logfile);
  deleteCanvasArray(nH2, cc_Harm1Ph, logfile);
  deleteCanvasArray(nH2, cc_Harm1Re, logfile);
  deleteCanvasArray(nH2, cc_Harm1Im, logfile);
  deleteCanvasArray(nH2, cc_Harm1Polar, logfile);
  deleteCanvasArray(nH2, cc_Harm1Polar_bgSubtracted, logfile);
  deleteCanvasArray(nH2, cc_Harm1Cartesian, logfile);
  deleteCanvasArray(nH2, cc_Harm1Cartesian_bgSubtracted, logfile);
  deleteCanvasArray(nH2, cc_Harm1Mod_vs_T0, logfile);
  deleteCanvasArray(nH2, cc_Timing_vs_T0, logfile);
  deleteCanvasArray(nH2, cc_Charge_vs_T0, logfile);
  deleteCanvasArray(nH2, cc_Charge_vs_Timing, logfile);
  deleteCanvasArray(nH2, cc_T0, logfile);
  deleteCanvasArray(nH2, cc_Charge, logfile);
  deleteCanvasArray(nH2, cc_Charge_T0Cut_TimingCut, logfile);
  deleteCanvasArray(nH2, cc_Charge_vs_T0_bgSubtracted, logfile);
  deleteCanvasArray(nH2, cc_Charge_bgSubtracted, logfile);
  deleteCanvasArray(nH2, cc_Charge_bgSubtracted_filtered, logfile);
  deleteCanvasArray(nH2, cc_Timing, logfile);
  deleteCanvasArray(nH2, cc_Timing_quantiles, logfile);
  return ;
}

#endif
