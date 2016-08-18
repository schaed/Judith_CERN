#ifndef MAKEMAP_HH
#define MAKEMAP_HH

#include "TF2.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/GSLMCIntegrator.h"
#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"

#define CLUSTERTYPE_SINGLE 1
#define CLUSTERTYPE_DOUBLE 2
#define CLUSTERTYPE_LSHAPE 3
#define CLUSTERTYPE_2X2 4
#define INTEGRATORTYPE_ADAPTIVE 0
#define INTEGRATORTYPE_GSLMC 1
#define INTEGRATOR_TOLERANCE 0.0001
#define DEFAULT_DOUBLE -1.
#define DEFAULT_UINT 0
#define CANVASSIZE 500

void style(){

  gStyle -> SetPaperSize(20, 20);
  gStyle -> SetPadTopMargin(0.15);
  gStyle -> SetPadRightMargin(0.15);
  gStyle -> SetPadBottomMargin(0.15);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetOptTitle(0);
  gStyle -> SetOptFit(1111);
  gStyle->SetTitleW(0.6);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetPalette(1);

}

int defineSurface(double &xMin, 
		  double &xMax, 
		  double &yMin, 
		  double &yMax,
		  const double pixX,
		  const double pixY,
		  const unsigned int clusterType){

  xMin = -pixX;
  yMin = -pixY;
  if(clusterType == CLUSTERTYPE_SINGLE){
    xMax = pixX;
    yMax = pixY;
  }
  else if(clusterType == CLUSTERTYPE_DOUBLE){
    xMax = pixX;
    yMax = 2.*pixY;
  }
  else if(clusterType == CLUSTERTYPE_LSHAPE || clusterType == CLUSTERTYPE_2X2){
    xMax = 2.*pixX;
    yMax = 2.*pixY;
  }
  else{
    cout << " - ERROR!!! - undefined cluster type clusterType = " << clusterType << endl;
    return 1;
  }

  return 0;
}

void defineBinning(unsigned int &nX, 
		   unsigned int &nY,
		   const double xMax, 
		   const double xMin, 
		   const double gridX,
		   const double yMax,
		   const double yMin, 
		   const double gridY){

  nX = (xMax - xMin) / gridX;
  nY = (yMax - yMin) / gridY;

  double diffX = nX - (xMax - xMin) / gridX;
  double diffY = nY - (yMax - yMin) / gridY;
  if(diffX != 0. || diffY != 0.){
    cout << " - WARNING!!! - map may contain Moiree patterns" << endl;
  }

  return ;  
}

TCanvas *getCanvas(const unsigned int clusterType,
		   const char *flag){
  TCanvas *cc = NULL;
  if(clusterType == CLUSTERTYPE_SINGLE){
    cc = new TCanvas(flag, flag, 0, 0, CANVASSIZE, CANVASSIZE);
  }
  else if(clusterType == CLUSTERTYPE_DOUBLE){
    cc = new TCanvas(flag, flag, 0, 0, CANVASSIZE, 2*CANVASSIZE);
  }
  else if(clusterType == CLUSTERTYPE_LSHAPE || clusterType == CLUSTERTYPE_2X2){
    cc = new TCanvas(flag, flag, 0, 0, 2*CANVASSIZE, 2*CANVASSIZE);
  }
  else{
    cout << " - ERROR!!! - undefined cluster type clusterType = " << clusterType << endl;
    return NULL;
  }
  return cc;
}

int getIntegralRanges(const double pixX,
		      const double pixY,
		      vector<double *> &integralRangeMin,
		      vector<double *> &integralRangeMax,
		      const unsigned int clusterType){

  // single pixel cluster
  double *integralRangeMin1 = new double[2];
  integralRangeMin1[0] = -pixX/2.;
  integralRangeMin1[1] = -pixY/2.;
  double *integralRangeMax1 = new double[2];
  integralRangeMax1[0] = pixX/2.;
  integralRangeMax1[1] = pixY/2.;
  integralRangeMin.push_back(integralRangeMin1);
  integralRangeMax.push_back(integralRangeMax1);
  if(clusterType == CLUSTERTYPE_SINGLE) return 0;

  // double pixels cluster
  double *integralRangeMin2 = new double[2];
  integralRangeMin2[0] = -pixX/2.;
  integralRangeMin2[1] = pixY/2.;
  double *integralRangeMax2 = new double[2];
  integralRangeMax2[0] = pixX/2.;
  integralRangeMax2[1] = 3.*pixY/2.;
  integralRangeMin.push_back(integralRangeMin2);
  integralRangeMax.push_back(integralRangeMax2);
  if(clusterType == CLUSTERTYPE_DOUBLE) return 0;

  // L-shaped pixels cluster
  double *integralRangeMin3 = new double[2];
  integralRangeMin3[0] = pixX/2.;
  integralRangeMin3[1] = -pixY/2.;
  double *integralRangeMax3 = new double[2];
  integralRangeMax3[0] = 3.*pixX/2.;
  integralRangeMax3[1] = pixY/2.;
  integralRangeMin.push_back(integralRangeMin3);
  integralRangeMax.push_back(integralRangeMax3);
  if(clusterType == CLUSTERTYPE_LSHAPE) return 0;

  // 2x2 pixels cluster
  double *integralRangeMin4 = new double[2];
  integralRangeMin4[0] = pixX/2.;
  integralRangeMin4[1] = pixY/2.;
  double *integralRangeMax4 = new double[2];
  integralRangeMax4[0] = 3.*pixX/2.;
  integralRangeMax4[1] = 3.*pixY/2.;
  integralRangeMin.push_back(integralRangeMin4);
  integralRangeMax.push_back(integralRangeMax4);
  if(clusterType == CLUSTERTYPE_2X2) return 0;

  // undefined cluster type
  cout << " - ERROR!!! - undefined cluster type clusterType = " << clusterType << endl;
  return 1;
}

double computeIntegral(TF2 *pdfArg,
		       double *integralRangeMin,
		       double *integralRangeMax,
		       const unsigned int integratorType){

  ROOT::Math::WrappedMultiTF1 wf1(*pdfArg);
  if(integratorType == INTEGRATORTYPE_ADAPTIVE){
    ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetFunction(wf1);
    ig.SetRelTolerance(INTEGRATOR_TOLERANCE);
    return ig.Integral(integralRangeMin, integralRangeMax);
  }
  else if(integratorType == INTEGRATORTYPE_GSLMC){
    ROOT::Math::GSLMCIntegrator ig; // much slower
    ig.SetFunction(wf1);
    ig.SetRelTolerance(INTEGRATOR_TOLERANCE);
    return ig.Integral(integralRangeMin, integralRangeMax);
  }
  else{
    cout << " - ERROR!!! - integrator not implemented for integratorType = " << integratorType << endl;
    return 0.;
  }

}
  
int buildMap(const double pixX, 
	     const double pixY, 
	     const double resX, 
	     const double resY, 
	     const double resErrX, 
	     const double resErrY, 
	     const unsigned int nX,
	     const unsigned int nY,
	     const double xMin,
	     const double xMax,
	     const double yMin,
	     const double yMax,
	     const unsigned int clusterType,
	     const unsigned int integratorType, 
	     TF2 *pdfArg, 
	     TF2 *pdfArgDerivativeX, 
	     TF2 *pdfArgDerivativeY, 
	     TH2F *map,
	     TH2F *mapErr){

  // calculating grid steps
  double dx = (xMax - xMin) / nX;
  double dy = (yMax - yMin) / nY;

  // calculating integration ranges
  vector<double *> integralRangeMin;
  vector<double *> integralRangeMax;
  if(getIntegralRanges(pixX,
		       pixY,
		       integralRangeMin,
		       integralRangeMax,
		       clusterType)){
    cout << " - ERROR!!! - cannot get integral ranges" << endl;
    return 1;
  }
  for(unsigned int iPdf=0; iPdf<integralRangeMin.size(); iPdf++){
    cout << " - pixel " << iPdf+1 << " x range = " << integralRangeMin[iPdf][0] << " : " << integralRangeMax[iPdf][0] << endl;
    cout << " - pixel " << iPdf+1 << " y range = " << integralRangeMin[iPdf][1] << " : " << integralRangeMax[iPdf][1] << endl;
  }

  // setting common parameters
  pdfArg -> SetParameter(1, resX);
  pdfArg -> SetParameter(3, resY);
  pdfArgDerivativeX -> SetParameter(1, resX);
  pdfArgDerivativeX -> SetParameter(3, resY);
  pdfArgDerivativeY -> SetParameter(1, resX);
  pdfArgDerivativeY -> SetParameter(3, resY);

  // looping over sensor surface
  for(unsigned int ix=0; ix<=nX; ix++){
    for(unsigned int iy=0; iy<=nY; iy++){
      double x0 = xMin+ ix * dx;
      double y0 = yMin+ iy * dy;
      double pdf = 0.;
      double pdfErrComponentX = 0.;
      double pdfErrComponentY = 0.;
      for(unsigned int iPdf=0; iPdf<integralRangeMin.size(); iPdf++){
	pdfArg -> SetParameter(0, x0);
	pdfArg -> SetParameter(2, y0);
	pdfArgDerivativeX -> SetParameter(0, x0);
	pdfArgDerivativeX -> SetParameter(2, y0);
	pdfArgDerivativeY -> SetParameter(0, x0);
	pdfArgDerivativeY -> SetParameter(2, y0);
      	pdf += computeIntegral(pdfArg,
			       integralRangeMin[iPdf],
			       integralRangeMax[iPdf],
			       integratorType);
	pdfErrComponentX += computeIntegral(pdfArgDerivativeX,
					    integralRangeMin[iPdf],
					    integralRangeMax[iPdf],
					    integratorType);
	pdfErrComponentY += computeIntegral(pdfArgDerivativeY,
					    integralRangeMin[iPdf],
					    integralRangeMax[iPdf],
					    integratorType);
      }
      double pdfErr = sqrt(pdfErrComponentX * pdfErrComponentX * resErrX * resErrX + pdfErrComponentY * pdfErrComponentY * resErrY * resErrY);
      map -> Fill(x0, y0, pdf);
      mapErr -> Fill(x0, y0, pdfErr);
    }
  }

  return 0;
}

int getPixelEdges(vector<TLine *> &edges,
		  const unsigned int clusterType,
		  const double pixX,
		  const double pixY){

  TLine *ll = NULL;

  // single pixel cluster
  ll = new TLine(-pixX/2., -pixY/2., -pixX/2., pixY/2);
  edges.push_back(ll);
  ll = new TLine(pixX/2., -pixY/2., pixX/2., pixY/2);
  edges.push_back(ll);
  ll = new TLine(-pixX/2., -pixY/2., pixX/2., -pixY/2);
  edges.push_back(ll);
  ll = new TLine(-pixX/2., pixY/2., pixX/2., pixY/2);
  edges.push_back(ll);
  if(clusterType == CLUSTERTYPE_SINGLE) return 0;

  // double pixels cluster
  ll = new TLine(-pixX/2., pixY/2., -pixX/2., 3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(-pixX/2., 3.*pixY/2., pixX/2., 3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(pixX/2., pixY/2., pixX/2., 3.*pixY/2);
  edges.push_back(ll);
  if(clusterType == CLUSTERTYPE_DOUBLE) return 0;

  // L-shaped pixels cluster
  ll = new TLine(pixX/2., -pixY/2., 3.*pixX/2., -pixY/2);
  edges.push_back(ll);
  ll = new TLine(3.*pixX/2., -pixY/2., 3.*pixX/2., pixY/2);
  edges.push_back(ll);
  ll = new TLine(pixX/2., pixY/2., 3.*pixX/2., pixY/2);
  edges.push_back(ll);
  if(clusterType == CLUSTERTYPE_LSHAPE) return 0;

  // 2x2 pixels cluster
  ll = new TLine(pixX/2., 3.*pixY/2., 3.*pixX/2., 3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(3.*pixX/2., pixY/2., 3.*pixX/2., 3.*pixY/2);
  edges.push_back(ll);
  if(clusterType == CLUSTERTYPE_2X2) return 0;

  // undefined cluster type
  cout << " - ERROR!!! - undefined cluster type clusterType = " << clusterType << endl;
  return 1;
}

int getAxes(vector<TLine *> &axes,
	    const unsigned int clusterType,
	    const double pixX,
	    const double pixY,
	    const double gridX,
	    const double gridY){
  
  TLine *ax = NULL;
  TLine *ay = NULL;
  if(clusterType == CLUSTERTYPE_SINGLE){
    ax = new TLine(-pixX, 0, pixX-gridX/2., 0);
    ay = new TLine(0, -pixY, 0, pixY-gridY/2.);
  }
  else if(clusterType == CLUSTERTYPE_DOUBLE){
    ay = new TLine(-pixX, 0, pixX-gridX/2., 0);
    ax = new TLine(0, -pixY, 0, 2.*pixY-gridY/2.);
  }
  else if(clusterType == CLUSTERTYPE_LSHAPE || clusterType == CLUSTERTYPE_2X2){
    ay = new TLine(-pixX, 0, 2.*pixX-gridX/2., 0);
    ax = new TLine(0, -pixY, 0, 2.*pixY-gridY/2.);
  }
  else{
    cout << " - ERROR!!! - undefined cluster type clusterType = " << clusterType << endl;
    return 1;
  }
  axes.push_back(ax);
  axes.push_back(ay);

  return 0;
}

#endif
