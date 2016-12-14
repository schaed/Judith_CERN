#include <iostream>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;


void fitGaussian(
    TH1D* hist,
    double& mean,
    double& sigma,
    double& max,
    double& background,
    bool display)
{
  TF1* gauss = new TF1("g1", "[0]*exp(-0.5 * ((x-[1])/[2])**2) + [3]");
  gauss->SetParNames("Constant", "Mean", "Sigma", "Offset");

  int maxPos = 1;
  max = 0;
  int risePos = 1;
  double rise = 0;
  int fallPos = 1;
  double fall = 0;

  for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
  {
    const double value = hist->GetBinContent(bin);
    if (value > max)
    {
      max = value;
      maxPos = bin;
    }

    const double valuePrev = hist->GetBinContent(bin > 1 ? bin-1 : bin);
    if (value - valuePrev > rise) {
      rise = value - valuePrev;
      risePos = bin;
    }

    const double valueNext = hist->GetBinContent(
        bin < hist->GetNbinsX() ? bin+1 : bin);
    if (value - valueNext > fall) {
      fall = value - valueNext;
      fallPos = bin;
    }
  }

  mean = hist->GetXaxis()->GetBinCenter(maxPos);
  sigma = hist->GetXaxis()->GetBinUpEdge(fallPos) -
      hist->GetXaxis()->GetBinLowEdge(risePos);

  //cout << mean << endl;
  //cout << sigma << endl;
  //sigma = 20.;
  const unsigned int initialWidth = 5;

  Int_t bgBin1 = hist->FindBin(mean - initialWidth * sigma);
  Int_t bgBin2 = hist->FindBin(mean + initialWidth * sigma);
  if (bgBin1 < 1) bgBin1 = 1;
  if (bgBin2 > hist->GetNbinsX()) bgBin2 = hist->GetNbinsX();
  background = hist->GetBinContent(bgBin1) + hist->GetBinContent(bgBin2);
  background /= 2.0;

  max -= background;

  gauss->SetRange(mean - initialWidth * sigma, mean + initialWidth * sigma);
  gauss->SetParameters(max, mean, sigma, background);
   gauss->SetParLimits(2, 0, 10 * sigma);
  hist->Fit("g1", "QR0");

  max = gauss->GetParameter(0);
  mean = gauss->GetParameter(1);
  sigma = gauss->GetParameter(2);
  background = gauss->GetParameter(3);

  gauss->SetRange(mean - 5 * sigma, mean + 5 * sigma);
  gauss->SetParameters(max, mean, sigma, background);
  hist->Fit("g1", "QR0");

  max = gauss->GetParameter(0);
  mean = gauss->GetParameter(1);
  sigma = gauss->GetParameter(2);
  background = gauss->GetParameter(3);

  if (display)
  {
    TCanvas* can = new TCanvas();
    gauss->SetLineColor(46);
    gauss->SetLineWidth(1);
    hist->Draw();
    gauss->Draw("SAME");
    can->Update();
    can->WaitPrimitive();
  }

  delete gauss;
}


void fitGaussian(
    TH1D* hist,
    double& mean,
    double& sigma,
    bool display) {
  double max = 0;
  double background = 0;
  fitGaussian(hist, mean, sigma, max, background, display);
}


void fitBox(
    TH1D* hist,
    double& mean,
    double& sigma,
    double& max,
    double& background,
    double sensorWidth,
    bool display)
{

  TF1 *errorfunc = new TF1("errorfunc","[5]+[0]*(1+TMath::Erf([1]*([2]-x)))*(1-TMath::Erf([3]*([4]-x)))");
  errorfunc->SetParNames("Height", "FallSlope","Fall", "RiseSlope", "Rise","Offset");

  //fit a gaussian for first approximation
  fitGaussian(hist, mean, sigma, max, background, display);
  
  //first order approximations
  // double maxBin = hist->GetMaximumBin();
  double rise = mean - sensorWidth/2;
  double fall = mean + sensorWidth/2;
  double height = max/4;
  double riseSlope = height*10;
  double fallSlope = height*10;
  double offset = 0.2;
  // double maxBin = 0;
  // double rise = -1000;
  // double fall = 3000;
  // double maxHeight = 40;
  // double riseSlope = 400;
  // double fallSlope = 400;
  errorfunc->SetParLimits(0, height-abs(0.5*height), height+abs(0.5*height));
  errorfunc->SetParLimits(4, rise-abs(1*rise), rise+abs(1*rise));
  errorfunc->SetParLimits(2, fall-abs(1*fall), fall+abs(1*fall));
  errorfunc->SetParLimits(3, 0.01*riseSlope, 100*riseSlope);
  errorfunc->SetParLimits(1, 0.01*fallSlope, 100*fallSlope);
  errorfunc->SetParLimits(5, 0., 0.5);


  hist->Fit("errorfunc", "");

  //max = gauss->GetParameter(0);
  mean = ( errorfunc->GetParameter(2) + errorfunc->GetParameter(4) ) / 2 ;
  sigma = abs ( errorfunc->GetParameter(4) - errorfunc->GetParameter(2) );
  //background = gauss->GetParameter(3);

  cout << "Resolution min: " << 1/(sqrt(2.)*errorfunc->GetParameter(1)) << endl;
  cout << "Resolution max: " << 1/(sqrt(2.)*errorfunc->GetParameter(3)) << endl;

  if (display)
  {
    TCanvas* can = new TCanvas();
    errorfunc->SetLineColor(46);
    errorfunc->SetLineWidth(1);
    hist->Draw();
    errorfunc->Draw("SAME");
    can->Update();
    can->WaitPrimitive();
    cout << "Write name of canvas " << endl;
    string str;
    cin >> str;
    const char *name = (str+".root").c_str();
    can->SaveAs(name);

  }

  delete errorfunc;
}


void fitBox(
    TH1D* hist,
    double& mean,
    double& sigma,
    double sensorWidth,
    bool display) {
  double max = 0;
  double background = 0;
  fitBox(hist, mean, sigma, max, background, sensorWidth, display);
}


void fit() {

  // TFile *_file1 = new TFile("map_clusterType_4_val_pixX_30.000000_pixY_30.000000_resX_7.400000_resY_8.600000_gridX_1.000000_gridY_1.000000_integratorType_0.root");
  // TCanvas *_c1 = (TCanvas*)_file1->Get("cMap");
  // TH2D *_hist1 = (TH2D*)_c1->GetPrimitive("map");
  // _hist1->RebinX(5);
  // _hist1->RebinY(5);
  // TH1D* _h1x = _hist1->ProjectionX("_px",5,25);
  // TH1D* _h1y = _hist1->ProjectionY("_py",5,25);
  // _h1x->Scale(0.0036);

  TFile *f0 = new TFile("efficiencyQuick_projectionX_Irr_30.root");
  TCanvas *c0 = (TCanvas*)f0->Get("c4");
  TH1D *h0 = (TH1D*)c0->GetPrimitive("DUTPlane0TrackResidualHitFine_px");
   double mean=0., sigma=3., sensorWidth=100.;
   fitBox(h0, mean, sigma, sensorWidth,true);


  TFile *f1 = new TFile("efficiencyQuick_projectionY_Irr_30.root");
  TCanvas *c1 = (TCanvas*)f1->Get("c5");
  TH1D *h1 = (TH1D*)c1->GetPrimitive("DUTPlane0TrackResidualHitFine_py");
  mean=0.; 
  sigma=0.; 
  fitBox(h1, mean, sigma, sensorWidth,true);

  return;
}
