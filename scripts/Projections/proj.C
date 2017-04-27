#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TObject.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

void average(TH1D* &proj, TH1D* &proj4){
  int bins = proj->GetNbinsX();
  double v1, v2, e1, e2,v,e,w1,w2;
  for (int i=0; i<bins; i++){
    v1 = proj->GetBinContent(i+1);
    v2 = proj4->GetBinContent(i+1);
    e1 = proj->GetBinError(i+1);
    e2 = proj4->GetBinError(i+1);
    w1 = 1/(e1*e1);
    w2 = 1/(e2*e2);
    v = (v1*w1+v2*w2)/(w1+w2);
    e = 1/sqrt(w1+w2);
    proj->SetBinContent(i+1,v);
    proj->SetBinError(i+1,e);
  }
}

void doproj(TH1D* proj, TH1D* proj2, double &mean, double &std, double &sist, string pos) {
  int bins = proj->GetNbinsX();
  TCanvas *can1 = new TCanvas("can1","can1",900,600);
  proj->Draw();
  can1->Update();
  can1->WaitPrimitive();
  for (int i=0; i<bins; i++) {
    proj->SetBinContent(i+1,proj->GetBinContent(i+2));
  }
  TCanvas *can2 = new TCanvas("can2","can2",900,600);
  proj->Draw();
  can2->Update();
  can2->WaitPrimitive();
  int RebinF=2;
  proj->Rebin(RebinF);
  proj2->Rebin(RebinF);
  proj->Scale(0.5);
  proj2->Scale(0.5);
  TCanvas *c1 = new TCanvas("c1","c1",900,600);
  proj->Draw();
  c1->Update();
  c1->WaitPrimitive();
   vector<double> values,sist_errors,stat_errors_down, stat_errors_up;
  for (int i=0; i<bins; i++) {
    if (proj->GetBinContent(i+1) > 1 ) values.push_back(1.);
    else values.push_back(proj->GetBinContent(i+1));
    sist_errors.push_back(proj2->GetBinContent(i+1));
    stat_errors_down.push_back(proj->GetBinError(i+1));
    if (proj->GetBinContent(i+1)+proj->GetBinError(i+1)>1) stat_errors_up.push_back(abs(1.-values[i]));
    else stat_errors_up.push_back(proj->GetBinError(i+1));
  }
  double vx=proj->GetXaxis()->GetBinLowEdge(1);
  double vx1=proj->GetXaxis()->GetBinLowEdge(2);
  double dx = vx1-vx;
  double *valuesmid = new double[values.size()];
  double *x = new double[values.size()];
  double *errorx = new double[values.size()];
  double *errorystup = new double[values.size()];
  double *errorystdown = new double[values.size()];
  double *errorysi = new double [values.size()];
  for (int i=0; i<values.size(); i++){
    valuesmid[i] = values.at(i);
    x[i] = vx+i*dx;
    errorx[i] = dx/2.;
    errorystdown[i] = stat_errors_down.at(i);
    errorystup[i] = stat_errors_up.at(i);
    errorysi[i] = sqrt(pow(abs(values.at(i)-sist_errors.at(i)),2)+pow(stat_errors_down.at(i),2)+pow(values.at(i)*0.002,2));
  }
  
  for (int i=5; i<10; i++) {
    mean += valuesmid[i];
    std += errorystdown[i];
    sist += errorysi[i] -errorystdown[i];
  } 
  mean /=5;
  std /=5;
  sist /=5;
  
  TGraphAsymmErrors *g1 = new TGraphAsymmErrors(values.size(),x,valuesmid,errorx,errorx,errorystdown,errorystup);
  TGraphErrors *g2 = new TGraphErrors(values.size(),x,valuesmid,errorx,errorysi);
  g1->SetLineColor(4);
  g2->SetLineColor(2);
  g2->SetTitle("");
  g2->GetXaxis()->SetRangeUser(-20,70);
  g2->GetYaxis()->SetTitle("Efficiency");
  const char* title = (pos+" position [#mum]").c_str();
  g2->GetXaxis()->SetTitle(title);
  TCanvas *c2 = new TCanvas;
  g2->Draw("AP");
  g1->Draw("SAMEP");
  c2->Update();
  c2->WaitPrimitive();
  cout << "Write name of canvas " << endl;
  string str;
  cin >> str;
  const char *name = (str+".root").c_str();
  c2->SaveAs(name);
  delete g1;
  delete g2;
  delete valuesmid;
  delete x;
  delete errorx;
  delete errorystdown;
  delete errorystup;
  delete errorysi;
  delete c2;
  return;
}

int proj() {

  TFile *f = new TFile("../CorrEfficiency/efficiencyQuick_projectionX_corr_UnIrr_50.root");
  TCanvas *c1 = (TCanvas*)f->Get("c6");
  TH1D *proj = (TH1D*)c1->GetPrimitive("DUTPlane0TrackResidualHitFine_px");
  TFile *f2 = new TFile("../CorrEfficiency/efficiencyQuick_projectionX_corr_UnIrr_50_err.root");
  TCanvas *c3 = (TCanvas*)f2->Get("c6");
  TH1D *proj2 = (TH1D*)c3->GetPrimitive("DUTPlane0TrackResidualHitFine_px");
  // TFile *f4 = new TFile("../CorrEfficiency/");
  // TCanvas *c4 = (TCanvas*)f4->Get("c6");
  // TH1D *proj4 = (TH1D*)c4->GetPrimitive("DUTPlane0TrackResidualHitFine_px");
  double mean = 0.;
  double std = 0.;
  double sist = 0.;
  
  //average(proj,proj4);
  doproj(proj, proj2, mean,std,sist, "X");

  cout << "Mean value of efficiency: " << mean << " +- " << std <<" (stat) +- " << sist << " (sist) " <<  endl;

  f = TFile::Open("../CorrEfficiency/efficiencyQuick_projectionY_corr_UnIrr_50.root");
  c1 = (TCanvas*)f->Get("c7");
  proj = (TH1D*)c1->GetPrimitive("DUTPlane0TrackResidualHitFine_py");
  f2 = TFile::Open("../CorrEfficiency/efficiencyQuick_projectionY_corr_UnIrr_50_err.root");
  c3 = (TCanvas*)f2->Get("c7");
  proj2 = (TH1D*)c3->GetPrimitive("DUTPlane0TrackResidualHitFine_py");
  // f4 = new TFile("../CorrEfficiency/");
  // c4 = (TCanvas*)f4->Get("c7");
  // proj4 = (TH1D*)c4->GetPrimitive("DUTPlane0TrackResidualHitFine_py");
  double mean1 = 0.;
  std = 0.;
  sist=0.;

  //average(proj,proj4);
  doproj(proj, proj2,mean1,std,sist,"Y");
  
  cout << "Mean value of efficiency: " << mean1 << " +- " << std <<" (stat) +- " << sist << " (sist) " <<  endl;

    
  return 0;
}
