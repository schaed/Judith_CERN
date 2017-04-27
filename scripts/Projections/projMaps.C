#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

void projMaps() {

  TFile *_file0 = new TFile("maps/map_clusterType_4_val_pixX_30.000000_pixY_29.000000_resX_8.500000_resY_10.500000_gridX_1.000000_gridY_1.000000_integratorType_0.root");
  TCanvas *_c0 = (TCanvas*)_file0->Get("cMap");
  TH2D *_hist0 = (TH2D*)_c0->GetPrimitive("map");
  _hist0->RebinX(5);
  _hist0->RebinY(5);
  TH1D* _h0x = _hist0->ProjectionX("_px",5,10);
  TH1D* _h0y = _hist0->ProjectionX("_py",5,10);
  
  TFile *_file1 = new TFile("maps/map_clusterType_4_val_pixX_30.000000_pixY_29.000000_resX_8.500000_resY_10.500000_gridX_1.000000_gridY_1.000000_integratorType_0.root");
  TCanvas *_c1 = (TCanvas*)_file1->Get("cMap");
  TH2D *_hist1 = (TH2D*)_c1->GetPrimitive("map"); 
  _hist1->RebinX(5);
  _hist1->RebinY(5);
  TH1D* _h1x = _hist1->ProjectionX("_px",5,10);
  TH1D* _h1y = _hist1->ProjectionY("_py",5,10);
  
  
  TFile *_file2 = new TFile("../CorrEfficiency/fitx.root");
  TCanvas *_c2 = (TCanvas*)_file2->Get("c1_n2");
  TH1D *_h2x = (TH1D*)_c2->GetPrimitive("DUTPlane0TrackResidualHitFine_px"); 

  TFile *_file3 = new TFile("../CorrEfficiency/fity.root");
  TCanvas *_c3 = (TCanvas*)_file3->Get("c1_n4");
  TH1D *_h2y = (TH1D*)_c3->GetPrimitive("DUTPlane0TrackResidualHitFine_py"); 
 
  TCanvas *canX = new TCanvas("canX","canX");
  canX->cd();
  _h0x->Scale(0.007);
  _h1x->Scale(0.007);
  _h0x->SetLineColor(3);
  _h0x->SetLineWidth(4);
  _h1x->SetLineColor(4);
  _h1x->SetLineWidth(4);
  _h0x->SetTitle("Projection X");
  _h0x->Draw("LP");
  _h1x->Draw("SAMELP");
  _h2x->Draw("SAMELP");
  canX->Update();
  canX->SaveAs("projXmap.root");
  
  TCanvas *canY = new TCanvas("canY","canY");
  canY->cd();
  _h0y->Scale(0.007);
  _h1y->Scale(0.007);
  _h0y->SetLineColor(3);
  _h0y->SetLineWidth(4);
  _h1y->SetLineColor(4);
  _h1y->SetLineWidth(4);
  _h0y->SetTitle("Projection Y");
  _h0y->Draw("LP");
  _h1y->Draw("SAMELP");
  _h2y->Draw("SAMELP");
  canY->Update();
  canY->SaveAs("projYmap.root");
  return;
  

}
