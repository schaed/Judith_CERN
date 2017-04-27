#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2D.h"
#include "TObject.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

Double_t ScaleX(Double_t x)
{
  Double_t v;
  v =  x - 35; // "linear scaling" function example
  return v;
}

Double_t ScaleY(Double_t y)
{
  Double_t v;
  v = y + 33; // "linear scaling" function example
  return v;
}

Double_t ScaleZ(Double_t z)
{
  Double_t v;
  v = 30 * z + 300; // "linear scaling" function example
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin()), // new Xmin
              Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Scale);
  return;
}

void ScaleYaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetYaxis(), Scale);
  return;
}

void ScaleZaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetZaxis(), Scale);
  return;
}

void drawline (TH2D * h1, Double_t xc, Double_t yc, TCanvas *c) {
 

  Double_t pixX = 25;
  Double_t pixY = 25;

  // h1->RebinX(5);
  // h1->RebinY(5);
  
 vector<TLine *> edges;
   
 TLine *ll  = NULL;


  // single pixel cluster
  ll = new TLine(xc +pixX/2.,yc +pixY/2.,xc +pixX/2.,yc -pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc -pixX/2., yc +pixY/2.,xc -pixX/2.,yc -pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc +pixX/2., yc +pixY/2.,xc -pixX/2., yc +pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc +pixX/2.,yc -pixY/2.,xc -pixX/2.,yc -pixY/2);
  edges.push_back(ll);
  

  // double pixels cluster
  ll = new TLine(xc +pixX/2., yc -pixY/2., xc +pixX/2., yc -3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc +pixX/2., yc -3.*pixY/2., xc -pixX/2., yc -3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc -pixX/2., yc -pixY/2., xc -pixX/2., yc -3.*pixY/2);
  edges.push_back(ll);
  

  // L-shaped pixels cluster
  ll = new TLine(xc -pixX/2., yc +pixY/2., xc -3.*pixX/2., yc +pixY/2);
  edges.push_back(ll);
   ll = new TLine(xc -3.*pixX/2.,yc +pixY/2., xc -3.*pixX/2.,yc -pixY/2);
   edges.push_back(ll);
   ll = new TLine(xc -pixX/2., yc -pixY/2., xc -3.*pixX/2., yc -pixY/2);
   edges.push_back(ll);
  

  // 2x2 pixels cluster
  ll = new TLine(xc -pixX/2., yc -3.*pixY/2., xc -3.*pixX/2., yc -3.*pixY/2);
  edges.push_back(ll);
  ll = new TLine(xc -3.*pixX/2.,yc -pixY/2., xc -3.*pixX/2., yc -3.*pixY/2);
  edges.push_back(ll);
  

  
  // // T-shaped pixels cluster
  // ll = new TLine(xc+3.*pixX/2., yc+3.*pixY/2., xc+5.*pixX/2.,yc+ 3.*pixY/2);
  // edges.push_back(ll);
  // ll = new TLine(xc+5.*pixX/2., yc+1.*pixY/2., xc+5.*pixX/2., yc+3.*pixY/2);
  // edges.push_back(ll);
  // ll = new TLine(xc+3.*pixX/2., yc+1.*pixY/2., xc+5.*pixX/2., yc+1.*pixY/2);
  // edges.push_back(ll);
  // ll = new TLine(xc+5.*pixX/2., yc-1.*pixY/2., xc+5.*pixX/2., yc+1.*pixY/2);
  // edges.push_back(ll);
  // ll = new TLine(xc+3*pixX/2., yc-1.*pixY/2., xc+5.*pixX/2., yc-1.*pixY/2);
  // edges.push_back(ll);
  //  ll = new TLine(xc+3.*pixX/2., yc-1.*pixY/2., xc+3.*pixX/2.,yc -3.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc+1.*pixX/2., yc-3.*pixY/2., xc+3.*pixX/2., yc-3.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc+1.*pixX/2., yc-3.*pixY/2., xc+1.*pixX/2., yc-1.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc-1.*pixX/2., yc-3.*pixY/2., xc+1.*pixX/2., yc-3.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc-1.*pixX/2., yc-3.*pixY/2., xc+-1.*pixX/2., yc-1.*pixY/2);
  //  edges.push_back(ll);
  // ll = new TLine(xc+1.*pixX/2., yc+3.*pixY/2.,xc+ 1.*pixX/2., yc+5.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc-1.*pixX/2., yc+3.*pixY/2.,xc -1.*pixX/2., yc+5.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc+3.*pixX/2., yc+3.*pixY/2.,xc+ 3.*pixX/2., yc+5.*pixY/2);
  //  edges.push_back(ll); 
  //  ll = new TLine(xc+1.*pixX/2., yc+5.*pixY/2.,xc+ 3.*pixX/2.,yc+ 5.*pixY/2);
  //  edges.push_back(ll);
  //  ll = new TLine(xc-1.*pixX/2., yc+5.*pixY/2., xc+1.*pixX/2., yc+5.*pixY/2);
  //  edges.push_back(ll);
  
  //DUTCS->RebinX(5);
  //DUTCS->RebinY(5);
    c->cd();
  h1->Draw("colz");

  for(unsigned int i=0; i<edges.size(); i++){
    edges[i] -> SetLineWidth(3);
    edges[i] -> Draw();
  }
  h1->SetTitle("");
  h1->GetXaxis()->SetTitle("X [#mum]");
  h1->GetXaxis()->SetTitleOffset(0.8);
  h1->GetYaxis()->SetTitle("Y [#mum]");
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetTitleSize(.05);
  h1->GetYaxis()->SetTitleSize(.05);
  //h1->GetZaxis()->SetRangeUser(200,1600);
  h1->GetZaxis()->SetRangeUser(1.,1.7);
  h1->GetXaxis()->SetRangeUser(-26,26.);
  h1->GetYaxis()->SetRangeUser(-26,26.);
  TLatex *lcp = new TLatex(0.65, 0.2, "25#mum Pitch" );
  lcp->SetNDC();
  lcp->SetTextSize(0.055);
  lcp->SetTextColor(1);
  lcp->Draw();

  TLatex *lcd = new TLatex(0.65, 0.27, "1e15 n_{eq}/cm^{2}" );
  lcd->SetNDC();
  lcd->SetTextSize(0.055);
  lcd->SetTextColor(1);
  lcd->Draw();

  
  TLatex *p1 = new TLatex(0.1, 0.2, "Pixel 1" );
  p1->SetNDC();
  p1->SetTextSize(0.055);
  p1->SetTextColor(1);
  p1->Draw();
  
  TLatex *p2 = new TLatex(0.1, 0.3, "Pixel 2" );
  p2->SetNDC();
  p2->SetTextSize(0.055);
  p2->SetTextColor(1);
  p2->Draw();
  
  TLatex *p3 = new TLatex(0.1, 0.4, "Pixel 3" );
  p3->SetNDC();
  p3->SetTextSize(0.055);
  p3->SetTextColor(1);
  p3->Draw();
  
  TLatex *p4 = new TLatex(0.1, 0.5, "Pixel 4" );
  p4->SetNDC();
  p4->SetTextSize(0.055);
  p4->SetTextColor(1);
  p4->Draw();
  
  TLatex *p5 = new TLatex(0.1, 0.4, "# of Clusters");
  p5->SetNDC();
  p5->SetNDC();
  p5->SetTextSize(0.050);
  p5->SetTextColor(1);
  p5->Draw();

  
  
  
  c->SetTickx();
  c->SetTicky();
  c->Update();
  return;

}


int plot_grid() {

  // TFile *f = new TFile("run_33234_tj_W3R13_M60_6V_3DRS_fft12_sync-analysis-result.root");
  // TH2D *DUTCS = (TH2D*)f->Get("Efficiency/DUTPlane0DUTClusterSize");
  // TH2D *DUTC = (TH2D*)f->Get("Efficiency/DUTPlane0DUTTime");
  // TH2D *DUTACC = (TH2D*)f->Get("Efficiency/DUTPlane0DUTCharge");
  // Double_t xc=-66;
  // Double_t yc=-21;
  
  // TCanvas *c1 = new TCanvas("c1","c1");
  // TCanvas *c2 = new TCanvas("c2","c2");
  // TCanvas *c3 = new TCanvas("c3","c3");
  // drawline(DUTCS,xc,yc,c1);
  // drawline(DUTC,xc,yc,c2);
  // drawline(DUTACC,xc,yc,c3);

  // TFile *f1 = new TFile("30efficiency.root");
  // TCanvas *can = (TCanvas*)f1->Get("c4");
  // TH2D *h1 = (TH2D*)can->GetPrimitive("DUTPlane0TrackResidualHitFine");
  // for (int i=50; i<100; i++) {
  //   h1->SetBinContent(91,i,0);
  //   h1->SetBinContent(i,91,0);
  //   //h1->SetBinContent(96,i,0);
  //   //h1->SetBinContent(i,96,0);
  // }
  
  
  
  
  // double xc=15;
  // double yc=15;
  // TCanvas *c4 = new TCanvas("c4","c4");
  // drawline(h1,xc,yc,c4);

  // TFile *f1 = new TFile("25charge.root");
  // TCanvas *can = (TCanvas*)f1->Get("Canvas_1");
  // TH2D *h1 = (TH2D*)can->GetPrimitive("DUTPlane0DUTCharge");
  // //h1->GetZaxis()->SetRangeUser(1.,1.7);
  
  
  // TCanvas *c4 = new TCanvas("c4","c4");
  // ScaleXaxis(h1, ScaleX);
  // ScaleYaxis(h1, ScaleY);
  // double xc=12.5;
  // double yc=12.5;
  // drawline(h1,xc,yc,c4);
  
  
  TFile *f1 = new TFile("25clustersize.root");
  TCanvas *can = (TCanvas*)f1->Get("Canvas_1");
  TH2D *h1 = (TH2D*)can->GetPrimitive("DUTPlane0DUTClusterSize");
  TCanvas *c4 = new TCanvas("c4","c4");
  ScaleXaxis(h1, ScaleX);
  ScaleYaxis(h1, ScaleY);
  double xc=12.5;
  double yc=12.5;
  drawline(h1,xc,yc,c4);

}
