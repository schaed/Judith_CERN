#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"


int clustersize() {

  TFile *f1 = new TFile("1dclustersize25_1.root");
  TCanvas *c1 = (TCanvas*)f1->Get("Canvas_1");
  TH1D *cl1 = (TH1D*)c1->GetPrimitive("DUTPlane0SizeDUT"); //DUTPlane0SizeDUT
  TFile *f2 = new TFile("1dclustersize30_1.root");
  TCanvas *c2 = (TCanvas*)f2->Get("Canvas_1");
  TH1D *cl2 = (TH1D*)c2->GetPrimitive("DUTPlane0SizeDUT");
  TFile *f3 = new TFile("1dclustersize50_unirr_1.root");
  TCanvas *c3 = (TCanvas*)f3->Get("Canvas_1");
  TH1D *cl3 = (TH1D*)c3->GetPrimitive("DUTPlane0SizeDUT");

  cl1->Sumw2();
  cl2->Sumw2();
  cl3->Sumw2();
  cl1->SetLineColor(2);
  cl1->SetLineWidth(3);
  cl2->SetLineColor(3);
  cl2->SetLineWidth(3);
  cl3->SetLineColor(4);
  cl3->SetLineWidth(3);
  cl1->Scale(1/1222.);
  cl2->Scale(1/2245.);
  cl3->Scale(1/3064.);
  cl3->GetXaxis()->SetTitle("Cluster Size");
  cl3->GetYaxis()->SetTitle("Events");
  cl3->GetXaxis()->SetTitleSize(0.055);
  cl3->GetYaxis()->SetTitleSize(0.055);
  cl3->GetXaxis()->SetTitleOffset(0.8);
  cl1->SetMarkerColor(cl1->GetLineColor());
  cl2->SetMarkerColor(cl2->GetLineColor());
  cl3->SetMarkerColor(cl3->GetLineColor());
  cl1->SetMarkerStyle(20);
  cl2->SetMarkerStyle(22);
  cl3->SetMarkerStyle(23);
  cl3->GetYaxis()->SetRangeUser(0.005,1.);
  cl3->GetYaxis()->SetTitleOffset(0.8);
  cl3->GetYaxis()->SetLabelSize(0.045);
  cl3->GetXaxis()->SetLabelSize(0.045);
  cl1->SetStats(0);
  cl2->SetStats(0);
  cl3->SetStats(0);
  cl3->SetTitle("");

  TLegend *leg=new TLegend(0.2,0.2,0.6,0.4);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(cl3,"50x50 #mum - Unirradiated");
  leg->AddEntry(cl2,"30x30 #mum - 1e15 n_{eq}/cm^{2}");
  leg->AddEntry(cl1,"25x25 #mum - 1e15 n_{eq}/cm^{2}");

  TLegend *leg1=new TLegend(0.3,0.3,0.7,0.5);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->AddEntry(cl3,"mean: 1.374");
  leg1->AddEntry(cl2,"mean: 1.345");
  leg1->AddEntry(cl1,"mean: 1.059");


  TCanvas *c0 = new TCanvas("c0","c0",900,600);
  cl3->Draw("E");
  cl1->Draw("SAME E");
  cl2->Draw("SAME E");
  leg->Draw("");
  leg1->Draw("");
  c0->SetLogy();
  c0->SetTickx();
  c0->SetTicky();
  c0->Update();

  return 0;
}
