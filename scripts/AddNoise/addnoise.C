#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include <string>
#include <sstream>

using namespace std;


string NumberToString ( int Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }


void addnoise() {
  
  TFile *file0 = new TFile("50signal_1.root");
  TCanvas *c0 = (TCanvas*)file0->Get("Canvas_1"); 
  TH1D *signal = (TH1D*)c0->GetPrimitive("DUTPlane0ChargeDUT");  //DUTPlane0ChargeDUT or hh
  TFile *file1 = new TFile("noiseM129_irr.root");
  TCanvas *c1 = (TCanvas*)file1->Get("prov0");
  TH1F *n = (TH1F*)c1->GetPrimitive("hrmsch0");
  int nbins = n->GetXaxis()->GetNbins();
  TH1F *noise = new TH1F("noise","title",nbins,-1000,9000);
  for (int i=1;i<=nbins;i++) {
      double y = n->GetBinContent(i);
      double x = n->GetXaxis()->GetBinCenter(i);
      double xnew = x*66203;  
      noise->Fill(xnew,y);
  }
  TCanvas *can =  new TCanvas("can","can",900,600);
    noise->Rebin(25);
  // noise->Draw();
  // can->Update();
  // can->WaitPrimitive();
  noise->Scale(0.0025);
  noise->SetLineColor(1);
  noise->SetLineStyle(2);
  noise->SetLineWidth(2);
  signal->SetLineWidth(2);
  signal->Rebin(20);
  TF1 *gaus = new TF1("gaus","gaus",-0.005,0.005);
  TF1 *landau = new TF1("landau","landau",0.,0.2);
  noise->Fit("gaus");
  signal->Fit("landau");
  can->cd();
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(noise,"noise");
  leg->AddEntry(signal,"hits");
  TLatex *lcp = new TLatex(0.2, 0.2, "50#mum Pitch" );
  lcp->SetNDC();
  lcp->SetTextSize(0.055);
  lcp->SetTextColor(1);
  TLatex *lcd = new TLatex(0.2, 0.27, "1e15 n_{eq}/cm^{2}" );
  lcd->SetNDC();
  lcd->SetTextSize(0.055);
  lcd->SetTextColor(1);
  string par1 = NumberToString(int(landau->GetParameter(1)));
  string str = "Signal mpv: "+par1+ " e";
  const char *str1 = str.c_str();
  TLatex *landst = new TLatex(0.2, 0.2, str1 );
  landst->SetNDC();
  landst->SetTextSize(0.05);
  landst->SetTextColor(4);
  par1 = NumberToString(int(gaus->GetParameter(2)));
  str = "Noise std dev: "+par1+ " e";
  const char *str2 = str.c_str();
  TLatex *gausst = new TLatex(0.2, 0.27, str2 );
  gausst->SetNDC();
  gausst->SetTextSize(0.05);
  gausst->SetTextColor(1);

  noise->SetTitle("");
  noise->GetXaxis()->SetTitle("Signal Size [e]");
  noise->GetYaxis()->SetTitle("Hits");
  noise->SetStats(0);
  signal->SetStats(0);
  noise->Draw();
  signal->Draw("same");
  leg->Draw();
  lcd->Draw();
  lcp->Draw();
  landst->Draw();
  gausst->Draw();
  can->SetTickx();
  can->SetTicky();
  can->Update();

   return;

}
