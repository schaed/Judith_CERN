{

  TCanvas *cc = new TCanvas("cc", "cc");

  TFile *file = TFile::Open("hCalc.root");
  h1 -> SetTitle("calculation");
  h1 -> SetLineColor(1);
  h1 -> SetLineWidth(2);
  h1 -> DrawNormalized();

  TFile *file = TFile::Open("hMC.root");
  h1 -> SetTitle("Monte Carlo");
  h1 -> SetLineColor(2);
  h1 -> SetLineWidth(2);
  h1 -> DrawNormalized("same");

  TLegend *leg = cc -> BuildLegend(0.1178161,0.7468354,0.3821839,0.8734177);
  leg -> SetFillColor(0);
  leg -> SetFillStyle(0);
  leg -> SetLineColor(0);
  leg -> Draw();

  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);

  cc -> Modified();
  cc -> Update();

  cc -> SaveAs("compare_calc_vs_MC.png");
  cc -> SaveAs("compare_calc_vs_MC.pdf");
  cc -> SaveAs("compare_calc_vs_MC.eps");
  cc -> SaveAs("compare_calc_vs_MC.root");
}
