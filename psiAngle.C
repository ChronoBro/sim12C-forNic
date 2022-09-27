{
  gROOT->Reset();

TStyle * Sty = new TStyle("MyStyle","MyStyle");
Sty->SetOptTitle(0);
Sty->SetOptStat(0);
//Sty->SetPalette(8,0);
Sty->SetCanvasColor(10);
Sty->SetCanvasBorderMode(0);
Sty->SetFrameLineWidth(3);
Sty->SetFrameFillColor(10);
Sty->SetPadColor(10);
Sty->SetPadTickX(1);
Sty->SetPadTickY(0);
Sty->SetPadBottomMargin(0.20);
Sty->SetPadTopMargin(0.05);
Sty->SetPadLeftMargin(0.2);
Sty->SetPadRightMargin(0.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kBlack);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.09,"xyz");
Sty->SetTitleOffset(1.1,"y");
Sty->SetTitleOffset(.9,"x");
Sty->SetTitleOffset(1.5,"z");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.04,"xyz");
Sty->SetNdivisions(5,"xyz");
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("psiAngle_7Be72","");

  double x[7]={-3,-2,-1,0,1,2,3};
  double ymid[7] = {0.3144,.1190,.0484,.0361,0.0484,.1190,.3144};
  double yhi[7] = {.2867,.1297,.0620,0.0433,0.0620,.1297,.2867};
  double ylow[7] = {.2505,.1035,.0935,.102,.0935,.1035,.2505};
  double ylowErr[7]={.002,.0036,.0022,.0119,.0025,.0036,.002};
  double xerr[7]={0.};

  TH2S frame("frame","",10,-4,4,10,0,.35);
  frame.GetXaxis()->SetTitle("\\ m");
  frame.GetYaxis()->SetTitle("\\rho^{\\ell}_{m,m}");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.GetXaxis()->SetNdivisions(8);
  frame.Draw();

  TGraph gmid(7,x,ymid);
  gmid.SetMarkerStyle(20);
  gmid.SetMarkerColor(2);
  gmid.SetLineColor(2);
  gmid.SetLineWidth(2);
  gmid.SetMarkerSize(1.5);
  gmid.Draw("PL");


  TGraph ghi(7,x,yhi);
  ghi.SetMarkerStyle(4);
  ghi.SetMarkerColor(4);
  ghi.SetLineColor(4);
  ghi.SetLineWidth(2);
  ghi.SetMarkerSize(1.5);
  ghi.SetLineStyle(2);
  ghi.Draw("PL");


  TGraph glow(7,x,ylow);
  glow.SetMarkerStyle(23);
  glow.SetMarkerColor(3);
  glow.SetLineColor(3);
  glow.SetLineWidth(2);
  glow.SetLineStyle(9);
  glow.SetMarkerSize(1.5);
  glow.Draw("PL");

  TGraphErrors glow2(7,x,ylow,xerr,ylowErr);
  glow2.SetMarkerStyle(23);
  glow2.SetMarkerColor(3);
  glow2.SetLineColor(3);
  glow2.SetLineWidth(2);
  glow2.SetMarkerSize(1.5);
  glow2.Draw("P");

  TLatex tt;
  tt.SetTextSize(.06);
  tt.SetTextColor(2);
  tt.DrawLatex(-2.8,.3,"2^{#circ}<#theta^{*}<4.5^{#circ}");
  tt.SetTextColor(3);
  tt.DrawLatex(-.3,.12,"#theta^{*}<2^{#circ}");
  tt.SetTextColor(4);
  tt.DrawLatex(-.4,.055,"#theta^{*}>4.5^{#circ}");


}
