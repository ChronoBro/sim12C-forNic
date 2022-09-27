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
Sty->SetPadTickY(1);
Sty->SetPadBottomMargin(.16);
Sty->SetPadTopMargin(.11);
Sty->SetPadLeftMargin(.19);
Sty->SetPadRightMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.6,"y");
Sty->SetTitleOffset(1.1,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();
 TCanvas canvas("thetaLab");


  TH2S frame("frame","",10,0,12,10,0,30000);
  frame.GetXaxis()->SetTitle("#theta_{lab} [deg]");
  frame.GetYaxis()->SetTitle("Counts");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH1I* histX = (TH1I*) file.Get("Be7/theta_7Be_2_el");
  histX->SetMarkerStyle(20);
  histX->SetMarkerSize(2.);
  histX->Draw("same PC");
  float Nx = histX->Integral();

  TLine ll;
  ll.SetLineStyle(2);
  ll.DrawLine(2.,0,2.,30000);
  ll.DrawLine(4.5,0,4.5,30000);

  TLatex tt;
  tt.DrawLatex(6,27000,"^{7}Be, J=7/2^{-}");

}
