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
Sty->SetPadTopMargin(.05);
Sty->SetPadLeftMargin(.16);
Sty->SetPadRightMargin(.05);
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
Sty->SetTitleOffset(1.1,"y");
Sty->SetTitleOffset(1.1,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

  TCanvas canvas("plm","",900,600);
 double overlap = .04;
TPad * pad1 = new TPad("pad1","",0.,0.,.5+overlap,1.);
 TPad * pad2 = new TPad("pad2","",.5-overlap,0,1.,1.);


pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);


pad1->Draw();
pad2->Draw();



 pad1->cd();
  TH2S frame1("frame1","",10,-1,1,10,-4,3);
  frame1.SetStats(kFALSE);
  frame1.GetXaxis()->SetTitle("cos#psi");
  frame1.GetYaxis()->SetTitle("P_{3}^{m}(cos#psi)");
  frame1.GetXaxis()->CenterTitle();
  frame1.GetYaxis()->CenterTitle();
  frame1.Draw();

  TF1* p0 = new TF1("p0",".5*(5.*x*x*x - 3.*x)*sqrt(7.)",-1,1);
  p0->SetLineColor(2);
  p0->Draw("same");

  TF1* p2 = new TF1("p2","15*x*(x*x-1.)*sqrt(7./120.)",-1,1);
  p2->SetLineColor(4);
  p2->SetLineStyle(2);
  p2->Draw("same");
  TLine ll;
  ll.SetLineStyle(2);
  ll.DrawLine(-1,0,1,0);

  TLatex tt;
  tt.SetTextColor(2);
  tt.DrawLatex(-.9,-2,"m=0");
  tt.SetTextColor(4);
  tt.DrawLatex(-.6,1.5,"m=2");

  TLatex text;
  text.SetNDC();
  text.SetTextSize(.07);
  text.DrawLatex(.19,.85,"(a)");

  pad2->cd();

  TH2S frame2("frame2","",10,-1,1,10,-4,3);
  frame2.SetStats(kFALSE);
  frame2.GetXaxis()->SetTitle("cos#psi");
  //frame2.GetYaxis()->SetTitle("P_{3}^{m}(cos#psi)");
  frame2.GetYaxis()->SetLabelSize(0.);
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.Draw();

  TF1* p1 = new TF1("p1","-1.5*sqrt(1-x*x)*(5.*x*x - 1.)*sqrt(7./12)",-1,1);
  p1->SetLineColor(2);
  p1->Draw("same");

  TF1* p3 = new TF1("p3","-15.*pow(1-x*x,1.5)*sqrt(7./120)",-1,1);
  p3->SetLineColor(4);
  p3->SetLineStyle(2);
  p3->Draw("same");
  ll.DrawLine(-1,0,1,0);

  tt.SetTextColor(2);
  tt.DrawLatex(-.1,1.5,"m=1");
  tt.SetTextColor(4);
  tt.DrawLatex(-.1,-3.2,"m=3");

  text.DrawLatex(.19,.85,"(b)");
}
