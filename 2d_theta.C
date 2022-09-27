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
Sty->SetPadTopMargin(0.1);
Sty->SetPadLeftMargin(0.2);
Sty->SetPadRightMargin(0.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kBlack);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(-1.5,"y");
Sty->SetTitleOffset(1.2,"x");
Sty->SetTitleOffset(1.5,"z");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.02,"xyz");
Sty->SetNdivisions(5,"xyz");
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();
  TCanvas canvas("2d_theta_7Be72","",600,1200);

double overlapy = .065;
 TPad *pad1 = new TPad("pad1","",0.,0.,1.,.5+overlapy);
 TPad *pad2 = new TPad("pad2","",0.,.5-overlapy/2.,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

pad1->Draw();
pad2->Draw();


 TLine ll;
 ll.SetLineWidth(3);



  TH2S frame("frame","",10,-1,1,10,0,360);
  frame.GetXaxis()->SetTitle("\\cos(\\psi)");
  frame.GetYaxis()->SetTitle("\\chi \\,\\,[deg]");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.GetYaxis()->SetNdivisions(4,kFALSE);
  frame.GetXaxis()->SetTicks("+-");
  frame.GetYaxis()->SetTicks("+-");
  frame.SetStats(kFALSE);


  TH2S frame2("frame2","",10,-1,1,10,0,360);
  //frame2.GetXaxis()->SetTitle("\\cos(\\psi)");
  frame2.GetXaxis()->SetLabelSize(0.);
frame2.GetYaxis()->SetTitle("\\chi \\,\\, [deg]");
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  //frame2.GetYaxis()->SetLabelSize(0.);
  frame2.GetYaxis()->SetNdivisions(4,kFALSE);
  frame2.SetStats(kFALSE);
  frame2.GetXaxis()->SetTicks("+-");
  frame2.GetYaxis()->SetTicks("+-");


  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH2I* histX_mid = (TH2I*) file.Get("Be7/theta3DeltaPhi_7Be_2_el_midAngle");
  TH2I* histX_hi = (TH2I*) file.Get("Be7/theta3DeltaPhi_7Be_2_el_hiAngle");


  float Nx_mid = histX_mid->Integral();
  float Nx_hi = histX_hi->Integral();


  TFile fileS("sim.root");
  TH2I* sr_mid = (TH2I*) fileS.Get("theta3DeltaPhi_midAngle_R");
  double Nsr_mid = sr_mid->Integral();
  TH2I* sp_mid = (TH2I*) fileS.Get("theta3DeltaPhi_midAngle_P");
  double Nsp_mid = sp_mid->Integral();



  TH2I* sr_hi = (TH2I*) fileS.Get("theta3DeltaPhi_hiAngle_R");
  double Nsr_hi = sr_hi->Integral();
  TH2I* sp_hi = (TH2I*) fileS.Get("theta3DeltaPhi_hiAngle_P");
  double Nsp_hi = sp_hi->Integral();


  TH2F* corr_mid = new TH2F("corr_mid","",50,-1,1,90,0,360);
  TH2F* corr_hi = new TH2F("corr_hi","",50,-1,1,90,0,360);



  for (int i=1;i<=50;i++)
    {
      for (int j=1;j<=90;j++)
	{
	double yr  = sr_mid->GetBinContent(i,j);
	double yp  = sp_mid->GetBinContent(i,j);
	double effic = yr/yp/Nsr_mid*Nsp_mid;
        
        double yexp = histX_mid->GetBinContent(i,j);

	if (yexp == 0 || yr == 0) corr_mid->SetBinContent(i,j,0.);
        else
        corr_mid->SetBinContent(i,j,yexp/effic);


	yr  = sr_hi->GetBinContent(i,j);
	yp  = sp_hi->GetBinContent(i,j);
	effic = yr/yp/Nsr_hi*Nsp_hi;
        
        yexp = histX_hi->GetBinContent(i,j);

	if (yexp == 0 || yr == 0) corr_hi->SetBinContent(i,j,0.);
        else
        corr_hi->SetBinContent(i,j,yexp/effic);
	}
    }

  TLatex ff;
  ff.SetNDC();
  ff.SetTextSize(.06);
  ff.DrawLatex(0.3,0.96,"(a)\\,\\, 2^{\\circ}<\\theta^{*}_{lab} < 4.5^{\\circ}");
  ff.DrawLatex(0.32,0.527,"(b)\\, \\,\\theta^{*}_{lab} > 4.5^{\\circ}");

  TLatex tt;
  //tt.SetNDC();
  pad1->cd();
  frame.Draw();
  corr_mid->Smooth();
  corr_mid->SetMinimum(0.);
  corr_mid->SetMaximum(110.);
  corr_mid->Draw("same zcol");
  //tt.DrawLatex(-.1,300,"1^{#circ}<#theta^{*}_{lab} < 4.5^{#circ}");

  ll.DrawLine(-1.,90,-.92,90);
  ll.DrawLine(0.92,90,1.,90);
  ll.DrawLine(-1.,180,-.92,180);
  ll.DrawLine(0.92,180,1.,180);
  ll.DrawLine(-1.,270,-.92,270);
  ll.DrawLine(0.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0.,0.,0.,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
  ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);
  pad2->cd();
  frame2.Draw();
  corr_hi->Smooth();
  corr_hi->SetMinimum(0.);
  corr_hi->SetMaximum(60.);
  corr_hi->Draw("same zcol");

  //tt.DrawLatex(-.1,300,"#theta^{*}_{lab} > 4.5^{#circ}");

  ll.DrawLine(-1.,90,-.92,90);
  ll.DrawLine(0.92,90,1.,90);
  ll.DrawLine(-1.,180,-.92,180);
  ll.DrawLine(0.92,180,1.,180);
  ll.DrawLine(-1.,270,-.92,270);
  ll.DrawLine(0.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0.,0.,0.,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
  ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);
}
