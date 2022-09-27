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
Sty->SetPadBottomMargin(0.17);
Sty->SetPadTopMargin(0.07);
Sty->SetPadLeftMargin(0.16);
Sty->SetPadRightMargin(0.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kBlack);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.2,"y");
Sty->SetTitleOffset(1.2,"x");
Sty->SetTitleOffset(1.5,"z");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.04,"xyz");
Sty->SetNdivisions(5,"xyz");
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();

 TCanvas canvas("ComThetaCM72");
  gPad->SetLogy();

  //TH2S frame("frame","",10,0,20,10,0,350000);
  TH2S frame("frame","",10,0,16,10,3000,350000);
  frame.GetXaxis()->SetTitle("#theta^{*}_{cm} [deg]");
  frame.GetYaxis()->SetTitle("d#sigma/d#Omega [rel]");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();


  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH1I* histX = (TH1I*) file.Get("Be7/thetaCM_7Be_2_el");
  histX->SetMarkerStyle(20);
  //histX->Draw("same P");
  float Nx = histX->Integral();


  TH1F* domega = new TH1F("domega","",125,0,40);
  TH1F* domega_R = new TH1F("domega_R","",125,0,40);
  TH1F* domega_P = new TH1F("domega_P","",125,0,40);



  TFile fs("sim.root");
  TH1I* theta_P = (TH1I*) fs.Get("theta_reactionCM_P");
  TH1I* theta_R = (TH1I*) fs.Get("theta_reactionCM_R");


  float N_P = theta_P.Integral();
  float N_R = theta_R.Integral();


  for (int i=1;i<=100;i++)
    {
      double y = histX->GetBinContent(i);
      double theta = histX->GetBinCenter(i);
      theta *= acos(-1.)/180.;

      double yomega = y/sin(theta);
      domega->SetBinContent(i,yomega);
      domega->SetBinError(i,sqrt(y)/sin(theta));

      //domega->SetBinContent(i,y);
      //domega->SetBinError(i,sqrt(y));
     
      y = theta_P.GetBinContent(i);
      domega_P ->SetBinContent(i,y*Nx/N_P/sin(theta));


      y = theta_R.GetBinContent(i);
      domega_R ->SetBinContent(i,y*Nx/N_R/sin(theta));

    }
  domega->SetMarkerStyle(20);
  domega->SetMarkerColor(4);
  domega->SetLineColor(4);
  domega->Draw("same P");
  double Nexp = domega->Integral();

  domega_R->SetLineColor(6);
  domega_R->Draw("same L");

  domega_P->SetLineColor(1);
  domega_P->SetLineWidth(3);
  //domega_P->Draw("same L");


  TH1F* domega_P_up = (TH1F*)domega_P->Clone();
  TH1F* domega_P_down = (TH1F*) domega_P->Clone();



  for (int i=1;i<100;i++)
    {
      float yy= domega_P->GetBinContent(i);
      float x = domega_P->GetBinCenter(i);
      float yy1 = yy+yy/3.*exp(-x/1.5);
      domega_P_up->SetBinContent(i,yy1);
      float yy2 = yy-yy/3.*exp(-x/1.5);
      domega_P_down->SetBinContent(i,yy2);

    }

  domega_P_up->Draw("same L");
  domega_P_down->Draw("same L");

  double  fx[300];
  double fy[300];
  int nnn = 0;
  for (int i=1;i<100;i++)
    {
      fx[nnn] = domega_P_down->GetBinCenter(i);
      fy[nnn] = domega_P_down->GetBinContent(i);
      nnn++;
    }
  for (int i=99;i>=1;i--)
    {

      fx[nnn] = domega_P_up->GetBinCenter(i);
      fy[nnn] = domega_P_up->GetBinContent(i);
      nnn++;
    }
  fx[nnn] = fx[0];
  fy[nnn] = fy[0];

  TGraph fg(nnn+1,fx,fy);
  fg.SetFillColor(1);
  fg.SetFillStyle(3013);
  fg.Draw("f");

  TLatex tt;
  tt.SetTextSize(.06);
  tt.DrawLatex(7.,200000.,"J=7/2^{-}, ^{} ^{7}Be");

  return;

  adist Adist;

  TH1F* domega0 = new TH1F("domega0","",125,0,40);
  TH1F* domega1 = new TH1F("domega1","",125,0,40);
  TH1F* domega2 = new TH1F("domega2","",125,0,40);


  for (int i=1;i<=100;i++)
    {
      double theta = histX->GetBinCenter(i);
      theta *= acos(-1.)/180.;
      double y0 = Adist.getDist(theta,47.,5.);
      double y1 = Adist.out1;
      double y2 = Adist.out2;

      //cout << y0 << " " << y1 << endl;

      domega0->SetBinContent(i,pow(y0,2));
      domega1->SetBinContent(i,pow(y1,2));
      domega2->SetBinContent(i,2.*pow(y2,2) + 3.*pow(y0,2));     
    }


  double N0 = domega0->Integral();
  double N1 = domega1->Integral();
  double N2 = domega2->Integral();


  domega0->Scale(Nexp/N0/2.6);
  domega1->Scale(Nexp/N1/1.3);
  domega2->Scale(Nexp/N2/1.2);

  domega0->SetLineColor(2);
  domega1->SetLineColor(4);
  domega2->SetLineColor(3);

  domega0->Draw("C same");
  domega1->Draw("C same");
  domega2->Draw("C same");

  TLatex tt;
  tt.DrawLatex(8.,250000,"L=47 hbar");
}
