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

 TCanvas canvas("ComPsi","",400,800);

 double overlap = .00;
TPad * pad1 = new TPad("pad1","",0.,0.,1.,.33 + overlap);
TPad * pad2 = new TPad("pad2","",0.,.33-overlap/2.,1.,.66+overlap/2.);
TPad * pad3 = new TPad("pad3","",0.,.66-overlap,1.,1.);

pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);
pad3->SetFillStyle(4000);

pad1->Draw();
pad2->Draw();
pad3->Draw();


  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

  pad3->cd();
  TH2S frame2("frame2","",10,-1,1,10,0,5000);
  frame2.GetXaxis()->SetTitle("cos(#psi)");
  frame2.GetYaxis()->SetTitle("Counts");
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.SetStats(kFALSE);
  frame2.Draw();


  TH1I* histX3 = (TH1I*) file.Get("Be7/theta3_7Be_2_el_midAngle");
  histX3->SetMarkerStyle(20);
  //histX3->Draw("same P");
  float Nx3 = histX3->Integral();



  TH1I* sr3 = (TH1I*) fileS.Get("theta3_midAngle_R");
  sr3->SetLineColor(2);
  float Nr3 = sr3->Integral();

  sr3->Scale(Nx3/Nr3); 
  //sr3->Draw("L same");


  TH1I* sp3 = (TH1I*) fileS.Get("theta3_midAngle_P");
  sp3->SetLineColor(4);
  float Np3 = sp3->Integral();

  sp3->Scale(Nx3/Np3); 
  sp3->SetLineWidth(3);
  //sp3->Draw("L same");



  sle SLE(4);
  SLE.clear();
  inverse Inverse(4);

  float sum33 = 0.;
  float sum32 = 0.;
  float sum31 = 0.;
  float sum30 = 0.;

  TH1I* corrected = (TH1I*) histX3->Clone();
  for (int i=1;i<=50;i++)
   {
   float yexp = histX3->GetBinContent(i);
   //cout << yexp << endl;
   float ysim = sr3->GetBinContent(i);
   float yprim = sp3->GetBinContent(i);

   float y = 0.;
   if (ysim == 0.) continue;
   else  y = yexp*yprim/ysim;


   corrected->SetBinContent(i,y);
   corrected->SetBinError(i,y/sqrt(yexp));
   float cosx = histX3->GetBinCenter(i);


   float y33 = 225.*pow(1.-pow(cosx,2),3)*7./720.;
   float y32 = 225.*pow(cosx,2)*pow(pow(cosx,2)-1.,2)*7./120.;
   float y31 = (1.-pow(cosx,2))*pow(5.*pow(cosx,2)-1.,2)*9./4.*7./12.;
   float y30 = pow(5.*pow(cosx,3)-3.*cosx,2)/4.*7.;


   sum33 += y33;
   sum32 += y32;
   sum31 += y31;
   sum30 += y30;
   


   //cout << y30 << " " << y31 << " " << y32 << " " << y33 << endl;


   SLE.Y[0] += y30*y;
   SLE.Y[1] += y31*y;
   SLE.Y[2] += y32*y;
   SLE.Y[3] += y33*y;
  
   //cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3]; endl;

   SLE.M[0][0] += y30*y30;
   SLE.M[1][0] += y30*y31;
   SLE.M[2][0] += y30*y32;
   SLE.M[3][0] += y30*y33;




   SLE.M[1][1] += y31*y31;
   SLE.M[2][1] += y31*y32;
   SLE.M[3][1] += y31*y33;

   SLE.M[2][2] += y32*y32;
   SLE.M[3][2] += y32*y33;

   SLE.M[3][3] += y33*y33;
   }



  SLE.M[0][1] = SLE.M[1][0];
  SLE.M[0][2] = SLE.M[2][0];
  SLE.M[0][3] = SLE.M[3][0];
  SLE.M[1][2] = SLE.M[2][1];
  SLE.M[1][3] = SLE.M[3][1];
  SLE.M[2][3] = SLE.M[3][2];

  for (int i=0;i<4;i++)
    for (int j=0;j<4;j++)
      Inverse.A[i][j] = SLE.M[i][j];
  SLE.solve();
  Inverse.solve();
  Inverse.S3000();

  cout << "integrals= " << sum33*2./50. << " " << sum32*2./50. << " " << 
    sum31*2./50. << " " <<  sum30*2./50. << endl;


  double tot = SLE.Y[0] + SLE.Y[1] + SLE.Y[2] + SLE.Y[3];
  for (int i=0;i<4;i++)
    {
      double mag = SLE.Y[i]/tot;
      double error = sqrt(Inverse.AI[i][i])/tot;
      if (i>0)
	{
	  mag /= 2.;
	  error /= 2.;
	}
      cout << i << " " << mag << " " << error << endl;
    }

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3] << endl;


  float sum = 0.;
  for (int i=0;i<4;i++) sum += SLE.Y[i];

  cout << SLE.Y[0]/sum << " " << SLE.Y[1]/sum << " " << SLE.Y[2]/sum << " " << SLE.Y[3]/sum << endl;

  corrected->SetMarkerColor(1);
  corrected->Draw("same P");


  TH1I* fit = (TH1I*) histX3->Clone();
  fit->SetMarkerStyle(0);
  fit->SetLineStyle(1);
  TH1I* fit33 =  (TH1I*) fit->Clone();
  TH1I* fit32 =  (TH1I*) fit->Clone();
  TH1I* fit31 =  (TH1I*) fit->Clone();
  TH1I* fit30 =  (TH1I*) fit->Clone();
  for (int i=1;i<=50;i++)
   {
   float cosx = histX3->GetBinCenter(i);


   float y33 = 225.*pow(1.-pow(cosx,2),3)*7./720.;
   float y32 = 225.*pow(cosx,2)*pow(pow(cosx,2)-1.,2)*7./120.;
   float y31 = (1.-pow(cosx,2))*pow(5.*pow(cosx,2)-1.,2)*9./4.*7./12.;
   float y30 = pow(5.*pow(cosx,3)-3.*cosx,2)/4.*7.;


   y33 *= SLE.Y[3];
   y32 *= SLE.Y[2];
   y31 *= SLE.Y[1];
   y30 *= SLE.Y[0];
   
   fit->SetBinContent(i,y33+y32+y31+y30);
   fit33->SetBinContent(i,y33);
   fit32->SetBinContent(i,y32);
   fit31->SetBinContent(i,y31);
   fit30->SetBinContent(i,y30);

   }

  fit->SetLineColor(2);
  fit->SetLineWidth(3);
  fit->Draw("SAME L");



  fit33->SetLineColor(4);
  fit33->Draw("same C");


  fit32->SetLineColor(3);
  fit32->Draw("same C");


  fit31->SetLineColor(6);
  fit31->Draw("same C");

  fit30->SetLineColor(7);
  fit30->Draw("same C");

  pad2.cd();

  TH2S frame3("frame3","",10,-4.5,4.5,10,0,0.4);
  frame3.SetStats(kFALSE);
  frame3.GetXaxis()->SetTitle("m");
  frame3.GetYaxis()->SetTitle("p_{m}");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.Draw();


  double xx[7] = {-3,-2,-1,0.,1.,2.,3.};
  double yy[7];
  yy[0] = SLE.Y[3]/2./sum;
  yy[1] = SLE.Y[2]/2./sum;
  yy[2] = SLE.Y[1]/2./sum;
  yy[3] = SLE.Y[0]/sum;  
  yy[4] = yy[2];
  yy[5] = yy[1];
  yy[6] = yy[0];


  TGraph g(7,xx,yy);
  g.Draw("B");


  TLatex tt;
  tt.DrawLatex(4.,2000,"j = 7/2^{-}, l = 3");
  tt.DrawLatex(.5,2000,"^{7}Be#rightarrow ^{3}He+#alpha,    2^{#circ}<#theta^{*}<4.5^{#circ}");


  pad1.cd();
  TH2S frame4("frame4","",10,-4.5,4.5,10,0,0.4);
  frame4.SetStats(kFALSE);
  frame4.GetXaxis()->SetTitle("m*");
  frame4.GetYaxis()->SetTitle("p_{m*}");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.Draw();

  double xxx[8]={-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5};
  double yyy[8];

  yyy[3] = 7./4.*SLE.Y[0]/sum;
  yyy[4] = yyy[3];

  yyy[2] = 7./5.*(SLE.Y[1]/sum - 3./7.*yyy[3]);
  yyy[5] = yyy[2];
  
  yyy[1] = 7./6.*(SLE.Y[2]/sum - 2./7*yyy[2]);
  yyy[6] = yyy[1];

  yyy[0] = SLE.Y[3]/sum - yyy[1]/7.;
  yyy[7] = yyy[0];

  for (int i=0;i<8;i++) yyy[i] /=2.;

  TGraph gg(8,xxx,yyy);
  gg.Draw("B");


}
