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

 TCanvas canvas("Psi_theta","",1200,600);



double overlapy = .06;
TPad *pad1 = new TPad("pad1","",0.,0.,.5+overlapy,1.);
TPad *pad2 = new TPad("pad2","",.5-overlapy/2.,0.,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

pad1->Draw();
pad2->Draw();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

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


  SLE.solve();

  float sum = SLE.Y[0] + SLE.Y[1] + SLE.Y[2] + SLE.Y[3];

  pad1->cd();

  TH2S frame1("frame1","",10,-4.5,4.5,10,0,0.4);
  frame1.SetStats(kFALSE);
  frame1.GetXaxis()->SetTitle("m");
  frame1.GetYaxis()->SetTitle("p_{m}");
  frame1.GetXaxis()->CenterTitle();
  frame1.GetYaxis()->CenterTitle();
  frame1.Draw();


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
  tt.DrawLatex(-1.3,.35,"2^{#circ}<#theta_{lab}^{*}<4.5^{#circ}");

  //***************************************************
  TH1I* histX4 = (TH1I*) file.Get("Be7/theta3_7Be_2_el_hiAngle");
  histX3->SetMarkerStyle(20);
  //histX3->Draw("same P");
  float Nx4 = histX4->Integral();



  TH1I* sr4 = (TH1I*) fileS.Get("theta3_hiAngle_R");
  sr4->SetLineColor(2);
  float Nr4 = sr4->Integral();

  sr4->Scale(Nx4/Nr4); 
  //sr3->Draw("L same");


  TH1I* sp4 = (TH1I*) fileS.Get("theta3_hiAngle_P");
  sp4->SetLineColor(4);
  float Np4 = sp4->Integral();

  sp4->Scale(Nx4/Np4); 
  sp4->SetLineWidth(3);

  SLE.clear();



  for (int i=1;i<=50;i++)
   {
   float yexp = histX4->GetBinContent(i);

   float ysim = sr4->GetBinContent(i);
   float yprim = sp4->GetBinContent(i);

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




   SLE.Y[0] += y30*y;
   SLE.Y[1] += y31*y;
   SLE.Y[2] += y32*y;
   SLE.Y[3] += y33*y;
  

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


  SLE.solve();

  sum = SLE.Y[0] + SLE.Y[1] + SLE.Y[2] + SLE.Y[3];

  pad2->cd();

  TH2S frame2("frame2","",10,-4.5,4.5,10,0,0.4);
  frame2.SetStats(kFALSE);
  frame2.GetXaxis()->SetTitle("m");
  frame2.GetYaxis()->SetLabelSize(0.);
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.Draw();



  yy[0] = SLE.Y[3]/2./sum;
  yy[1] = SLE.Y[2]/2./sum;
  yy[2] = SLE.Y[1]/2./sum;
  yy[3] = SLE.Y[0]/sum;  
  yy[4] = yy[2];
  yy[5] = yy[1];
  yy[6] = yy[0];


  TGraph g4(7,xx,yy);
  g4.Draw("B");


  tt.DrawLatex(-1,.35,"#theta_{lab}^{*}>4.5^{#circ}");


}
