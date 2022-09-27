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
Sty->SetPadBottomMargin(0.24);
Sty->SetPadTopMargin(0.05);
Sty->SetPadLeftMargin(0.22);
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

 TCanvas canvas("ComPsi7Be72","",400,800);

 double overlap = .00;
TPad * pad1 = new TPad("pad1","",0.,0.,1.,.33 + overlap+.02);
TPad * pad2 = new TPad("pad2","",0.,.33-overlap/2.,1.,.66+overlap/2.+.02);
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
  TH2S frame2("frame2","",10,-1,1,10,0,9000);
  frame2.GetXaxis()->SetTitle("\\cos(\\psi)");
  frame2.GetYaxis()->SetTitle("Counts");
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.SetStats(kFALSE);
  frame2.Draw();

  TLatex tt;
  tt.SetTextSize(.08);
  tt.SetNDC();
  tt.DrawLatex(.25,.85,"(a)");


  

 

  TH1I* histX3 = (TH1I*) file.Get("Be7/theta3_7Be_2_el");
  histX3->SetMarkerStyle(20);
  //histX3->Draw("same P");
  float Nx3 = histX3->Integral();


  TH1I* sr3 = (TH1I*) fileS.Get("theta3_R");
  sr3->SetLineColor(2);
  float Nr3 = sr3->Integral();

  sr3->Scale(Nx3/Nr3); 
  //sr3->Draw("L same");


  TH1I* sp3 = (TH1I*) fileS.Get("theta3_P");
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



  double sig[4];

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
      sig[i] = error;
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
  fit33->SetLineStyle(2);
  fit33->Draw("same C");


  fit32->SetLineColor(3);
  fit32->SetLineStyle(2);
  fit32->Draw("same C");


  fit31->SetLineColor(6);
  fit31->SetLineStyle(2);
  fit31->Draw("same C");

  fit30->SetLineColor(7);
  fit30->SetLineStyle(2);
  fit30->Draw("same C");


  TLatex jt;
  jt.SetTextColor(7);
  jt.SetTextSize(.075);
  jt.DrawLatex(.4,400,"0");
  jt.SetTextColor(6);
  jt.DrawLatex(-.02,1000,"1");
  jt.SetTextColor(3);
  jt.DrawLatex(-.4,2100,"2");
  jt.SetTextColor(4);
  jt.DrawLatex(-.30,4000,"3");


  pad2->cd();

  TH2S frame3("frame3","",10,-4.5,4.5,10,0,0.4);
  frame3.SetStats(kFALSE);
  frame3.GetXaxis()->SetTitle("\\ m");
  frame3.GetYaxis()->SetTitle("\\rho^{\\ell}_{m,m}");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.GetXaxis()->SetNdivisions(9);
  frame3.GetXaxis()->SetTicks("+-");
  frame3.GetXaxis()->SetTitleOffset(1.1);
  frame3.Draw();
  tt.DrawLatex(.25,.85,"(b)");



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


  float AAA = 0.;
  float sigAAA = 0.;
  for (int i=0;i<7;i++) 
    {
     AAA += yy[i]*(3.*pow(xx[i],2)-12.)/15.;
     //sigAAA += pow(syyy[i]*(3.*pow(xxx[i],2)-63./4.)/(42./2.),2);
     cout << xx[i] << " " << (3.*pow(xx[i],2)-12.)/15. << endl;
    }

  cout << "alignment = " << AAA << endl; 


  /*
  TLatex tt;
  tt.DrawLatex(4.,2000,"j = 7/2^{-}, l = 3");
  tt.DrawLatex(.5,2000,"^{7}Be#rightarrow ^{3}He+#alpha,    2^{#circ}<#theta^{*}<4.5^{#circ}");
  */

  pad1->cd();
  TH2S frame4("frame4","",10,-4.5,4.5,10,0,0.4);
  frame4.SetStats(kFALSE);
  frame4.GetXaxis()->SetTitleOffset(1.3);
  frame4.GetXaxis()->SetTitle("\\ m^*");
  frame4.GetYaxis()->SetTitle("\\rho^{s^*}_{m^*,m^*}");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.GetXaxis()->SetNdivisions(9,kFALSE);
  frame4.GetXaxis()->SetTicks("+-");
  frame4.GetXaxis()->SetLabelSize(0.);
  frame4.Draw();

  TLatex ft;
  ft.SetTextSize(.06);
  ft.DrawLatex(-3.7,-.07,"- #frac{7}{2}");
  ft.DrawLatex(-2.7,-.07,"- #frac{5}{2}");
  ft.DrawLatex(-1.7,-.07,"- #frac{3}{2}");
  ft.DrawLatex(-0.7,-.07,"- #frac{1}{2}");
  ft.DrawLatex(0.4,-.07,"#frac{1}{2}");
  ft.DrawLatex(1.4,-.07,"#frac{3}{2}");
  ft.DrawLatex(2.4,-.07,"#frac{5}{2}");
  ft.DrawLatex(3.4,-.07,"#frac{7}{2}");
  tt.DrawLatex(.25,.85,"(c)");

  double xxx[8]={-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5};
  double yyy[8];
  double syyy[8];

  yyy[3] = 7./4.*SLE.Y[0]/sum;
  yyy[4] = yyy[3];
  syyy[3] = 7./4.*sig[0];
  syyy[4] = syyy[3];

  yyy[2] = 7./5.*(SLE.Y[1]/sum - 3./7.*yyy[3]);
  yyy[5] = yyy[2];
  syyy[2] = sqrt(pow(7./5.*sig[1],2)+pow(3./7.*syyy[3],2));
  syyy[5] = syyy[2];  

  yyy[1] = 7./6.*(SLE.Y[2]/sum - 2./7*yyy[2]);
  yyy[6] = yyy[1];
  syyy[1] = sqrt(pow(7./6.*sig[2],2)+pow(2./7.*syyy[2],2));
  syyy[6] = syyy[1];

  yyy[0] = SLE.Y[3]/sum - yyy[1]/7.;
  yyy[7] = yyy[0];
  syyy[0] = sqrt(pow(sig[3],2)+pow(syyy[1],2));

  for (int i=0;i<8;i++) yyy[i] /=2.;
  for (int i=0;i<8;i++) syyy[i] /=2.;

  TGraph gg(8,xxx,yyy);
  gg.Draw("B");

  cout << " 1/2 " << yyy[3] << " " <<syyy[3]<< endl;
  cout << " 3/2 " << yyy[2] << " " <<syyy[2]<<endl;
  cout << " 5/2  " << yyy[1] << " " << syyy[1] << endl;
  cout << " 7/2 " << yyy[0] << " " << syyy[0] << endl;

  double gx[8] = {-7./2.,-5./2.,-3./2.,-1./2.,1./2.,3./2.,5./2.,7./2.};
  double gy[8] = {2./3.,8./21.,2./21.,0.,0.,2./21,8./21.,2./3.};


  for (int i=0;i<8;i++) gy[i] /= 2.5;
  TGraph ggg(8,gx,gy);
  ggg.SetMarkerStyle(21);
  ggg.SetMarkerColor(2);
  ggg.SetLineColor(2);
  //ggg.Draw("PL");


  //alignment

  float AA = 0.;
  float sigAA = 0.;
  for (int i=0;i<8;i++) 
    {
     AA += yyy[i]*(3.*pow(xxx[i],2)-63./4.)/(42./2.);
     sigAA += pow(syyy[i]*(3.*pow(xxx[i],2)-63./4.)/(42./2.),2);
     cout << xxx[i] << " " << (3.*pow(xxx[i],2)-63./4.)/(42./2.) << endl;
    }

  cout << "alignment = " << AA << " " << sqrt(sigAA) << endl; 

}
