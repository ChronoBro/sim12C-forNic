{
  gROOT->Reset();

  TCanvas canvas("ComPsiMid","",600,800);
  canvas.Divide(1,2);



  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

  canvas.cd(1);
  TH2S frame2("frame2","",10,-1,1,10,0,4700);
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

   floay y = 0.;
   if (ysim == 0.) continue;
   else  y = yexp*yprim/ysim;


   corrected->SetBinContent(i,y);

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
  SLE.solve();

  cout << "integrals= " << sum33*2./50. << " " << sum32*2./50. << " " << 
    sum31*2./50. << " " <<  sum30*2./50. << endl;


  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3] << endl;


  float sum = 0.;
  for (int i=0;i<4;i++) sum += SLE.Y[i];

  cout << SLE.Y[0]/sum << " " << SLE.Y[1]/sum << " " << SLE.Y[2]/sum << " " << SLE.Y[3]/sum << endl;

  corrected->SetMarkerColor(1);
  corrected->Draw("same P");


  TH1I* fit = (TH1I*) corrected->Clone();
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
  fit->Draw("SAME C");


  fit33->SetLineColor(4);
  fit33->Draw("same C");


  fit32->SetLineColor(3);
  fit32->Draw("same C");


  fit31->SetLineColor(6);
  fit31->Draw("same C");

  fit30->SetLineColor(7);
  fit30->Draw("same C");

  canvas.cd(2);

  TH2S frame3("frame3","",10,-.5,5.5,10,0,2500);
  frame3.SetStats(kFALSE);
  frame3.GetXaxis()->SetTitle("|m_{l}|");
  frame3.GetXaxis()->CenterTitle();
  frame3.Draw();


  double xx[5] = {0.,1.,2.,3.,4.};
  double yy[5];
  for (int i=0;i<4;i++) yy[i] = SLE.Y[i];
  yy[4] = 0.;

  TGraph g(5,xx,yy);
  g.Draw("B");


  TLatex tt;
  tt.DrawLatex(4.,2000,"j = 7/2^{-}, l = 3");
  tt.DrawLatex(.5,2000,"^{7}Be#rightarrow ^{3}He+#alpha,    2^{#circ}<#theta^{*}<4.5^{#circ}");

}
