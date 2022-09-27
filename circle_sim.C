{
  gROOT->Reset();
  TCanvas canvas("circle","",800,800);
  canvas.Divide(2,2);  
canvas.cd(1);

 TH2S frame1("frame1","",10,-2,12,10,-6,6);
 frame1.GetXaxis()->SetTitle("v_{Z} [cm/ns]");
 frame1.GetYaxis()->SetTitle("v_{X} [cm/ns]");
 frame1.SetStats(kFALSE);
 frame1.Draw();

  TFile f("sim.root");

  TH2I* vzvx_H2 = (TH2I*) f.Get("vzvx_H2");
  vzvx_H2->Draw("cont same");

  //TH1I* vzvx_T = f.Get("Li6/vzvx_6Li_da_el_T");
  //vzvx_T->Draw("cont same");
  double v = 3.33767;
  double v1 = 3.31892;
  double v2 = 4.94074;
  double x[361];
  double y[361];
  double x2[361];
  double y2[361];
  for (int i=0;i<361;i++)
    {
      double theta = (double)i*acos(-1.)/180.;

      x[i] = (v + v1*cos(theta))/(1.+ v*v1*cos(theta)/pow(30.,2));
      y[i] = v1*sin(theta)*sqrt(1.-pow(v/30.,2))/(1.+ v*v1*cos(theta)/pow(30.,2));


      x2[i] = (v + v2*cos(theta))/(1.+ v*v2*cos(theta)/pow(30.,2));
      y2[i] = v2*sin(theta)*sqrt(1.-pow(v/30.,2))/(1.+ v*v2*cos(theta)/pow(30.,2));


    }
  TGraph gcircle(361,x,y);
  gcircle.Draw("C");
  TGraph gcircle2(361,x2,y2);
  gcircle2.SetLineStyle(2);
  gcircle2.Draw("C");



  TLine l;
  l.SetLineStyle(2);
  l.DrawLine(-2,0,14,0);
  l.DrawLine(0.,-6,0.,6); 

  TArrow arrow;
  arrow.SetAngle(30);
  arrow.SetFillColor(1);
  arrow.SetLineWidth(3);
  arrow.DrawArrow(0.,0.,4.818,0.,.03,"|>");

  //TEllipse ellipse(4.818,0.,4.818,4.818);
  //ellipse.SetFillStyle(4100);
  //ellipse.Draw("same");

  TLatex text;
  text.SetTextSize(.06);
  text.DrawLatex(9,3.5,"d");
  text.DrawLatex(-1.5,-1,"^{9}Be");


  canvas.cd(2);
  TH2S frame2("frame2","",10,-6,6,10,-6,6);
  frame2.GetXaxis()->SetTitle("v_{Y} [cm/ns]");
  frame2.GetYaxis()->SetTitle("v_{X} [cm/ns]");
  frame2.SetStats(kFALSE);
  frame2.Draw();

  TH2I *  swap_H2 = new TH2I("swap_H2","",100,-7,7,100,-7,7);
  TH2I* vxvy_H2 = (TH2I*) f.Get("vxvy_H2");
  for (int i=1;i<=100;i++)
    for (int j=1;j<=100;j++)
      {
        double yy = vxvy_H2->GetBinContent(i,j);
        swap_H2->SetBinContent(j,i,yy);
      }
  swap_H2->Draw("col same");

  l.DrawLine(-6,.0.,6.,0.);
  l.DrawLine(0.,-6,0.,6);
  text.DrawLatex(3,3,"d");
  canvas.cd(3);
  frame1.Draw();
  TH2I* vzvx_He4 = (TH2I*) f.Get("vzvx_He4");
  vzvx_He4->Draw("cont same");
  //ellipse.Draw("same");
  gcircle.Draw("same");
  gcircle2.Draw("same");


  l.DrawLine(-2,0,14,0);
  l.DrawLine(0.,-7,0.,7); 
  arrow.DrawArrow(0.,0.,4.818,0.,.03,"|>");

  vzvx_T->Draw("cont same");


  text.DrawLatex(9,2.,"#alpha");
  text.DrawLatex(-1.5,-1,"^{9}Be");


  canvas.cd(4);
  frame2.Draw();


  TH2I *  swap_He4 = new TH2I("swap_He4","",100,-7,7,100,-7,7);

  TH2I* vxvy_He4 = (TH2I*) f.Get("vxvy_He4");

  for (int i=1;i<=100;i++)
    for (int j=1;j<=100;j++)
      {
        double yy = vxvy_He4->GetBinContent(i,j);
        swap_He4->SetBinContent(j,i,yy);
      }
  swap_He4->Draw("col same");

  l.DrawLine(-6,.0.,6.,0.);
  l.DrawLine(0.,-3,0.,6);
  text.DrawLatex(2,2.,"#alpha");
  text.DrawLatex(-4,-4,"E/A=36.6 MeV ^{6}Li+^{9}Be");
  text.DrawLatex(-4,-5,"^{6}Li#rightarrow ^{2}H+#alpha");
}
