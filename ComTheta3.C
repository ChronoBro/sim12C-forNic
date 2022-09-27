{
  gROOT->Reset();

  TCanvas canvas("ComTheta3","",600,800);
  canvas.Divide(1,3);

  TLatex tt;
  tt.SetTextSize(.07);
  tt.SetNDC();

  canvas.cd(1);



  TH2S frame1("frame1","",10,-1,1,10,0,1000);
  frame1.GetXaxis()->SetTitle("cos(#psi_{^{2}H})");
  frame1.SetStats(kFALSE);
  frame1.Draw();
  tt.DrawLatex(.2,.8,"cos(#Delta#phi) < -0.5");
  TFile file("/home/Carbon8/daq/corr_6Li.root");
  TH1I* histX = (TH1I*) file.Get("Li6/theta3_6Li_1_el_center");
  histX->SetMarkerStyle(20);
  histX->Draw("same P");
  float Nx = histX->Integral();

  TFile fileS("sim.root");
  TH1I* sr = (TH1I*) fileS.Get("theta3_R_center");
  sr->SetLineColor(2);
  float Nr = sr->Integral();

  sr->Scale(Nx/Nr); 
  sr->Draw("L same");


  float chi2_center = 0.;

  for (int i=1;i<=50;i++)
    {
      chi2_center += pow(histX->GetBinContent(i)-sr->GetBinContent(i),2);
    }

  cout << chi2_center << endl;
  canvas.cd(2);
  TH2S frame2("frame2","",10,-1,1,10,0,1500);
  frame2.GetXaxis()->SetTitle("cos(#psi_{^{2}H})");
  frame2.SetStats(kFALSE);
  frame2.Draw();
  tt.DrawLatex(.2,.8,"-0.5 < cos(#Delta#phi) < 0.5");

  TH1I* histX_side = (TH1I*) file.Get("Li6/theta3_6Li_1_el_side");
  histX_side->SetMarkerStyle(20);
  histX_side->Draw("same P");
  float Nx_side = histX_side->Integral();


  TH1I* sr_side = (TH1I*) fileS.Get("theta3_R_side");
  sr_side->SetLineColor(2);
  float Nr_side = sr_side->Integral();

  sr_side->Scale(Nx_side/Nr_side); 
  sr_side->Draw("L same");

  float chi2_side = 0.;
  for (int i=1;i<=50;i++)
    {
      chi2_side += pow(histX_side->GetBinContent(i)-sr_side->GetBinContent(i),2);
    }
  cout << chi2_side << endl;

  canvas.cd(3);
  TH2S frame3("frame3","",10,-1,1,10,0,1000);
  frame3.GetXaxis()->SetTitle("cos(#psi_{^{2}H})");
  frame3.SetStats(kFALSE);
  frame3.Draw();
  tt.DrawLatex(.2,.8,"cos(#Delta#phi) > 0.5");

  TH1I* histX_away = (TH1I*) file.Get("Li6/theta3_6Li_1_el_away");
  histX_away->SetMarkerStyle(20);
  histX_away->Draw("same P");
  float Nx_away = histX_away->Integral();


  TH1I* sr_away = (TH1I*) fileS.Get("theta3_R_away");
  sr_away->SetLineColor(2);
  float Nr_away = sr_away->Integral();

  sr_away->Scale(Nx_away/Nr_away); 
  sr_away->Draw("L same");

  float chi2_away = 0.;
  for (int i=1;i<=50;i++)
    {
      chi2_away += pow(histX_away->GetBinContent(i)-sr_away->GetBinContent(i),2);
    }

  cout << chi2_center << " " << chi2_side << " " << chi2_away << " " << chi2_center+chi2_side+chi2_away << endl;

  return;

  TH1I* sp = (TH1I*) fileS.Get("theta");
  sp->SetLineColor(4);
  float Np = sp->Integral();

  sp->Scale(Nx/Np); 
  sp->Draw("L same");

  return;

  double xx[100];
  double yy[100];
  int N = 0;;
  for (int i=1;i<200;i++)
    {
      float x = sr->GetBinCenter(i);


      if (x < 7.7) continue;
      if (x > 11.1) break;


      float y = sr->GetBinContent(i);

      float y1 = histX->GetBinContent(i);
      cout << x << " " << y / y1 << endl;

      xx[N] = x;
      yy[N] = y/y1;
      N++;
    }

  TGraph g(N,xx,yy);
  g.Draw("A*");
  g.Fit("pol1");


}
