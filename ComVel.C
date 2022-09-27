{
  gROOT->Reset();

  TH2S frame("frame","",10,5,10,10,0,60000);
  frame.GetXaxis()->SetTitle("V_{lab} [cm/ns]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH1I* histX = (TH1I*) file.Get("Be7/vel_7Be_2_el");
  histX->SetMarkerStyle(20);
  histX->Draw("same P");
  float Nx = histX->Integral();

  TFile fileS("sim.root");
  TH1I* sr = (TH1I*) fileS.Get("velR");
  sr->SetLineColor(2);
  float Nr = sr->Integral();

  sr->Scale(Nx/Nr); 
  sr->Draw("L same");


  return;
  TH1I* sp = (TH1I*) fileS.Get("vel");
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
