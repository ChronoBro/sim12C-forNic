{
  gROOT->Reset();

  TH2S frame("frame","",10,0,360,10,0,5000);
  frame.GetXaxis()->SetTitle("#DeltaPhi [deg]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TFile file("/home/Carbon8/daq/corr_6Li.root");
  TH1I* histX = (TH1I*) file.Get("Li6/deltaPhi_6Li_1_el");
  histX->SetMarkerStyle(20);
  histX->Draw("same P");
  float Nx = histX->Integral();

  TFile fileS("sim.root");
  TH1I* sr = (TH1I*) fileS.Get("deltaPhi_R");
  sr->SetLineColor(2);
  float Nr = sr->Integral();

  sr->Scale(Nx/Nr); 
  sr->Draw("L same");

  //return;
  /*
  TH1I* sp = (TH1I*) fileS.Get("vel");
  sp->SetLineColor(4);
  float Np = sp->Integral();

  sp->Scale(Nx/Np); 
  sp->Draw("L same");
  */
  //return;

  double xx[90];
  double yy[90];
  int N = 0;;
  for (int i=1;i<90;i++)
    {
      float x = sr->GetBinCenter(i);


      float y = sr->GetBinContent(i);

      float y1 = histX->GetBinContent(i);
      cout << x << " " << y / y1 << endl;

      xx[N] = x;
      yy[N] = y1/y*3000;
      N++;
    }

  TGraph g(N,xx,yy);
  g.SetLineColor(6);
  g.Draw("L");
  //g.Fit("pol1");


}
