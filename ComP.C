{
  gROOT->Reset();
  TCanvas canvas("P","",600,800);
  canvas.Divide(1,3);
  canvas.cd(1);
  TH2S frame("frame","",10,1600,2800,10,0,40000);
  frame.GetXaxis()->SetTitle("P_{z} [MeV/c]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH1I* histZ = (TH1I*) file.Get("Be7/Ppara_7Be_2_el");
  histZ->SetMarkerStyle(20);
  histZ->Draw("same P");
  float Nz = histZ->Integral();

  TFile fileS("sim.root");
  TH1I* sr = (TH1I*) fileS.Get("Ppara");
  sr->SetLineColor(2);
  float Nr = sr->Integral();

  sr->Scale(Nz/Nr); 
  sr->Draw("L same");

  canvas.cd(2);

  TH2S frame2("frame2","",10,-500,500,10,0,4000);
  frame2.GetXaxis()->SetTitle("P_{x} [MeV/c]");
  frame2.SetStats(kFALSE);
  frame2.Draw();


  TH1I* histX = (TH1I*) file.Get("Be7/Px_7Be_2_el");
  histX->SetMarkerStyle(20);
  histX->Draw("same P");
  float Nx = histX->Integral();


  TH1I* srx = (TH1I*) fileS.Get("Px");
  srx->SetLineColor(2);
  float Nrx = srx->Integral();

  srx->Scale(Nx/Nrx); 
  srx->Draw("L same");



  canvas.cd(3);

  TH2S frame3("frame3","",10,-500,500,10,0,4000);
  frame3.GetXaxis()->SetTitle("P_{y} [MeV/c]");
  frame3.SetStats(kFALSE);
  frame3.Draw();


  TH1I* histY = (TH1I*) file.Get("Be7/Py_7Be_2_el");
  histY->SetMarkerStyle(20);
  histY->Draw("same  P");
  float Ny = histY->Integral();


  TH1I* sry = (TH1I*) fileS.Get("Py");
  sry->SetLineColor(2);
  float Nry = sry->Integral();

  sry->Scale(Ny/Nry); 
  sry->Draw("L same");

}
