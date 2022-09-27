{
  gROOT->Reset();
  gROOT->SetStyle("Pub");
  TCanvas canvas("ComAngle");

  canvas.Divide(2,2);


  TFile file("/home/Carbon8/daq/corr_6Li.root");
   TFile fileS("sim.root");


  canvas.cd(1);
  TH2S frame1("frame1","",10,-1,1,10,0,2000);
  frame1.GetXaxis()->SetTitle("cos(#theta_{d-#alpha})");
  frame1.GetYaxis()->SetTitle("Counts");
  frame1.GetXaxis()->CenterTitle();
  frame1.GetYaxis()->CenterTitle();
  frame1.SetStats(kFALSE);
  frame1.Draw();

  TLatex text;
  text.SetTextSize(.08);
  text.DrawLatex(-.2,600,"2^{#circ}>#theta_{plf}>6^{#circ}");


  TH1I* histX1 = (TH1I*) file.Get("Li6/theta3_6Li_1_peak1_el");
  histX1->SetMarkerStyle(20);
  histX1->Draw("same P");
  float Nx1 = histX1->Integral();


  TH1I* sr1 = (TH1I*) fileS.Get("theta3_R_peak1");
  sr1->SetLineColor(2);
  float Nr1 = sr1->Integral();

  sr1->Scale(Nx1/Nr1); 
  sr1->Draw("L same");


  TH1I* sp1 = (TH1I*) fileS.Get("theta3_P");
  sp1->SetLineColor(4);
  float Np1 = sp1->Integral();

  sp1->Scale(Nx1/Np1); 
  sp1->SetLineWidth(3);
  sp1->Draw("L same");

canvas.cd(2);
  TH2S frame2("frame2","",10,-1,1,10,0,1000);
  frame2.GetXaxis()->SetTitle("cos(#theta_{d-#alpha})");
  frame2.GetYaxis()->SetTitle("Counts");
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.SetStats(kFALSE);
  frame2.Draw();

  text.DrawLatex(-.2,300,"6^{#circ}>#theta_{plf}>10^{#circ}");

  TH1I* histX2 = (TH1I*) file.Get("Li6/theta3_6Li_1_peak2_el");
  histX2->SetMarkerStyle(20);
  histX2->Draw("same P");
  float Nx2 = histX2->Integral();


  TH1I* sr2 = (TH1I*) fileS.Get("theta3_R_peak2");
  sr2->SetLineColor(2);
  float Nr2 = sr2->Integral();

  sr2->Scale(Nx2/Nr2); 
  sr2->Draw("L same");


  TH1I* sp2 = (TH1I*) sp1->Clone();
  sp2->SetLineColor(4);
  float Np2 = sp2->Integral();

  sp2->Scale(Nx2/Np2); 
  sp2->SetLineWidth(3);
  sp2->Draw("L same");

  canvas.cd(3);


  TH2S frame3("frame3","",10,0,360,10,0,1200);
  frame3.GetXaxis()->SetTitle("#phi_{d-#alpha} [deg]");
  frame3.SetStats(kFALSE);
  frame3.GetYaxis()->SetTitle("Counts");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.Draw();

  TH1I* histX3 = (TH1I*) file.Get("Li6/deltaPhi_6Li_1_peak1_el");
  histX3->SetMarkerStyle(20);
  histX3->Draw("same P");
  float Nx3 = histX3->Integral();


  TH1I* sr3 = (TH1I*) fileS.Get("deltaPhi_R_peak1");
  sr3->SetLineColor(2);
  float Nr3 = sr3->Integral();

  sr3->Scale(Nx3/Nr3); 
  sr3->Draw("L same");


  double xx3[90];
  double yy3[90];
  int N3 = 0;;
  for (int i=1;i<90;i++)
    {
      float x = sr3->GetBinCenter(i);


      float y = sr3->GetBinContent(i);

      float y1 = histX3->GetBinContent(i);

      xx3[N3] = x;
      yy3[N3] = y1/y*400;
      N3++;
    }

  TGraph g3(N3,xx3,yy3);
  g3.SetLineColor(4);
  g3.SetLineWidth(2);
  g3.Draw("L");

  canvas.cd(4);


  TH2S frame4("frame4","",10,0,360,10,0,600);
  frame4.GetXaxis()->SetTitle("#phi_{d-#alpha} [deg]");
  frame4.SetStats(kFALSE);
  frame4.GetYaxis()->SetTitle("Counts");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.Draw();

  TH1I* histX4 = (TH1I*) file.Get("Li6/deltaPhi_6Li_1_peak2_el");
  histX4->SetMarkerStyle(20);
  histX4->Draw("same P");
  float Nx4 = histX4->Integral();


  TH1I* sr4 = (TH1I*) fileS.Get("deltaPhi_R_peak2");
  sr4->SetLineColor(2);
  float Nr4 = sr4->Integral();

  sr4->Scale(Nx4/Nr4); 
  sr4->Draw("L same");



  double xx4[90];
  double yy4[90];
  int N4 = 0;;
  for (int i=1;i<90;i++)
    {
      float x = sr4->GetBinCenter(i);


      float y = sr4->GetBinContent(i);

      float y1 = histX4->GetBinContent(i);

      xx4[N4] = x;
      yy4[N4] = y1/y*200;
      N4++;
    }

  TGraph g4(N4,xx4,yy4);
  g4.SetLineColor(4);
  g4.SetLineWidth(2);
  g4.Draw("L");


}
