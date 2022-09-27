{
  gROOT->Reset();

  TCanvas canvas("ComChi","",600,900);
  canvas.Divide(1,1);
  canvas->cd(1);

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

  TH2S frame3("frame3","",10,0,360,10,-700,4000);
  frame3.GetXaxis()->SetTitle("#chi [deg]");
  frame3.GetYaxis()->SetTitle("Counts");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.SetStats(kFALSE);
  frame3.Draw();
  

  TLatex tt;
  tt.DrawLatex(50,3500,"^{7}Be#rightarrow ^{3}He+#alpha");
  tt.DrawLatex(250,3500,"j=7/2^{-}, l=3");


  TH1I* deltaPhi = (TH1I*) file.Get("deltaPhi_7Be_2_el_midAngle");
  deltaPhi->SetMarkerStyle(20);
  //deltaPhi->Draw("same P");
  float Nexp2 = deltaPhi->Integral();

  TH1I* deltaPhi_sim = (TH1I*) fileS.Get("deltaPhi_midAngle_R");
  float Nsim2 = deltaPhi_sim->Integral();
  deltaPhi_sim->Scale(Nexp2/Nsim2);
  deltaPhi_sim->SetLineColor(2);
  //deltaPhi_sim->Draw("L same");



  TH1I* deltaPhi_sim_P = (TH1I*) fileS.Get("deltaPhi_midAngle_P");
  float Nsim2_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P->Scale(Nexp2/Nsim2_P);
  deltaPhi_sim_P->SetLineColor(4);
  deltaPhi_sim_P->SetLineWidth(3);
  //deltaPhi_sim_P->Draw("same L");

  sle SLE(7);
  SLE.clear();

  TH1I* deltaPhi_corr = (TH1I*)deltaPhi_sim->Clone();
  for (int i=1;i<=90;i++)

    {
      float yexp = deltaPhi->GetBinContent(i);
      float ysec = deltaPhi_sim->GetBinContent(i);
      float yprim = deltaPhi_sim_P->GetBinContent(i);

      cout << yexp << " " << ysec << " " << yprim << endl;

      float y = yexp/ysec*yprim;
      deltaPhi_corr->SetBinContent(i,y);
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(2.*x);
      float y2 = sin(2.*x);
      float y3 = cos(4.*x);
      float y4 = sin(4.*x);
      float y5 = cos(6.*x);
      float y6 = sin(6.*x);
      
      SLE.Y[0] += y0*y;
      SLE.Y[1] += y1*y;
      SLE.Y[2] += y2*y;
      SLE.Y[3] += y3*y;
      SLE.Y[4] += y4*y;
      SLE.Y[5] += y5*y;
      SLE.Y[6] += y6*y;

      SLE.M[0][0] += y0*y0;
      SLE.M[0][1] += y0*y1;
      SLE.M[0][2] += y0*y2;
      SLE.M[0][3] += y0*y3;
      SLE.M[0][4] += y0*y4;
      SLE.M[0][5] += y0*y5;
      SLE.M[0][6] += y0*y6;

      SLE.M[1][1] += y1*y1;
      SLE.M[1][2] += y1*y2;
      SLE.M[1][3] += y1*y3;
      SLE.M[1][4] += y1*y4;
      SLE.M[1][5] += y1*y5;
      SLE.M[1][6] += y1*y6;
   
      SLE.M[2][2] += y2*y2;
      SLE.M[2][3] += y2*y3;
      SLE.M[2][4] += y2*y4;
      SLE.M[2][5] += y2*y5;
      SLE.M[2][6] += y2*y6;


      SLE.M[3][3] += y3*y3;
      SLE.M[3][4] += y3*y4;
      SLE.M[3][5] += y3*y5;
      SLE.M[3][6] += y3*y6;


      SLE.M[4][4] += y4*y4;
      SLE.M[4][5] += y4*y5;
      SLE.M[4][6] += y4*y6;

      SLE.M[5][5] += y5*y5;
      SLE.M[5][5] += y5*y6;

      SLE.M[6][6] += y6*y6;

    }

  SLE.M[1][0] = SLE.M[0][1];
  SLE.M[2][0] = SLE.M[0][2];
  SLE.M[3][0] = SLE.M[0][3];
  SLE.M[4][0] = SLE.M[0][4];
  SLE.M[5][0] = SLE.M[0][5];
  SLE.M[6][0] = SLE.M[0][6];


  SLE.M[2][1] = SLE.M[1][2];
  SLE.M[3][1] = SLE.M[1][3];
  SLE.M[4][1] = SLE.M[1][4];
  SLE.M[5][1] = SLE.M[1][5];
  SLE.M[6][1] = SLE.M[1][6];


  SLE.M[3][2] = SLE.M[2][3];
  SLE.M[4][2] = SLE.M[2][4];
  SLE.M[5][2] = SLE.M[2][5];
  SLE.M[6][2] = SLE.M[2][6];


  SLE.M[4][3] = SLE.M[3][4];
  SLE.M[5][3] = SLE.M[3][5];
  SLE.M[6][3] = SLE.M[3][6];


  SLE.M[5][4] = SLE.M[4][5];
  SLE.M[6][4] = SLE.M[4][6];

  SLE.M[6][5] = SLE.M[5][6];

  SLE.solve();

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3]
       << " " << SLE.Y[4] << " " << SLE.Y[5] << " " << SLE.Y[6] << endl;


  float sum=0.;
  for (int i=0;i<7;i++)sum+= SLE.Y[i];


  cout << SLE.Y[0]/SLE.Y[0] << " " << SLE.Y[1]/SLE.Y[0] << " " << SLE.Y[2]/SLE.Y[0]<< " " << SLE.Y[3]/SLE.Y[0]
       << " " << SLE.Y[4]/SLE.Y[0] << " " << SLE.Y[5]/SLE.Y[0] << " " << SLE.Y[6]/SLE.Y[0] << endl;

  deltaPhi_corr->SetLineColor(2);
  deltaPhi_corr->SetMarkerStyle(20);
  deltaPhi_corr->Draw("same P");

  TH1I* fit = (TH1I*) deltaPhi->Clone();
  TH1I* fit1 = (TH1I*) fit->Clone();
  TH1I* fit2 = (TH1I*) fit->Clone();
  TH1I* fit3 = (TH1I*) fit->Clone();

  for (int i=1;i<=90;i++)

    {
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(2.*x);
      float y2 = sin(2.*x);
      float y3 = cos(4.*x);
      float y4 = sin(4.*x);
      float y5 = cos(6.*x);
      float y6 = sin(6.*x);

      fit->SetBinContent(i,y0*SLE.Y[0]+y1*SLE.Y[1]+y2*SLE.Y[2]+y3*SLE.Y[3]
      +y4*SLE.Y[4] + y5*SLE.Y[5] + y6*SLE.Y[6]);

      fit1->SetBinContent(i,y1*SLE.Y[1]+y2*SLE.Y[2]);
      fit2->SetBinContent(i,y3*SLE.Y[3]+y4*SLE.Y[4]);

      fit3->SetBinContent(i,y5*SLE.Y[5]+y6*SLE.Y[6]);



    }

  fit->SetLineWidth(3);
  fit->SetLineColor(2);
  fit->Draw("same C");

  fit1->SetLineColor(4);
  fit1->Draw("same C");
  fit2->SetLineColor(6);
  fit2->Draw("same C");
  fit3->SetLineColor(3);
  fit3->Draw("same C");

    TLine ll;
  ll.SetLineStyle(9);
    ll.SetLineWidth(3);
    ll.DrawLine(0.,0.,360.,0.);


    ll.SetLineColor(7);
    ll.SetLineWidth(1);
    ll.DrawLine(0.,SLE.Y[0],360.,SLE.Y[0]);

}
