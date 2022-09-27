{
  gROOT->Reset();

  TCanvas canvas("ComChi","",600,900);
  canvas.Divide(1,2);
  canvas->cd(1);

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

  TH2S frame3("frame3","",10,0,360,10,0,2500);
  frame3.GetXaxis()->SetTitle("#chi [deg]");
  frame3.GetYaxis()->SetTitle("Counts");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.SetStats(kFALSE);
  frame3.Draw();
  

  TLatex tt;
  tt.DrawLatex(50,2300,"^{7}Be#rightarrow ^{3}He+#alpha");
  tt.DrawLatex(250,2300,"j=7/2^{-}, l=3");

  tt.SetTextColor(2);
  tt.DrawLatex(210,1700,"Y(#psi>0)");
  tt.SetTextColor(4);
  tt.DrawLatex(160,900,"Y(#psi<0)");


  TH1I* deltaPhi_plus = (TH1I*) file.Get("deltaPhi_7Be_2_el_plus");
  TH1I* deltaPhi = (TH1I*) file.Get("deltaPhi_7Be_2_el");
  TH1I* deltaPhi_minus = (TH1I*) file.Get("deltaPhi_7Be_2_el_minus");
  float Nexp2 = deltaPhi->Integral();



  TH1I* deltaPhi_sim = (TH1I*) fileS.Get("deltaPhi_R");
  TH1I* deltaPhi_sim_plus = (TH1I*) fileS.Get("deltaPhi_plus_R");
  TH1I* deltaPhi_sim_minus = (TH1I*) fileS.Get("deltaPhi_minus_R");
  float Nsim2 = deltaPhi_sim->Integral();
  deltaPhi_sim_plus->Scale(Nexp2/Nsim2);
  deltaPhi_sim_minus->Scale(Nexp2/Nsim2);

  TH1I* deltaPhi_sim_P = (TH1I*) fileS.Get("deltaPhi_P");
  TH1I* deltaPhi_sim_P_plus = (TH1I*) fileS.Get("deltaPhi_plus_P");
  TH1I* deltaPhi_sim_P_minus = (TH1I*) fileS.Get("deltaPhi_minus_P");
  float Nsim2_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P_plus->Scale(Nexp2/Nsim2_P);
  deltaPhi_sim_P_minus->Scale(Nexp2/Nsim2_P);

  sle SLE(7);
  SLE.clear();

  TH1I* deltaPhi_corr_plus = (TH1I*)deltaPhi_sim->Clone();
  TH1I* deltaPhi_corr_minus = (TH1I*)deltaPhi_sim->Clone();
  TH1I* deltaPhi_corr_sub = (TH1I*)deltaPhi_sim->Clone();
  for (int i=1;i<=90;i++)

    {
      float yexp_plus = deltaPhi_plus->GetBinContent(i);
      float ysec_plus = deltaPhi_sim_plus->GetBinContent(i);
      float yprim_plus = deltaPhi_sim_P_plus->GetBinContent(i);


      float y_plus = yexp_plus/ysec_plus*yprim_plus;
      deltaPhi_corr_plus->SetBinContent(i,y_plus);


      float yexp_minus = deltaPhi_minus->GetBinContent(i);
      float ysec_minus = deltaPhi_sim_minus->GetBinContent(i);
      float yprim_minus = deltaPhi_sim_P_minus->GetBinContent(i);

      float y_minus = yexp_minus/ysec_minus*yprim_minus;
      deltaPhi_corr_minus->SetBinContent(i,y_minus);

      float y = y_plus-y_minus;

      deltaPhi_corr_sub->SetBinContent(i,y_plus-y_minus);

     
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = sin(x);
      float y3 = cos(3.*x);
      float y4 = sin(3.*x);
      float y5 = cos(5.*x);
      float y6 = sin(5.*x);
      
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

  deltaPhi_corr_plus->SetMarkerStyle(20);
  deltaPhi_corr_plus->SetMarkerColor(2);

  deltaPhi_corr_plus->Draw("Same P");

  deltaPhi_corr_minus->SetMarkerStyle(20);
  deltaPhi_corr_minus->SetMarkerColor(4);


  deltaPhi_corr_minus->Draw("Same P");

  canvas.cd(2);
  TH2S frame4("frame4","",10,0,360,10,-1000,1000);
  frame4.GetXaxis()->SetTitle("#chi [deg]");
  frame4.GetYaxis()->SetTitle("Y(#psi>0) - Y(#psi<0)");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.SetStats(kFALSE);
  frame4.Draw();

  cout << "her" << endl;
  deltaPhi_corr_sub->SetMarkerStyle(20);
  deltaPhi_corr_sub->Draw("Same P");


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


  TH1I* fit = (TH1I*) deltaPhi->Clone();
  TH1I* fit1 = (TH1I*) fit->Clone();
  TH1I* fit2 = (TH1I*) fit->Clone();
  TH1I* fit3 = (TH1I*) fit->Clone();

  for (int i=1;i<=90;i++)

    {
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = sin(x);
      float y3 = cos(3.*x);
      float y4 = sin(3.*x);
      float y5 = cos(5.*x);
      float y6 = sin(5.*x);

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
