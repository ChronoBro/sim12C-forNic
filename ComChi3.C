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
Sty->SetPadBottomMargin(0.15);
Sty->SetPadTopMargin(0.05);
Sty->SetPadLeftMargin(0.19);
Sty->SetPadRightMargin(0.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kBlack);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);

Sty->SetLineWidth(3);
Sty->SetLabelSize(0.05,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.05,"xyz");
Sty->SetTitleOffset(2.,"y");
Sty->SetTitleOffset(1.3,"x");
Sty->SetTitleOffset(1.5,"z");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.04,"xyz");
Sty->SetNdivisions(5,"xyz");
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();
  TCanvas canvas("ComChi3","",600,600);

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");


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


  cout << "bad robot " << Nexp2 << " " << Nsim2 << " " << Nsim2_P << endl;

  sle SLE(4);
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

      if (i < 5) cout << i << " " << y << " " << y_plus << " " << y_minus <<
	" " << yexp_minus << " " << ysec_minus << " " << yprim_minus << endl;

      deltaPhi_corr_sub->SetBinContent(i,y_plus-y_minus);
      deltaPhi_corr_sub->SetBinError(i,y_plus/sqrt(yexp_plus)+y_minus/sqrt(yexp_minus));
     
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = cos(3.*x);
      float y3 = cos(5.*x);
      
      SLE.Y[0] += y0*y;
      SLE.Y[1] += y1*y;
      SLE.Y[2] += y2*y;
      SLE.Y[3] += y3*y;

      SLE.M[0][0] += y0*y0;
      SLE.M[0][1] += y0*y1;
      SLE.M[0][2] += y0*y2;
      SLE.M[0][3] += y0*y3;

      SLE.M[1][1] += y1*y1;
      SLE.M[1][2] += y1*y2;
      SLE.M[1][3] += y1*y3;

   
      SLE.M[2][2] += y2*y2;
      SLE.M[2][3] += y2*y3;



      SLE.M[3][3] += y3*y3;
      
    }

  deltaPhi_corr_plus->SetMarkerStyle(20);
  deltaPhi_corr_plus->SetMarkerColor(2);

  //deltaPhi_corr_plus->Draw("Same P");

  deltaPhi_corr_minus->SetMarkerStyle(20);
  deltaPhi_corr_minus->SetMarkerColor(4);


  //deltaPhi_corr_minus->Draw("Same P");


  TH2S frame4("frame4","",10,0,360,10,-1000,1000);
  frame4.GetXaxis()->SetTitle("#chi [deg]");
  frame4.GetYaxis()->SetTitle("W(#psi>0,#chi) - W(#psi<0,#chi)");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.SetStats(kFALSE);
  frame4.Draw();


  deltaPhi_corr_sub->SetMarkerStyle(20);
  deltaPhi_corr_sub->Draw("Same P");


  SLE.M[1][0] = SLE.M[0][1];
  SLE.M[2][0] = SLE.M[0][2];
  SLE.M[3][0] = SLE.M[0][3];


  SLE.M[2][1] = SLE.M[1][2];
  SLE.M[3][1] = SLE.M[1][3];


  SLE.M[3][2] = SLE.M[2][3];

  SLE.solve();

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3]
        << endl;


  float sum=0.;
  for (int i=0;i<4;i++)sum+= SLE.Y[i];


  cout << SLE.Y[0]/SLE.Y[0] << " " << SLE.Y[1]/SLE.Y[0] << " " << SLE.Y[2]/SLE.Y[0]<< " " << SLE.Y[3]/SLE.Y[0]
       <<  endl;


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
      float y2 = cos(3.*x);
      float y3 = cos(5.*x);

      fit->SetBinContent(i,y0*SLE.Y[0]+y1*SLE.Y[1]+y2*SLE.Y[2]+y3*SLE.Y[3]);


      fit1->SetBinContent(i,y1*SLE.Y[1]);
      fit2->SetBinContent(i,y2*SLE.Y[2]);

      fit3->SetBinContent(i,y3*SLE.Y[3]);



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
