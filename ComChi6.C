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
Sty->SetPadTickY(1);
Sty->SetPadBottomMargin(.16);
Sty->SetPadTopMargin(.07);
Sty->SetPadLeftMargin(.2);
Sty->SetPadRightMargin(.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.01,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(1.6,"y");
Sty->SetTitleOffset(1.2,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();
 gStyle->SetErrorX(0.);



  TCanvas canvas("ComChi6_7Be72","",600,900);

 double overlap1 = .04;
  TPad * pad1 = new TPad("pad1","",0.,0.,1.,.5 + overlap1);
  TPad * pad2 = new TPad("pad2","",0.,.5-overlap1,1.,1.);
  pad1->SetFillStyle(4000);
  pad2->SetFillStyle(4000);
  pad1->Draw();
  pad2->Draw();

  pad2->cd();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TFile fileS("sim.root");

  TH2S frame3("frame3","",10,0,360,10,-700,4000);
  //frame3.GetXaxis()->SetTitle("#chi [deg]");
  frame3.GetXaxis()->SetLabelSize(0.);
  frame3.GetYaxis()->SetTitle("W_{even}");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.SetStats(kFALSE);
  frame3.Draw();
  

  TLatex tt;
  tt.SetTextSize(.08);
  tt.SetNDC();
  tt.DrawLatex(.25,.85,"(a)");
  //tt.DrawLatex(50,3500,"^{7}Be#rightarrow ^{3}He+#alpha");
  //tt.DrawLatex(250,3500,"j=7/2^{-}, l=3");


  TH1I* deltaPhi = (TH1I*) file.Get("deltaPhi_7Be_2_el");
  deltaPhi->SetMarkerStyle(20);
  //deltaPhi->Draw("same P");
  float Nexp2 = deltaPhi->Integral();

  TH1I* deltaPhi_sim = (TH1I*) fileS.Get("deltaPhi_R");
  float Nsim2 = deltaPhi_sim->Integral();
  deltaPhi_sim->Scale(Nexp2/Nsim2);
  deltaPhi_sim->SetLineColor(2);
  //deltaPhi_sim->Draw("L same");



  TH1I* deltaPhi_sim_P = (TH1I*) fileS.Get("deltaPhi_P");
  float Nsim2_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P->Scale(Nexp2/Nsim2_P);
  deltaPhi_sim_P->SetLineColor(4);
  deltaPhi_sim_P->SetLineWidth(3);
  //deltaPhi_sim_P->Draw("same L");

  sle SLE(5);
  SLE.clear();

  inverse Inverse(5);

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
      float y1 = cos(x);
      float y2 = cos(2.*x);
      float y3 = cos(4.*x);
      float y4 = cos(6.*x);
      
      SLE.Y[0] += y0*y;
      SLE.Y[1] += y1*y;
      SLE.Y[2] += y2*y;
      SLE.Y[3] += y3*y;
      SLE.Y[4] += y4*y;

      SLE.M[0][0] += y0*y0;
      SLE.M[0][1] += y0*y1;
      SLE.M[0][2] += y0*y2;
      SLE.M[0][3] += y0*y3;
      SLE.M[0][4] += y0*y4;


      SLE.M[1][1] += y1*y1;
      SLE.M[1][2] += y1*y2;
      SLE.M[1][3] += y1*y3;
      SLE.M[1][4] += y1*y4;
   
      SLE.M[2][2] += y2*y2;
      SLE.M[2][3] += y2*y3;
      SLE.M[2][4] += y2*y4;



      SLE.M[3][3] += y3*y3;
      SLE.M[3][4] += y3*y4;

      SLE.M[4][4] += y4*y4;

    }

  SLE.M[1][0] = SLE.M[0][1];
  SLE.M[2][0] = SLE.M[0][2];
  SLE.M[3][0] = SLE.M[0][3];
  SLE.M[4][0] = SLE.M[0][4];



  SLE.M[2][1] = SLE.M[1][2];
  SLE.M[3][1] = SLE.M[1][3];
  SLE.M[4][1] = SLE.M[1][4];


  SLE.M[3][2] = SLE.M[2][3];
  SLE.M[4][2] = SLE.M[2][4];

  SLE.M[4][3] = SLE.M[3][4];

  for (int i=0;i<5;i++)
    for (int j=0;j<5;j++)
      Inverse.A[i][j] = SLE.M[i][j];

  SLE.solve();
  Inverse.solve();

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << " " << SLE.Y[3]
       << " " << SLE.Y[4] << endl;


  float sum=0.;
  for (int i=0;i<5;i++)sum+= SLE.Y[i];




  deltaPhi_corr->SetLineColor(2);
  deltaPhi_corr->SetMarkerStyle(20);
  deltaPhi_corr->SetMarkerSize(1.5);
  deltaPhi_corr->Draw("same P");

  TH1I* fit = (TH1I*) deltaPhi->Clone();
  TH1I* fit1 = (TH1I*) fit->Clone();
  TH1I* fit2 = (TH1I*) fit->Clone();
  TH1I* fit3 = (TH1I*) fit->Clone();
  TH1I* fit4 = (TH1I*) fit->Clone();

  for (int i=1;i<=90;i++)

    {
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = cos(2.*x);
      float y3 = cos(4.*x);
      float y4 = cos(6.*x);

      fit->SetBinContent(i,y0*SLE.Y[0]+y1*SLE.Y[1]+y2*SLE.Y[2]+y3*SLE.Y[3]
+y4*SLE.Y[4]);


      fit1->SetBinContent(i,y1*SLE.Y[1]);
      fit2->SetBinContent(i,y2*SLE.Y[2]);

      fit3->SetBinContent(i,y3*SLE.Y[3]);
      fit4->SetBinContent(i,y4*SLE.Y[4]);



    }

  fit->SetLineWidth(3);
  fit->SetLineColor(2);
  fit->Draw("same C");

  fit1->SetLineColor(4);
  fit1->SetLineStyle(2);
  fit1->Draw("same C");
  fit2->SetLineColor(6);
  fit2->SetLineStyle(2);
  fit2->Draw("same C");
  fit3->SetLineColor(3);
  fit3->SetLineStyle(2);
  fit3->Draw("same C");

  fit4->SetLineColor(5);
  fit4->SetLineStyle(2);
  fit4->Draw("same C");

    TLine ll;
    ll.SetLineStyle(9);
    ll.SetLineWidth(3);

    ll.DrawLine(0.,0.,360.,0.);


    ll.SetLineColor(8);
    ll.SetLineWidth(3);
    ll.SetLineStyle(2);
    ll.DrawLine(0.,SLE.Y[0],360.,SLE.Y[0]);

    pad1->cd();


  TH1I* deltaPhi_plus = (TH1I*) file.Get("deltaPhi_7Be_2_el_plus");
  TH1I* deltaPhi_minus = (TH1I*) file.Get("deltaPhi_7Be_2_el_minus");





  TH1I* deltaPhi_sim_plus = (TH1I*) fileS.Get("deltaPhi_plus_R");
  TH1I* deltaPhi_sim_minus = (TH1I*) fileS.Get("deltaPhi_minus_R");

  deltaPhi_sim_plus->Scale(Nexp2/Nsim2);
  deltaPhi_sim_minus->Scale(Nexp2/Nsim2);


  TH1I* deltaPhi_sim_P_plus = (TH1I*) fileS.Get("deltaPhi_plus_P");
  TH1I* deltaPhi_sim_P_minus = (TH1I*) fileS.Get("deltaPhi_minus_P");
  float Nsim22_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P_plus->Scale(Nexp2/Nsim2_P);
  deltaPhi_sim_P_minus->Scale(Nexp2/Nsim2_P);




  sle SLEE(3);
  SLEE.clear();

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
      deltaPhi_corr_sub->SetBinError(i,y_plus/sqrt(yexp_plus)+y_minus/sqrt(yexp_minus));
     
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = cos(3.*x);
      float y3 = cos(5.*x);
      

      SLEE.Y[0] += y1*y;
      SLEE.Y[1] += y2*y;
      SLEE.Y[2] += y3*y;

      SLEE.M[0][0] += y1*y1;
      SLEE.M[0][1] += y1*y2;
      SLEE.M[0][2] += y1*y3;

   
      SLEE.M[1][1] += y2*y2;
      SLEE.M[1][2] += y2*y3;



      SLEE.M[2][2] += y3*y3;
      
    }

  deltaPhi_corr_plus->SetMarkerStyle(20);
  deltaPhi_corr_plus->SetMarkerColor(2);

  //deltaPhi_corr_plus->Draw("Same P");

  deltaPhi_corr_minus->SetMarkerStyle(20);
  deltaPhi_corr_minus->SetMarkerColor(4);


  //deltaPhi_corr_minus->Draw("Same P");


  TH2S frame4("frame4","",10,0,360,10,-1000,1000);
  frame4.GetXaxis()->SetTitle("#chi [deg]");
  //frame4.GetYaxis()->SetTitle("W(#psi>0,#chi) - W(#psi<0,#chi)");
  frame4.GetYaxis()->SetTitle("W_{odd}");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.SetStats(kFALSE);
  frame4.Draw();


  tt.DrawLatex(.25,.85,"(b)");

  deltaPhi_corr_sub->SetMarkerStyle(20);
  deltaPhi_corr_sub->SetLineColor(1);
  deltaPhi_corr_sub->Draw("Same P");


  SLEE.M[1][0] = SLEE.M[0][1];
  SLEE.M[2][0] = SLEE.M[0][2];
  SLEE.M[2][1] = SLEE.M[1][2];




  SLEE.solve();

  cout << "solve" << endl;
  cout << SLEE.Y[0] << " " << SLEE.Y[1] << " " << SLEE.Y[2] << endl;

  cout << "hen" << endl;
  sum=0.;
  for (int i=0;i<3;i++)sum+= SLEE.Y[i];


  cout << SLEE.Y[0]/SLEE.Y[0] << " " << SLEE.Y[1]/SLEE.Y[0] << " " << SLEE.Y[2]/SLEE.Y[0]<<  endl;

  cout << "cow"<< endl;

  TH1I* ffit = (TH1I*) deltaPhi->Clone();
  TH1I* ffit1 = (TH1I*) ffit->Clone();
  TH1I* ffit2 = (TH1I*) ffit->Clone();
  TH1I* ffit3 = (TH1I*) ffit->Clone();

  for (int i=1;i<=90;i++)

    {
      float x = deltaPhi->GetBinCenter(i);
      x = x/180.*3.14159;

      float y0 = 1.;
      float y1 = cos(x);
      float y2 = cos(3.*x);
      float y3 = cos(5.*x);

      ffit->SetBinContent(i,y1*SLEE.Y[0]+y2*SLEE.Y[1]+y3*SLEE.Y[2]);


      ffit1->SetBinContent(i,y1*SLEE.Y[0]);
      ffit2->SetBinContent(i,y2*SLEE.Y[1]);

      ffit3->SetBinContent(i,y3*SLEE.Y[2]);



    }

  ffit->SetLineWidth(3);
  ffit->SetLineColor(2);
  ffit->Draw("same C");

  ffit1->SetLineColor(4);
  ffit1->SetLineStyle(2);
  ffit1->Draw("same C");
  ffit2->SetLineColor(6);
  ffit2->SetLineStyle(2);
  ffit2->Draw("same C");
  ffit3->SetLineColor(3);
  ffit3->SetLineStyle(2);
  ffit3->Draw("same C");

  //    TLine ll;
  ll.SetLineStyle(9);
    ll.SetLineWidth(3);
    ll.SetLineColor(1);
    ll.DrawLine(0.,0.,360.,0.);


    ll.SetLineColor(1);
    ll.SetLineStyle(2);
    ll.SetLineWidth(3);
  
    //ll.DrawLine(0.,SLEE.Y[0],360.,SLEE.Y[0]);




}
