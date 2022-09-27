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
Sty->SetPadBottomMargin(.25);
Sty->SetPadTopMargin(.05);
Sty->SetPadLeftMargin(.2);
Sty->SetPadRightMargin(.05);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.07,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.07,"xyz");
Sty->SetTitleOffset(2,"y");
Sty->SetTitleOffset(1.1,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();
  TCanvas canvas("ComTheta","",600,900);
  canvas.Divide(2,4);
  canvas->cd(1);

  TH2S frame("frame","",10,0,14,10,0,28000);
  frame.GetXaxis()->SetTitle("#theta_{lab} [deg]");
  frame.SetStats(kFALSE);
  frame.Draw();

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH1I* histX = (TH1I*) file.Get("Be7/theta_7Be_2_el");
  histX->SetMarkerStyle(20);
  histX->Draw("same P");
  float Nx = histX->Integral();

  TFile fileS("sim.root");
  TH1I* sr = (TH1I*) fileS.Get("thetaR");
  sr->SetLineColor(2);
  float Nr = sr->Integral();

  sr->Scale(Nx/Nr); 
  sr->Draw("L same");



  TH1I* sp = (TH1I*) fileS.Get("theta");
  sp->SetLineColor(4);
  float Np = sp->Integral();

  sp->Scale(Nx/Np);
  sp->SetLineWidth(3); 
  sp->Draw("L same");


  canvas.cd(3);
  TH2S frame3("frame3","",10,-1,1,10,0,10000);
  frame3.GetXaxis()->SetTitle("cos(#psi_{Beam})");
  frame3.SetStats(kFALSE);
  frame3.Draw();


  TH1I* histX3 = (TH1I*) file.Get("Be7/theta3_7Be_2_el");
  histX3->SetMarkerStyle(20);
  histX3->Draw("same P");
  float Nx3 = histX3->Integral();


  TH1I* sr3 = (TH1I*) fileS.Get("theta3_R");
  sr3->SetLineColor(2);
  float Nr3 = sr3->Integral();

  sr3->Scale(Nx3/Nr3); 
  sr3->Draw("L same");


  TH1I* sp3 = (TH1I*) fileS.Get("theta3_P");
  sp3->SetLineColor(4);
  float Np3 = sp3->Integral();

  sp3->Scale(Nx3/Np3); 
  sp3->SetLineWidth(3);
  sp3->Draw("L same");

  /*
  TH1F* prim = new TH1F("prim","",50,-1,1);

  TH1F* p22 = new TH1F("p22","",50,-1,1);
  TH1F* p21 = new TH1F("p21","",50,-1,1);
  TH1F* p20 = new TH1F("p20","",50,-1,1);
  TH1F* pp = new TH1F("pp","",50,-1,1);

  sle SLE(3);
  SLE.clear();

  


  for (int i=1;i<=50;i++)
    {
      float yexp = histX3->GetBinContent(i);
      float ysim = sr3->GetBinContent(i);
      float y = 2000.*yexp/ysim;
      prim->SetBinContent(i,y);
      float cosx = histX3->GetBinCenter(i);
      float theta = acos(cosx);
      float sinx = sin(theta);

      float y20 = pow((3.*pow(cosx,2)-1.)/2.,2);
      float y21 = pow(3.*sinx*cosx,2);
      float y22 = pow(3.*sinx*sinx,2);

      y20 /= 2./5.;
      y21 /= 12./5.;
      y22 /= 48./5.;

      SLE.Y[0] += y20*y;
      SLE.Y[1] += y21*y;
      SLE.Y[2] += y22*y;

      SLE.M[0][0] += y20*y20;
      SLE.M[1][0] += y20*y21;
      SLE.M[2][0] += y20*y22;
      SLE.M[1][1] += y21*y21;
      SLE.M[2][1] += y21*y22;
      SLE.M[2][2] += y22*y22;


      p22->SetBinContent(i,y22*2784);
      p21->SetBinContent(i,y21*901.);
      p20->SetBinContent(i,y20*263);
      pp->SetBinContent(i,263*y20+901.*y21+2784.*y22);
    }
  SLE.M[0][1] = SLE.M[1][0];
  SLE.M[0][2] = SLE.M[2][0];
  SLE.M[1][2] = SLE.M[2][1];
  SLE.solve();

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << endl;
  prim->SetLineColor(6);
  prim->Draw("same L");

  //p22->Draw("same L");
  //p21->Draw("same L");
  //p20->Draw("same L");
  pp->SetLineColor(3);
  pp->Draw("same L");
  prim->Draw("same");

  */

  canvas.cd(4);
  TH2S frame4("frame4","",10,0,360,0,10,5000);
  frame4.GetXaxis()->SetTitle("#phi_{Beam} [deg]");
  frame4.GetYaxis()->SetTitle("Counts");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.SetStats(kFALSE);
  frame4.Draw();
  

  TH1I* deltaPhi = (TH1I*) file.Get("Be7/deltaPhi_7Be_2_el");
  deltaPhi->SetMarkerStyle(20);
  deltaPhi->Draw("same P");
  float Nexp4 = deltaPhi->Integral();

  TH1I* deltaPhi_sim = (TH1I*) fileS.Get("deltaPhi_R");
  float Nsim4 = deltaPhi_sim->Integral();
  deltaPhi_sim->Scale(Nexp4/Nsim4);
  deltaPhi_sim->SetLineColor(2);
  deltaPhi_sim->Draw("L same");


  
  TH1I* deltaPhi_sim_P = (TH1I*) fileS.Get("deltaPhi_P");
  float Nsim4_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P->Scale(Nexp4/Nsim4_P);
  deltaPhi_sim_P->SetLineColor(4);
  deltaPhi_sim_P->SetLineWidth(2);
  deltaPhi_sim_P->Draw("L same");



  canvas.cd(5);
  TH2S frame5("frame5","",10,-1,1,10,0,10000);
  frame5.GetXaxis()->SetTitle("cos(#psi_{Q})");
  frame5.SetStats(kFALSE);
  frame5.Draw();


  TH1I* histX5 = (TH1I*) file.Get("Be7/theta3Q_7Be_2_el");
  histX5->SetMarkerStyle(20);
  histX5->Draw("same P");
  float Nx5 = histX5->Integral();


  TH1I* sr5 = (TH1I*) fileS.Get("cosPsiQ_R");
  sr5->SetLineColor(2);
  float Nr5 = sr5->Integral();

  sr5->Scale(Nx5/Nr5); 
  sr5->Draw("L same");


  TH1I* sp5 = (TH1I*) fileS.Get("cosPsiQ_P");
  sp5->SetLineColor(4);
  float Np5 = sp5->Integral();

  sp5->Scale(Nx5/Np5); 
  sp5->Draw("L same");


  canvas.cd(6);


  TH2S frame6("frame6","",10,0,360,10,0,6000);
  frame6.GetXaxis()->SetTitle("#chi_{Q} [deg]");
  frame6.GetYaxis()->SetTitle("Counts");
  frame6.GetXaxis()->CenterTitle();
  frame6.GetYaxis()->CenterTitle();
  frame6.SetStats(kFALSE);
  frame6.Draw();



  TH1I* histX6 = (TH1I*) file.Get("Be7/chiQ_7Be_2_el");
  histX6->SetMarkerStyle(20);
  histX6->Draw("same P");
  float Nx6 = histX6->Integral();


  TH1I* sr6 = (TH1I*) fileS.Get("chiQ_R");
  sr6->SetLineColor(2);
  float Nr6 = sr6->Integral();

  sr6->Scale(Nx6/Nr6); 
  sr6->Draw("L same");


  TH1I* sp6 = (TH1I*) fileS.Get("chiQ_P");
  sp6->SetLineColor(4);
  float Np6 = sp6->Integral();

  sp6->Scale(Nx6/Np6); 
  sp6->Draw("L same");



  canvas.cd(7);


  TH2S frame7("frame7","",10,-1,1,10,0,10000);
  frame7.GetXaxis()->SetTitle("cos(#psi_{QP}) [deg]");
  frame7.GetYaxis()->SetTitle("Counts");
  frame7.GetXaxis()->CenterTitle();
  frame7.GetYaxis()->CenterTitle();
  frame7.SetStats(kFALSE);
  frame7.Draw();



  TH1I* histX7 = (TH1I*) file.Get("Be7/theta3QP_7Be_2_el");
  histX7->SetMarkerStyle(20);
  histX7->Draw("same P");
  float Nx7 = histX7->Integral();


  TH1I* sr7 = (TH1I*) fileS.Get("cosPsiQP_R");
  sr7->SetLineColor(2);
  float Nr7 = sr7->Integral();

  sr7->Scale(Nx7/Nr7); 
  sr7->Draw("L same");


  TH1I* sp7 = (TH1I*) fileS.Get("cosPsiQP_P");
  sp7->SetLineColor(4);
  float Np7 = sp7->Integral();

  sp7->Scale(Nx7/Np7); 
  sp7->Draw("L same");

  TH1I* spp7 = (TH1I*)histX7->Clone();
  for (int i=1;i<=50;i++)
    {
      float y = histX7->GetBinContent(i);
      spp7->SetBinContent(i,y/sr7->GetBinContent(i)*sp7->GetBinContent(i));
    }

	spp7->SetLineColor(6);
	spp7->SetLineWidth(3);
      spp7->Draw("same L");

      cout << "RMS = " << spp7->GetRMS() << endl;

  canvas.cd(8);
  TH2S frame8("frame8","",10,0,360,10,0,6000);
  frame8.GetXaxis()->SetTitle("#chi_{Q} [deg]");
  frame8.GetYaxis()->SetTitle("Counts");
  frame8.GetXaxis()->CenterTitle();
  frame8.GetYaxis()->CenterTitle();
  frame8.SetStats(kFALSE);
  frame8.Draw();



  TH1I* histX8 = (TH1I*) file.Get("Be7/chiQP_7Be_2_el");
  histX8->SetMarkerStyle(20);
  histX8->Draw("same P");
  float Nx8 = histX8->Integral();

  
  TH1I* sr8 = (TH1I*) fileS.Get("chiQP_R");
  sr8->SetLineColor(2);
  float Nr8 = sr8->Integral();

  sr8->Scale(Nx8/Nr8); 
  sr8->Draw("L same");


  TH1I* sp8 = (TH1I*) fileS.Get("chiQP_P");
  sp8->SetLineColor(4);
  float Np8 = sp8->Integral();

  sp8->Scale(Nx8/Np8); 
  sp8->Draw("L same");


}
