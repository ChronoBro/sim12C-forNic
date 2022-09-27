{
  gROOT->Reset();


  TCanvas canvas("ComTheta","",600,900);
  canvas.Divide(1,3);
  canvas->cd(1);

  TH2S frame("frame","",10,0,14,10,0,10000);
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


  canvas.cd(2);
  TH2S frame2("frame2","",10,-1,1,10,0,3000);
  frame2.GetXaxis()->SetTitle("cos(#psi_{^{2}H})");
  frame2.SetStats(kFALSE);
  frame2.Draw();


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

  canvas.cd(3);
  TH2S frame3("frame3","",10,0,360,0,10,1800);
  frame3.GetXaxis()->SetTitle("#phi_{d} [deg]");
  frame3.GetYaxis()->SetTitle("Counts");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.SetStats(kFALSE);
  frame3.Draw();
  

  TH1I* deltaPhi = (TH1I*) file.Get("Be7/deltaPhi_7Be_2_el");
  deltaPhi->SetMarkerStyle(20);
  deltaPhi->Draw("same P");
  float Nexp2 = deltaPhi->Integral();

  TH1I* deltaPhi_sim = (TH1I*) fileS.Get("deltaPhi_R");
  float Nsim2 = deltaPhi_sim->Integral();
  deltaPhi_sim->Scale(Nexp2/Nsim2);
  deltaPhi_sim->SetLineColor(2);
  deltaPhi_sim->Draw("L same");


  
  TH1I* deltaPhi_sim_P = (TH1I*) fileS.Get("deltaPhi_P");
  float Nsim2_P = deltaPhi_sim_P->Integral();
  deltaPhi_sim_P->Scale(Nexp2/Nsim2_P);
  deltaPhi_sim_P->SetLineColor(4);
  deltaPhi_sim_P->SetLineWidth(2);
  deltaPhi_sim_P->Draw("L same");
  
  return;



  TH1I* deltaPhi_sim_P = (TH1I*)deltaPhi_sim->Clone();
  TH1I* prim = (TH1I*)deltaPhi_sim->Clone();
  for (int i=1;i<=90;i++)
    {
      deltaPhi_sim_P->SetBinContent(i,Nexp2/90.);
      prim->SetBinContent(i,1000.*deltaPhi->GetBinContent(i)/deltaPhi_sim->GetBinContent(i));
    }
  deltaPhi_sim_P->SetLineColor(4);
  deltaPhi_sim_P->SetLineWidth(3);
  deltaPhi_sim_P->Draw("same L");

  prim->SetLineColor(6);
  prim->Draw("same L");

  TF1 * funct = new TF1("funct","[0] + [1]*exp(-pow(x/[2],2)/2.) + [1]*exp(-pow((x-360)/[2],2)/2.) + [3]*exp(-pow((x-180)/[4],2)/2.) + [5]*exp(-pow((x-90)/[6],2)/2.)+[5]*exp(-pow((x-270)/[6],2)/2.)",0,360);
  funct->SetParameter(0,1000);
  funct->SetParameter(1,400);
  funct->SetParameter(2,20);
  funct->SetParameter(3,200);
  funct->SetParameter(4,20);
  funct->SetParameter(5,100);
  funct->SetParameter(6,20);

  prim->Fit(funct);

}
