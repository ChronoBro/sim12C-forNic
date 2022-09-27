#include "Profile.h"
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH2S.h"
#include "TH1I.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "coul.h"
#include "sle.h"
int main()
{
  const double temp = 3.2;

  coul Coul;
  bool me = 1;
  Profile * profile1 = new Profile();
  Profile * profile2 = new Profile();
  double r0 = 1.45;

  //first resonance
  int l0 = 1;
  int z1_0 = 2;
  int z2_0 = 2;  
  double a0 = r0*(pow(3.,1./3.)+pow(3.,1./3.));
  double mu0 = 3.*3./6.;
  double gamma0 = .84;
  double Er0 = 2.54;
  double P0 = profile1->P_l(Er0,a0,l0,mu0,z1_0,z2_0); 
   double width = 2.*gamma0*P0;
  cout << "width= " << width << endl;
  double dE = .02;
  double shift11 = profile1->Shift(Er0-dE,a0,l0,mu0,z1_0,z2_0); 
  double shift12 = profile1->Shift(Er0+dE,a0,l0,mu0,z1_0,z2_0);
  double dShift = (shift12-shift11)/2./dE;
  double shift1 = profile1->Shift(Er0,a0,l0,mu0,z1_0,z2_0);

  cout << "P0 = " << P0 << endl;
  cout << "dS0/dE = " << dShift << endl; 
  cout << "corred  width= " << width/(1.+gamma0*dShift) << endl;   
  double gammaW0 = 62.34/mu0/pow(a0,2);
  cout << "Wig limit = " << gammaW0 << endl; 

  profile1->Er = Er0;
  profile1->gamma = gamma0;
  profile1->l = l0;

  //second resonance
  gamma0 = 1.96;
  Er0 = 5.25;
  l0 = 0;
  //double shift21 = profile2->Shift(Er0-dE,a0,l0,mu0,z1_0,z2_0); 


  profile2->Er = Er0;
  profile2->gamma = gamma0;
  profile2->l = l0;
  double shift2 = profile1->Shift(Er0,a0,l0,mu0,z1_0,z2_0);


  TFile ifile("/home/Carbon8/daq/corr_7Be.root");
  TH1I * et = (TH1I*) ifile.Get("Be6/Et_6Be_he3he3");


  sle SLE(3);
  SLE.clear();
  for (int i=1;i<200;i++)
    {
      double E = et->GetBinCenter(i);
      double P = profile1->P_l(E,a0,profile1->l,mu0,z1_0,z2_0);
      double gamma = 2.*P*profile1->gamma;
      double delta = -profile1->gamma*(profile1->Shift(E,a0,profile1->l,mu0,z1_0,z2_0)-shift1); 

      double y1 = gamma/(pow(E-profile1->Er-delta,2)+pow(gamma/2.,2));

      //cout << E << " "<< delta << endl;


      P = profile2->P_l(E,a0,profile2->l,mu0,z1_0,z2_0);
      gamma = 2.*P*profile2->gamma;
      delta = -profile2->gamma*(profile2->Shift(E,a0,profile2->l,mu0,z1_0,z2_0)-shift2); 

      double y2 = gamma/(pow(E-profile2->Er-delta,2)+pow(gamma/2.,2));



      double b = E*E*exp(-E/temp);


      double eff = .196 - E*9.194e-3;
      if (eff < 0.) eff = 0.;

      y1 *= eff;
      y2 *= eff;
      b*= eff;

      double yx = et->GetBinContent(i);

      SLE.Y[0] += yx*y1;
      SLE.Y[1] += yx*y2;
      SLE.Y[2] += yx*b;
      SLE.M[0][0]  += y1*y1;
      SLE.M[1][1]  += y2*y2;
      SLE.M[2][2] += b*b;
      SLE.M[0][1] += y1*y2;
      SLE.M[0][2] += y1*b;
      SLE.M[1][2] += y2*b;

    }

  SLE.M[2][1] = SLE.M[1][2];
  SLE.M[2][0] = SLE.M[0][2];
  SLE.M[1][0] = SLE.M[0][1];

  SLE.solve();

  cout << SLE.Y[0] << " " << SLE.Y[1] << " " << SLE.Y[2] << endl; 

  TFile fout("out.root","RECREATE");
  TCanvas canvas("out");
  TH2S frame("frame","",10,0,20,10,0,5000);
  frame.Draw();
  et->Draw("same");

  TH1I* ot = (TH1I*) et->Clone();
  TH1I* oy1 = (TH1I*) et->Clone();
  TH1I* oy2 = (TH1I*) et->Clone();
  TH1I* ob = (TH1I*) et->Clone();


  double chi = 0;
  for (int i=1;i<200;i++)
    {
      double E = et->GetBinCenter(i);
      double P = profile1->P_l(E,a0,profile1->l,mu0,z1_0,z2_0);
      double gamma = 2.*P*profile1->gamma;
      double delta = -profile1->gamma*(profile1->Shift(E,a0,profile1->l,mu0,z1_0,z2_0)-shift1); 

      double y1 = gamma/(pow(E-profile1->Er-delta,2)+pow(gamma/2.,2));


      P = profile2->P_l(E,a0,profile2->l,mu0,z1_0,z2_0);
      gamma = 2.*P*profile2->gamma;
      delta = -profile2->gamma*(profile2->Shift(E,a0,profile2->l,mu0,z1_0,z2_0)-shift2); 

      double y2 = gamma/(pow(E-profile2->Er-delta,2)+pow(gamma/2.,2));

      double b = E*E*exp(-E/temp);


      double eff = .196 - E*9.194e-3;
      if (eff < 0.) eff = 0.;

      y1 *= eff;
      y2 *= eff;
      b*= eff;

      double y = SLE.Y[0]*y1 + SLE.Y[1]*y2 + SLE.Y[2]*b;
      //cout << E << " " << y << endl;

      ot->SetBinContent(i,y);
      oy1->SetBinContent(i,SLE.Y[0]*y1);
      oy2->SetBinContent(i,SLE.Y[1]*y2);
      ob->SetBinContent(i,SLE.Y[2]*b);


      chi += pow(y-et->GetBinContent(i),2);
    }
   ot->Draw("L same");
   oy1->SetLineColor(2);
   oy1->Draw("L same");

   oy2->SetLineColor(4);
   oy2->Draw("L same");
   

   ob->SetLineColor(6);
      ob->Draw("L same");

      TLine line;
      line.SetLineStyle(2);
      line.DrawLine(profile1->Er,0,profile1->Er,5000);
      line.DrawLine(profile2->Er,0,profile2->Er,5000);

canvas.Write();
fout.Write();
  
 cout << "chisq + " << chi << endl;
  
  
}
