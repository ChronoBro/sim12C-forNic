#include "plf.h"
#include <iostream>
#include "TH1S.h"
#include "TH1I.h"
#include "TH2S.h"
#include "TFile.h"
#include "frag.h"
#include "decay.h"
#include "loss.h"
#include "matrix.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

Double_t Gaus3(Double_t * x, Double_t * par)
{
  Double_t result;

  result = par[6]*TMath::Gaus(x[0],par[0],par[1]);
  result += par[7]*TMath::Gaus(x[0],par[2],par[3]);
  result += par[8]*TMath::Gaus(x[0],par[4],par[5]);

  return result;
}



int main()
{


  float const alpha = 9.23/180.*acos(-1.);

  bool einstein = 1; //switch for newtonian(0) or relativistic(1) 

  string LossHe("Helium.loss");
  string LossH("Hydrogen.loss");
  string LossLi("Lithium.loss");
  string H3thres("H3.thres");
  string He4thres("He4.thres");
  CLoss Ploss(LossLi,7.);
 

  float mass_target=12;
  int throwawayCounts=0;
  int excludedCounts=0;
  int const Nevents = 5000000;
  float thickness = 9.60;  //mg/cm2
  //thickness = 0.01; //for resolution determination
  float targetSize = 0.4;//0.4;//0.5; // spotsize at target in cm
  float spotSizeAtHira = 0.4;//0.4;//0.4;//0.5; //spotsize at detector in cm

  //adding for resolution determination
  //targetSize = 0.;
  //spotSizeAtHira = 0.;

  
  float CsiRes = .085;
  //CsiRes = 0.; //adding line for resolution determination
  //float CsiRes = .085/2.;
  float momAccept = .000; // A1900 momentum acceptance 


  //float thickness = 0.;
  //float const targetSize = .0;
  //float const spotSizeAtHira = 0.;
  //
  //float momAccept = 0.;

  //float EPA0 = 23.938;
  float EPA0 = 24.;
  //energy loss in 1/2 of target
  float epa2_mean = Ploss.getEout(EPA0*7.,thickness/2.)/7.;
  float Pbeam2_mean = sqrt(pow((epa2_mean+931.478)*7.,2)-pow(7.*931.478,2));


  //reaction center of mass velocity
  float Ecm = (epa2_mean+931.478)*7. + mass_target*931.478;
  float pcm = Pbeam2_mean;
  float vCM = pcm/Ecm*30.;

  cout << "Pbeam2_mean=" << Pbeam2_mean << endl; 
  cout << "vCM= " << vCM << endl;

  CPlf plf(thickness);
  //CFrag frag1(1.,3.,LossH,CsiRes,thickness,H3thres); //triton 
  CFrag frag1(1.,1.,LossH,CsiRes,thickness,H3thres);  //proton
  //CFrag frag1(1.,2.,LossH,CsiRes,thickness,H3thres); //going to try 6Li, I'm an idiot for not making a separate file or save for the two but its too late now
  //CFrag frag2(2.,4.,LossHe,CsiRes,thickness,He4thres); // alpha
  CFrag frag2(2.,6.,LossHe,CsiRes,thickness,He4thres); // 6He



  CDecay decay(&frag1,&frag2,einstein,vCM);


  TH1I hist_Solution("NSolution","",20,0,20);
  TH1I hist_ISolution("ISolution","",20,0,20);
  TH1F hist_cosPsi_R("cosPsi_R","",50,-1,1);
  hist_cosPsi_R.GetXaxis()->SetTitle("\\cos\\psi \\,(reconstructed)");

  TH1F hist_cosPsi_lowAngle_R("cosPsi_lowAngle_R","",50,-1,1);
  TH1F hist_cosPsi_midAngle_R("cosPsi_midAngle_R","",50,-1,1);
  TH1F hist_cosPsi_hiAngle_R("cosPsi_hiAngle_R","",50,-1,1);

  TH1I hist_cosPsi_R_peak1("cosPsi_R_peak1","",50,-1,1);
  TH1I hist_cosPsi_R_peak2("cosPsi_R_peak2","",50,-1,1);
  TH1I hist_cosPsi_R_center("cosPsi_R_center","",50,-1,1);
  TH1I hist_cosPsi_R_side("cosPsi_R_side","",50,-1,1);
  TH1I hist_cosPsi_R_away("cosPsi_R_away","",50,-1,1);
  TH1F hist_cosPsi_P("cosPsi_P","",50,-1,1);
  hist_cosPsi_P.GetXaxis()->SetTitle("\\cos\\psi \\,(primary)");


  TH1F hist_cosPsi_lowAngle_P("cosPsi_lowAngle_P","",50,-1,1);
  TH1F hist_cosPsi_midAngle_P("cosPsi_midAngle_P","",50,-1,1);
  TH1F hist_cosPsi_hiAngle_P("cosPsi_hiAngle_P","",50,-1,1);


  TH1I hist_cosPsiQ_R("cosPsiQ_R","",50,-1,1);
  TH1I hist_cosPsiQ_P("cosPsiQ_P","",50,-1,1);
  TH1I hist_cosPsiQP_R("cosPsiQP_R","",50,-1,1);
  TH1I hist_cosPsiQP_P("cosPsiQP_P","",50,-1,1);
  TH1F hist_chi_R("chi_R","",90,0,360.);
  hist_chi_R.GetXaxis()->SetTitle("\\chi \\, (reconstructed) \\, [deg]");
  TH1I hist_chi_midAngle_R("chi_midAngle_R","",90,0,360.);
  TH1I hist_chi_hiAngle_R("chi_hiAngle_R","",90,0,360.);
  TH1I hist_chi_R_peak1("chi_R_peak1","",90,0,360.);
  TH1I hist_chi_R_peak2("chi_R_peak2","",90,0,360.);
  TH1I hist_chiQ_R("chiQ_R","",90,0,360.);
  TH1I hist_chiQ_P("chiQ_P","",90,0,360.);
  TH1I hist_chiQP_R("chiQP_R","",90,0,360.);
  TH1I hist_chiQP_P("chiQP_P","",90,0,360.);
  TH1F hist_chi_P("chi_P","",90,0,360.);
  hist_chi_P.GetXaxis()->SetTitle("\\chi \\, (primary) \\, [deg]");
  TH1I hist_chi_midAngle_P("chi_midAngle_P","",90,0,360.);
  TH1I hist_chi_hiAngle_P("chi_hiAngle_P","",90,0,360.);
  TH1I hist_ExTarget("ExTarget","",200,-15,15);
  hist_ExTarget.GetXaxis()->SetTitle("\\ E^*_{Target} \\, (reconstructed)\\, [MeV]");
  hist_ExTarget.GetYaxis()->SetTitle("counts");
  TH1I hist_EPA("EPA","",100,40,70);

  TH1I hist_chi_CsI_exclude("chi_CsI_exclude","",90,0,360.);
  TH1I hist_chi_Ring_exclude("chi_Ring_exclude","",90,0,360.);
  TH1I hist_chi_beamhole("chi_beamhole","",90,0,360.);
  TH1I hist_chi_wide("chi_wide","",90,0,360.);
  TH1I hist_chi_lucky("chi_lucky","",90,0,360.);
  TH1I hist_chi_brick("chi_brick","",90,0,360.);

  TH1I cosPsi_CsI_excludeT("cosPsi_CsI_excludeT","",50,-1,1);
  TH1I cosPsi_Ring_exclude("cosPsi_Ring_exclude","",50,-1,1);
  TH1I cosPsi_beamholeT("cosPsi_beamholeT","",50,-1,1);
  TH1I cosPsi_wideT("cosPsi_wideT","",50,-1,1);
  TH1I cosPsi_luckyT("cosPsi_luckyT","",50,-1,1);
  TH1I cosPsi_brickTinAout("cosPsi_brickTinAout","",50,-1,1);
  TH1I cosPsi_brickToutAin("cosPsi_brickToutAin","",50,-1,1);
  TH1I cosPsi_brickTinAin("cosPsi_brickTinAin","",50,-1,1);
  TH1I cosPsi_brickToutAout("cosPsi_brickToutAout","",50,-1,1);
  TH1I cosPsi_CsI_exclude("cosPsi_CsI_exclude","",50,-1,1);
  TH1I cosPsi_brick("cosPsi_brick","",50,-1,1);

  TH1I cosPsi_CsI_excludeA("cosPsi_CsI_excludeA","",50,-1,1);
  // TH1I cosPsi_Ring_excludeA("cosPsi_Ring_excludeA","",50,-1,1);
  TH1I cosPsi_beamholeA("cosPsi_beamholeA","",50,-1,1);
  TH1I cosPsi_wideA("cosPsi_wideA","",50,-1,1);
  TH1I cosPsi_luckyA("cosPsi_luckyA","",50,-1,1);
  TH1F CsIresTest("CsIresTest","",500,0,100);
  TH1F CsIresTest2("CsIresTest2","",500,0,100);
  //  TH1I cosPsi_brickA("cosPsi_brickA","",50,-1,1);



  TH1I hist_chi_plus_R("chi_plus_R","",90,0,360.);
  TH1I hist_chi_plus_P("chi_plus_P","",90,0,360.);

  TH1I hist_chi_minus_R("chi_minus_R","",90,0,360.);
  TH1I hist_chi_minus_P("chi_minus_P","",90,0,360.);

  TH1I hist_Ek1("Ek1","",100,0,2.5);
  TH1I hist_Ek2("Ek2","",100,0,2.5);
  TH1I hist_cosTheta12("cosTheta12","",20,-1,1);
  TH1I hist_Ppara("Ppara","",500,1000,4500);
  TH1I hist_Px("Px","",500,-500,500);
  TH1I hist_Py("Py","",500,-500,500);

  TH2I hist_vzvx_He4("vzvx_He4","",200,-2,14,100,-7,7);
  TH2I hist_vzvx_He3("vzvx_He3","",200,-2,14,100,-7,7);
  TH2I hist_vzvx_T("vzvx_T","",200,-2,14,100,-7,7);
  TH2I hist_vxvy_He4("vxvy_He4","",100,-7,7,100,-7,7);
  TH2I hist_vxvy_He3("vxvy_He3","",100,-7,7,100,-7,7);


  TH2F cosPsi_Chi_P("cosPsi_Chi_P","",50,-1,1,90,0,360);
  cosPsi_Chi_P.GetXaxis()->SetTitle("\\cos\\psi \\,\\,(primary)");
  cosPsi_Chi_P.GetYaxis()->SetTitle("\\chi \\,\\, (primary) [deg]");
  TH2F cosPsi_Chi_R("cosPsi_Chi_R","",50,-1,1,90,0,360);
  cosPsi_Chi_R.GetXaxis()->SetTitle("\\cos\\psi \\,\\,(reconstructed)");
  cosPsi_Chi_R.GetYaxis()->SetTitle("\\chi \\,\\, (reconstructed) [deg]");

  TH2F cosPsi_Chi_Rus("cosPsi_Chi_Rus","",50,-1,1,90,0,360);
  cosPsi_Chi_R.GetXaxis()->SetTitle("\\cos\\psi \\,\\,(reconstructed)");
  cosPsi_Chi_R.GetYaxis()->SetTitle("\\chi \\,\\, (reconstructed) [deg]");
  TH2F cosPsi_Chi_S2("cosPsi_Chi_S2","",50,-1,1,90,0,360);
  cosPsi_Chi_R.GetXaxis()->SetTitle("\\cos\\psi \\,\\,(reconstructed)");
  cosPsi_Chi_R.GetYaxis()->SetTitle("\\chi \\,\\, (reconstructed) [deg]");
  TH2F cosPsi_Chi_RusS2("cosPsi_Chi_RusS2","",50,-1,1,90,0,360);
  cosPsi_Chi_R.GetXaxis()->SetTitle("\\cos\\psi \\,\\,(reconstructed)");
  cosPsi_Chi_R.GetYaxis()->SetTitle("\\chi \\,\\, (reconstructed) [deg]");

  TH2F cosPsi_Chi_beamhole("cosPsi_Chi_beamhole","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_lucky("cosPsi_Chi_lucky","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_CsI_exclude("cosPsi_Chi_CsI_exclude","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_Ring_exclude("cosPsi_Chi_Ring_exclude","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_brick("cosPsi_Chi_brick","",50,-1,1,90,0,360);


  TH2F hist_theta_check("theta_check","",300,0,15,300,0,15);
  TH2F hist_phi_check("phi_check","",90,-180,180,90,-180,180);
  TH2F hist_chi_check("chi_check","",90,0,360,90,0,360);
  TH2F hist_cosPsi_check("cosPsi_check","",50,-1,1,50,-1,1);

  TH2F cosPsi_Chi_smallAngle_P("cosPsi_Chi_smallAngle_P","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_smallAngle_R("cosPsi_Chi_smallAngle_R","",50,-1,1,90,0,360);

  TH2F cosPsi_Chi_midAngle_P("cosPsi_Chi_midAngle_P","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_midAngle_R("cosPsi_Chi_midAngle_R","",50,-1,1,90,0,360);

  TH2F cosPsi_Chi_hiAngle_P("cosPsi_Chi_hiAngle_P","",50,-1,1,90,0,360);
  TH2F cosPsi_Chi_hiAngle_R("cosPsi_Chi_hiAngle_R","",50,-1,1,90,0,360);
   

  TH1I hist_vel_P("vel_P","vel",200,0.,15.);
  TH1F hist_theta_P("theta_P","theta",125,0,40);
  hist_theta_P.GetXaxis()->SetTitle("\\theta^*_{lab} \\,\\, (primary) \\,\\, [deg]");
  TH1F hist_theta_R("theta_R","theta",125,0,40);
  hist_theta_R.GetXaxis()->SetTitle("\\theta^*_{lab} \\,\\, (reconstructed) \\,\\, [deg]");
  //TH1I hist_thetaCM_P("thetaCM_P","thetaCM",125,0,40);
  TH1I hist_phi_P("phi_P","phi",125,0,360);
  hist_phi_P.GetXaxis()->SetTitle("\\phi^*_{lab} \\.\\. (primary) \\,\\, [deg]");

  TH1I hist_vel_R("vel_R","vel",200,0.,15.);
  TH1F hist_theta_reactionCM_R("theta_reactionCM_R","theta",125,0,40);
  hist_theta_reactionCM_R.GetXaxis()->SetTitle("\\theta^*_{cm} \\,\\, (reconstructed) \\,\\, [deg]");
  TH1F hist_theta_reactionCM_R_xsection("theta_reactionCM_R_xsection","theta",125,0,40);
  hist_theta_reactionCM_R.GetXaxis()->SetTitle("\\theta^*_{cm} \\,\\, (reconstructed) \\,\\, [deg]");
  TH1F hist_theta_reactionCM_P("theta_reactionCM_P","theta",125,0,40);
  hist_theta_reactionCM_P.GetXaxis()->SetTitle("\\theta^*_{cm} \\,\\, (primary) \\,\\, [deg]");
  TH1F hist_theta_reactionCM_P_xsection("theta_reactionCM_P_xsection","theta",125,0,40);
  TH1F hist_theta_reactionCM_Detector_axis_P("theta_reactionCM_Detector_axis_P","theta",125,0,40);
  hist_theta_reactionCM_P.GetXaxis()->SetTitle("\\theta^*_{cm} \\,\\, (primary) \\,\\, [deg]");
  TH1I hist_phi_R("phi_R","phi",125,0,360);
  hist_phi_R.GetXaxis()->SetTitle("\\phi^*_{lab} \\.\\. (reconstructed) \\,\\, [deg]");
  TH1I hist_theta_emission_R("thetaEmissionR","thetaEmission",100,0,180);
  TH1F hist_thetaInput("thetaInput","thetaInput",125,0,40);
  
  TH1F hist_theta_test("theta_test","",125,0,40);
  TH1F hist_theta_test2("theta_test2","",125,0,40);


  TH1I hist_E_triton_R("E-triton","E-triton",200,0,200.);//stuff i added
  hist_E_triton_R.GetXaxis()->SetTitle("TKE\\, of \\, Triton \\,\\, (reconstructed)\\,\\ [MeV]");
 TH1I hist_E_alpha_R("E-alpha","E-alpha",200,0,200.);//stuff i added
 hist_E_alpha_R.GetXaxis()->SetTitle("TKE\\, of \\, Alpha \\,\\, (reconstructed)\\,\\ [MeV]");//stuff i added
  TH1I hist_ET_R("ET_R","ET_R",500,0,5.);
  hist_ET_R.GetXaxis()->SetTitle("E_{T} \\,\\, (reconstructed) \\,\\, [MeV]");
  TH1I hist_ET_P("ET_P","ET_P",500,0,5.);
  hist_ET_P.GetXaxis()->SetTitle("E_{T} \\,\\, (primary) \\,\\, [MeV]");
  TH1I hist_Ex_R("Ex_R","Ex",1800,-3,15.);
  TH1I hist_Ex_R_RUS("Ex_R_Rus","Ex_Rus",700,-3,10.);
  TH1I hist_Ex_R_S2("Ex_R_S2","Ex_S2",700,-3,10.);

  
  TH2S xy("xy","xy",250,-20,20,250,-13,13);
  xy.GetXaxis()->SetTitle("\\theta_{lab} \\cos(\\phi) \\,\\, [deg]");
  xy.GetYaxis()->SetTitle("\\theta_{lab} \\sin(\\phi) \\,\\, [deg]");
  TH2S xy_t("xy_t","xy_t",250,-20,20,250,-13,13);
  TH2S xy_a("xy_a","xy_a",250,-20,20,250,-13,13);
  TH2S xy_hole("xy_hole","xy_hole",250,-20,20,250,-13,13);
  TH2F cosPsi_Chi_hole("cosPsi_Chi_hole","",50,-1,1,90,0,360);
  TH2S xy_test("xy_test","xy_test",250,-20,20,250,-13,13);
  TH2F beamAtRus("beamAtRus","beamAtRus",250,-20,-20,250,-13,-13);
  TH2F beamAtTarget("beamAtTarget","beamAtTarget",250,-20,-20,250,-13,-13);
  
  float mean =0.;
  float sig = 0.;

  int RusC = 0;
  int S2C=0;
  int RusS2C=0;


  int Ndet = 0;
  int NdetClean = 0;

  //generating rotation matrix to observe effect of beam misalignment on the data

  //uses rodrigues formula for rotation
  // k is axis of rotation

  float x1;
  float y1;
 
  
  // matrix  K(3,3);
  // matrix  K2(3,3); 
  // matrix  R(3,3); 

  float x=3.5;
  float y=-1.5;
  
  float L[3] ={0,0,354.}; //center of S2
  float B[3] = {x,y,354.}; // beam vector
  
  float Lmag = 354;
  float Bmag = sqrt(x*x+y*y+354*354);
  float LdotB=L[0]*B[0]+L[1]*B[1]+L[2]*B[2];

    
  float theta_off = acos(LdotB/(Lmag*Bmag));
  cout << "theta_off = " << theta_off*180./3.14159 << endl;
  
  float cross[3] = {B[1]*L[2]-L[1]*B[2],B[2]*L[0]-L[2]*B[0],B[0]*L[1]-L[0]*B[1]};
  float cross_mag = sqrt(pow(cross[0],2.)+pow(cross[1],2.)+pow(cross[2],2.));
  float k[3] = {cross[0]/cross_mag,cross[1]/cross_mag,cross[2]/cross_mag};
  


  // K.C[0][1]=-1.*k[2];
  // K.C[0][2]=k[1];
  // K.C[1][2]=-1.*k[0];
  // K.C[1][0]=-1.*K.C[0][1];
  // K.C[2][0] = -1.*K.C[0][2];
  // K.C[2][1]=-1.*K.C[1][2];
  // K.C = K.mult(K.I);
  

  MatrixXd K(3,3);
  MatrixXd K2(3,3);
  MatrixXd R(3,3);
  MatrixXd I(3,3);
  // K.C[0][1]=-1.*k[2];
  // K.C[0][2]=k[1];
  // K.C[1][2]=-1.*k[0];
  // K.C[1][0]=-1.*K.C[0][1];
  // K.C[2][0] = -1.*K.C[0][2];
  // K.C[2][1]=-1.*K.C[1][2];

  K << 0., -k[2], k[1],
    k[2], 0, -k[0],
    -k[1], k[0], 0;

  I << 1,0,0,
    0,1,0,
    0,0,1;

  K2 = K*K;

  R = I + sin(theta_off)*K + (1-cos(theta_off))*K2;

  cout << "R = \n" <<R << endl;


  
  // K2.C = K.C;
  // K2.C = K2.mult(K2.C);
  
  // R.C = R.I;
  // K.C = K.scale(sin(theta_off));
  // K2.C = K2.scale(1-cos(theta_off));
  // R.C = R.add(K.C);
  // //cout << "here1" << endl;
  // R.C = R.add(K2.C);
  //  abort();
  int counting=0;
  int counts1 =0;
  int counts2 = 0;
  int counts =0;
  int is_triton;

  double midGate_lo = 19; //degrees of gate
  double midGate_hi = 21;

  
  TF1 * f2 = new TF1("f2",Gaus3,0,180,9);
  Double_t params[9] = {6.888,2.641,10,.866,16,4.5,23274.,10972.,10000.};
  f2->SetParameters(params);


  
  //going to sample FRESCO disitrubion to see effects of detector system 1/22/2018 (damn still working on this...)
  TFile * FRESCO = new TFile("~/Programs/fitFRESCO/fGen/test.root");
  TGraph * AngDist = (TGraph*)FRESCO->Get("expDistInel")->Clone("AngDist");
  TGraph * inelFRESCO = (TGraph*)FRESCO->Get("inelFRESCO")->Clone("inelFRESCO");
  
  int Nangle = 1801;
  double thetaFit[Nangle];
  double normalize = AngDist->Integral(0,-1);

  cout << "AngDist Integral = " << normalize << endl;

  double integrated[Nangle];
  integrated[0]=0.;
  double convoluted[Nangle];
  double integrate = 0;
  double h=0.1;
  convoluted[0] = 0.;
  double convolute = 0.;
  double sigma = 0.6;
  //generating probability distribution
  for(int q=1;q<Nangle;q++)
    {
      thetaFit[q] = q*h;
      integrate += AngDist->Eval(thetaFit[q])*h;
      integrated[q] = integrate/normalize;

    }

  for(int q=1;q<Nangle;q++)
    {
      //integrated[q] = integrated[q]/integrate;
      convolute=0.;
      for(int t=0;t<q;t++)
	{
	  convolute += 1/sqrt(2*3.1415927)*exp(-0.5*pow(thetaFit[q]-thetaFit[t],2.)/sigma)*AngDist->Eval(thetaFit[t])*h;
	}
      //cout << convolute << endl;
      convoluted[q] = convolute/sin(thetaFit[q]*3.14159/180.);
    }

  

  TGraph * probDist = new TGraph(Nangle,integrated,thetaFit);
  probDist->SetName("probDist");

  TGraph * AngDistConvoluted = new TGraph(Nangle,thetaFit,convoluted);
  AngDistConvoluted->SetName("AngDistConvoluted");

  double normalize2 = normalize/AngDistConvoluted->Integral(0,-1);

  for(int q=1;q<Nangle;q++)
    {
      convoluted[q] = convoluted[q]/normalize2;
    }
  

  for (int i=0;i<Nevents;i++)
    {
      // cout << "theta_off = " << theta_off << endl;


      if (i%3000 == 0) cout << '\xd' <<i << flush;//endl;
      if (i==3)
	{
      cout << "i=" << i << endl;
	}
      float dthick = thickness*plf.ran.Rndm();
      //float dthick = thickness/2.;


      float ET = 2.185;//0.7117;//<-for 6Li // 2.185; <-for 7Li
      ET = 1.2646; // <-7Li (IAS)
      const float c = 30.;

      hist_ET_P.Fill(ET);

      //initial energy of beam

      float Pbeam0 = sqrt(pow((EPA0+931.478)*7.,2)-pow(931.478*7.,2));

      //momentum acceptance 
      float deltaPbeam = Pbeam0*momAccept*(1.-2.*plf.ran.Rndm());
      float Pbeam1 = Pbeam0 + deltaPbeam;
      float EPA1 = sqrt(pow(Pbeam1,2)+ pow(931.478*7.,2)) - 931.478*7.;
      EPA1 /= 7.;
      

      //energy loss in target
      float EPA2 = Ploss.getEout(EPA1*7.,thickness-dthick);
      EPA2 /= 7.;
      //EPA2 = 24.; //adding this line for resolution determination
      float Pbeam2 = sqrt(pow((EPA2+931.478)*7.,2)-pow(931.478*7.,2));
      float Etot = (EPA2+931.478)*7.+931.478*mass_target;
      //float Etot = (EPA2+931.478)*6.+mass_target*931.478;
      hist_EPA.Fill(EPA2);


      
      //find com velocity
      float velocityRCOM = Pbeam2/Etot*c;
      //cout << velocityRCOM << endl;
      float gamma = 1./sqrt(1.-pow(velocityRCOM/c,2));
      //projectile and target momentum in com
      float p = mass_target*931.478*velocityRCOM/c*gamma;
      float EtotRCOM = sqrt(pow(p,2)+pow(7.*931.478,2)) + 
	sqrt(pow(p,2)+pow(mass_target*931.478,2));

      //after interaction
      //EtotRCOM -= ET + 2.468;//1.6052;//1.5866;//2.468;//1.5866;  //for alpha + t
      EtotRCOM -= ET + 9.9754;//1.6052;//1.5866;//2.468;//1.5866;  //for 6He + p
      


      //new projectile and target momentum in RCM
      float pafter = p;
      int tries = 0;
      for(;;)
	{
	  float y =  sqrt(pow(pafter,2)+pow(7.*931.478,2)) + sqrt(pow(pafter,2)+pow(mass_target*931.478,2)); //+pow(7.*931.478,2)) + 
	  
	  float dy = pafter/sqrt(pow(pafter,2)+pow(7.*931.478,2)) +pafter/sqrt(pow(pafter,2)+pow(mass_target*931.478,2)); //+pow(7.*931.478,2)) + 
	  

	float dp = (y-EtotRCOM)/dy*.99;
        pafter -= dp;
        if (fabs(dp) < .0001) break;
	tries++;
        if (tries > 50) break;
	}       

      //find angle in RCOM
      const float temp1 = .5*acos(-1.)/180.; 
      const float temp2 = .1*acos(-1.)/180.;
     float theta;


     for (;;)
       {
         if (plf.ran.Rndm() < .31)
            {
              if (plf.ran.Rndm() < .7)
	        {
                theta = plf.ran.Gaus(2.4/180.*3.14159,.4/180.*3.14159);
	        }

              else theta = plf.ran.Gaus(2.8/180.*3.14159,.6/180.*3.14159);

            }
	 else if (plf.ran.Rndm() < .03)
         theta = plf.ran.Gaus(4.79/180.*3.14159,.6/180.*3.14159);
         else if(plf.ran.Rndm() < .016)
            {
            //float theta1= plf.rotan.Gaus(0.001,.6/180.*3.14159);
            //float theta2= plf.ran.Gaus(0.001,.6/180.*3.14159);
            //theta = pow(theta1,2) + pow(theta2,2);
	      // cout << "in" << endl;
            theta = .28/180.*3.14159*sqrt(-2.*log(plf.ran.Rndm()+1.e-15));
	    //cout << " out " << theta << endl;
            }
         else if (plf.ran.Rndm() < .08) theta = plf.ran.Gaus(5.6/180.*3.14159,.55/180.*3.14159);
         else 
           {
             if(plf.ran.Rndm()>0.2)theta = sqrt(temp1)*sqrt(-log(plf.ran.Rndm()+1.e-32));
             else theta = sqrt(temp2)*sqrt(-log(plf.ran.Rndm()+1.e-32));
           }
         if (theta > 0.)break;
       }

     theta *= 9./7.*1.5;

     float theta_check;
     float rando = plf.ran.Rndm();

     
     //7Li @ 24 MeV/A scattering off 9Be @ 24 MeV/A
     for(;;)
       {
	 if(rando < 0.05 ) theta_check = plf.ran.Gaus(3,2);
	 if(rando > 0.05 && rando < 0.2) theta_check = plf.ran.Gaus(12,2.3);
	 if(rando > 0.2 && rando<0.35)theta_check =plf.ran.Gaus(7,2);
	 if(rando > 0.35 && rando <0.55) theta_check = plf.ran.Gaus(5.5,2.3);
	 if(rando > 0.55 && rando < 0.8)theta_check = plf.ran.Gaus(16,2.);
	 if(rando > 0.8) theta_check = plf.ran.Gaus(21,3.1);
	 theta=theta_check*3.1415927/180.;
	 if(theta_check>0) break;
      }

     //     theta = acos(1.-2.*plf.ran.Rndm());

     //attempting to use FRESCO disitrbution input here 1/25/2018

     // for(;;){
     // 	 theta = probDist->Eval(rando)*3.1415927/180.;
     // 	 if(theta < 25) break;
     //   }
     theta = probDist->Eval(rando)*3.1415927/180.;
     hist_thetaInput.Fill(theta*180./3.14159,1./sin(theta));
  
     
     //theta = plf.ran.Gaus(theta,sigma.*3.1415927/180.); //adding in detector resolution of angle
     // for(int q=0;q<Nangle;q++){
     //   if(rando < integrated[q] && q>0){
     // 	 theta = thetaFit[q-1]*3.14159/180.;
     // 	 break;
     // 	   }
     //   }

     
     
     //7Li @ 24 MeV/A scattering off 12C
     // for(;;)
     //   {
     // 	 if(rando < 0.15) theta_check = plf.ran.Gaus(2,1.2);
     // 	 if(rando > 0.15 && rando<0.25)theta_check =plf.ran.Gaus(5,3);
     // 	 if(rando>0.25 && rando < 0.30)theta_check = plf.ran.Gaus(16,0.6);
     // 	 if(rando > 0.30 && rando < 0.40)theta_check = plf.ran.Gaus(13,0.6);
     // 	 if(rando>0.4 && rando < 0.55)theta_check = plf.ran.Gaus(15,0.8);
     // 	 if(rando > 0.55 && rando < 0.80)theta_check = plf.ran.Gaus(8,2);
     // 	 if(rando > 0.80 && rando < 1)theta_check = plf.ran.Gaus(9,0.95);
     // 	 theta=theta_check*3.1415927/180.;
     // 	 if(theta_check >0) break;
     // }
     


     //     for(int i=0;i<
     //theta = f2->GetRandom()*3.14159/180.;
     //theta = 15*3.14159/180.;
     

     //7Li @ 24 MeV/A on 9Be transfer to 6Li & 10Be 
     /*
     for(;;)
       {
	 if(rando < 0.05 ) theta_check = plf.ran.Gaus(1.5,1.2);
	 if(rando > 0.05 && rando < 0.15) theta_check = plf.ran.Gaus(2.5,2);
	 if(rando > 0.15 && rando<0.25)theta_check =plf.ran.Gaus(3.9,2);
	 if(rando > 0.25 && rando < 0.5)theta_check = plf.ran.Gaus(5,1.3);
	 if(rando > 0.5 && rando < 0.8)theta_check = plf.ran.Gaus(8,1.5);
	 if(rando > 0.8) theta_check = plf.ran.Gaus(10,1.7);
	 theta=theta_check*3.1415927/180.;
	 if(theta_check>0) break;
       }
     */     

     //theta = acos(1.-2.*plf.ran.Rndm()); //isotropic
     //     plf.frame->phi= acos(-1.)*2.*plf.ran.Rndm();


     
     //cout << "original theta = " << theta*180./3.14159 << endl;
     // cout << "original phi = " << plf.frame->phi << endl;

     //potentially rotation of velocity here
     
     
     //see if effeciency correction is better w/o rotation
     //rotation code below
     //it doesn't belong here! understand why Dan!

     /*
     float * help = new float[3];
     help[0]=sin(theta)*cos(plf.frame->phi);
     help[1]=sin(theta)*sin(plf.frame->phi);
     help[2]=cos(theta);
     help = R.vect_mult(help);
     theta = atan2(sqrt(pow(help[0],2.)+pow(help[1],2.)),help[2]);
     plf.frame->phi = atan2(help[1],help[0]);
     delete help;
     */

     double thetaOG = theta;
     //float  * help = new float[3];

     
     // help[0]=sin(theta)*cos(plf.frame->phi);
     // help[1]=sin(theta)*sin(plf.frame->phi);
     // help[2]=cos(theta);
     // //help = R.vect_mult(help);
     //r = R*r;
     // // cout << r[0] << endl;
     // // cout << r[1] << endl;
     // // cout << r[2] << endl;
     //theta = atan2(sqrt(pow(r[0],2.)+pow(r[1],2.)),r[2]);
     //plf.frame->phi = atan2(r[1],r[0]);


     //cout << "new theta = " << theta*180./3.14159 << endl;
     //cout << "new phi = " << plf.frame->phi << endl;
     
     


   


     //  scattered projectile in RCOM
     float pparaRCOM = pafter*cos(theta);
     float pperpRCOM = pafter*sin(theta);
     float Etotafter = sqrt(pow(pafter,2)+ pow(7.*931.478,2));
     // now in lab
     float ppara = (pparaRCOM + Etotafter*velocityRCOM/c)*gamma;
     float pperp = pperpRCOM;



     //momentum in lab
     float pp = sqrt(pow(ppara,2)+pow(pperp,2));
     //total energy in lab
     float Epp = sqrt(pow(pp,2)+pow(7.*931.478,2));//+pow(6.*931.478,2));//+pow(7.*931.478,2)); //<-change this line back to 7Li to have velocities similar 7Li (reproduce original plots)
     // projectile velocity in lab
     plf.frame->velocity = pp/Epp*c; 
     //angle in lab
     plf.frame->theta = acos(ppara/pp);
    


     //find Q
     plf.Q->v[0] = pperp*cos(plf.frame->phi);
     plf.Q->v[1] = pperp*sin(plf.frame->phi);
     plf.Q->v[2] = ppara - Pbeam2;



     plf.Q->velocity = sqrt(pow(plf.Q->v[0],2) + pow(plf.Q->v[1],2)
			     + pow(plf.Q->v[2],2));
     plf.Q->theta = acos(plf.Q->v[2]/plf.Q->velocity);
     plf.Q->phi = plf.frame->phi;


						    

     //find QP
     plf.QP->v[0] = -cos(plf.Q->theta)*cos(plf.Q->phi);
     plf.QP->v[1] = -cos(plf.Q->theta)*sin(plf.Q->phi);
     plf.QP->v[2] = sin(plf.Q->theta);
     plf.QP->velocity = 1.;
     plf.QP->theta = acos(plf.QP->v[2]);
     plf.QP->phi = plf.Q->phi;




     // target velocity in the cm
     float pparaRCOM_t = -pafter*cos(theta);
     float pperpRCOM_t = -pafter*sin(theta);
     float Etotafter_t = sqrt(pow(pafter,2) + pow(mass_target*931.478,2));
     // now in lab
     float ppara_t = (pparaRCOM_t + Etotafter_t*velocityRCOM/c)*gamma;
     float pperp_t = pperpRCOM_t;
     //momentum in lab
     float pp_t = sqrt(pow(ppara_t,2)+pow(pperp_t,2));
     //total energy in lab
     float Epp_t = sqrt(pow(pp_t,2)+pow(mass_target*931.478,2));

      // beam spot at target
      float rTarget = sqrt(plf.ran.Rndm())*targetSize;
      double theta2 = 2.*plf.pi*plf.ran.Rndm();

      

      float xTarget = rTarget*cos(theta2);
      float yTarget = rTarget*sin(theta2);

      xTarget = plf.ran.Gaus(0,targetSize/2.355);
      yTarget = plf.ran.Gaus(0,targetSize/2.355);

      beamAtTarget.Fill(xTarget,yTarget);
      
      //beam spot at hira
      rTarget = sqrt(plf.ran.Rndm())*spotSizeAtHira;
      theta2 = 2.*plf.pi*plf.ran.Rndm();

      //beam divergence x = 26.7 mrad
      //beam divergence y = 5.5 mrad

      float divX = 0.0267; //radians
      float divY = 0.0055; //radians

      //adding this for resolution determination
      //divX = 0.;
      //divY= 0.;
      
      divX = plf.ran.Gaus(0,divX);
      divY = plf.ran.Gaus(0,divY);
      
      float xSpot = rTarget*cos(theta2);
      float ySpot = rTarget*sin(theta2);

      float waistToDet = 14.9-2.5;

      float majAxis = targetSize/2. + sin(divY)*waistToDet;
      float minAxis = targetSize/2. + sin(divX)*waistToDet;

      xSpot = xTarget + tan(divX)*waistToDet;
      ySpot = yTarget + tan(divY)*waistToDet;

 


     //cout << "in " << endl;     
     //cout << (EPA2+931.478)*7. + mass_target*931.478 - Epp - Epp_t - ET -1.5866 << endl;
     //cout << Epp - 7.*931.478 + ET << endl;
     //cout << Epp_t - 7.*931.478 << endl;


      
      beamAtRus.Fill(xSpot,ySpot);
      
      float deltax = xSpot-xTarget;
      float deltay = ySpot-yTarget;

      
      
      float deltar = sqrt(pow(deltax,2)+pow(deltay,2));
      float intdirTheta = atan(deltar/14.9);
      float intdirPhi = atan2(deltay,deltax);



      plf.BeamDivergence(intdirTheta,intdirPhi);
      if (thickness > 0.)  plf.MultiScat(1.-dthick/thickness);

      //need to do rotation here if I want to do it properly 1/25/2016

      Vector3d r(plf.frame->v[0],plf.frame->v[1],plf.frame->v[2]);
      r = R*r; //removing lines for resolution determination
      plf.frame->v[0]=r[0];
      plf.frame->v[1]=r[1];
      plf.frame->v[2]=r[2];
      //plf.frame->theta = atan2(sqrt(pow(r[0],2.)+pow(r[1],2.)),r[2]);
      //plf.frame->phi = atan2(r[1],r[0]);
   
      
      //plf.QP = plf.frame;

      //transform into reaction CM frame to compare with FRESCO inputs
      float velTran[3];
      double thetaTran = plf.frame->theta;
      double phiTran = plf.frame->phi;
      double velocity2 = plf.frame->velocity;
      velTran[0] = velocity2*sin(thetaTran)*cos(phiTran);
      velTran[1] = velocity2*sin(thetaTran)*sin(phiTran);
      velTran[2] = velocity2*cos(thetaTran);

      velTran[2] = (velTran[2] - vCM)/(1-velTran[2]*vCM/pow(30.,2.));

      float magV = sqrt(pow(velTran[0],2.)+pow(velTran[1],2.)+pow(velTran[2],2.));
    
      float thetaCM = atan(pcm*sin(thetaTran)/(pcm*cos(thetaTran)-7*931.478*pcm/Ecm));
      
      hist_theta_reactionCM_P.Fill(thetaCM*180./3.1415927);//(theta*180./acos(-1.));
      if(theta > 0.0084)
	hist_theta_reactionCM_P_xsection.Fill(thetaCM*180./3.1415927,1./sin(theta));

      

      hist_theta_reactionCM_Detector_axis_P.Fill(plf.frame->theta*180./3.1415927);
      


      hist_vel_P.Fill(plf.frame->velocity);
      hist_theta_P.Fill(plf.frame->theta*180./plf.pi);

      hist_phi_P.Fill(plf.frame->phi*180./plf.pi);


      //now beak up plf (I don't believe this does anything currently)
      plf.isotropic();



      
      decay.Mode(ET);
      decay.Mode(ET,plf.frame->phi);  //for some reason deltaPhi (Chi) calculated in decay.Mode was not matching up with the reconstruction which used decay.getErel(Recon), so I changed it to calculate input Chi from decay.ErelReal() and now they match. Should figure out what's wrong since my effeciencies do not match my data quantitatively (qualitatively things line up) 3/8/2016, //I believe the issue is that decay.Mode doesn't affect the variables that decay.getErel uses, this is rather confirmed by the fact that ErelP called belowed never gets used.If I can calculate the momentum of the fragments using Psi and Chi then I should be able to use decay.getErel without issues
      //bug solved, was an issue converting from degrees to radians...
      
      //decay.Mode(ET,plf.frame->phi,alpha); 


    
      
      // float * temp = R.vect_mult(plf.frame->v);
      // plf.frame->v[0]=temp[0];
      // plf.frame->v[1]=temp[1];
      // plf.frame->v[2]=temp[2];


      

      frag1.AddVelocity(plf.frame->v);
      frag2.AddVelocity(plf.frame->v);

    

      


      float ErelP = decay.getErelReal(); //removing this line to see if code stil works without it and the additions to decay.mode, it does just need to modify the decay.partCM
      //this line gets the cosPsi to check out correctly 

      //removing the block below brings back inconsistency of chi calculations
      





      //need these lines uncommented if want consisted Chi calculations with an isotropic angular distribution
      
      float try_chi = atan2(decay.partCM[0]->v[1],decay.partCM[0]->v[0]);
      if(try_chi < 0) try_chi += 2*3.1415927;
      if(try_chi > 2*3.14159) try_chi -= 2*3.1415927;
      try_chi = try_chi - plf.frame->phi;
      if(try_chi < 0) try_chi += 2*3.1415927;
      if(try_chi > 2*3.14159) try_chi -= 2*3.1415927;
      //cout << try_chi*180/3.1415927 << endl;
      decay.deltaPhi = try_chi;
      

      /* this is the code called below for reconstructing chi
      float phiDecay = atan2(decay.partCM[0]->v[1],decay.partCM[0]->v[0]);
      if (phiDecay < 0.) phiDecay += 2.*acos(-1.);
      if (phiDecay > 2.*3.14159) phiDecay -= 2.*3.14151927;
      float deltaPhi = phiDecay -  decay.plfRecon->phi;
      if(deltaPhi < 0.) deltaPhi += 2.*acos(-1.);
      if(deltaPhi > 2.*3.14159) deltaPhi -= 2.*3.1415927;
      */

      //frustrating... 8/16/2017
      // double theta_rot = 0.; //in degrees
      // theta_rot = theta_rot*3.1415927/180.;
      // double x1,y1,z1;
      // double psi = decay.theta;
      // x1 = cos(decay.deltaPhi)*sin(psi);
      // y1 = sin(decay.deltaPhi)*sin(psi);
      // z1 = cos(psi);

      // //cout << i << " " << deltaPhi << endl;
      // //cout << i << " " << psi << endl;
      
      // double x2,y2,z2;
      // //for rotation about y-axis(in reaction-plane)
      // x2 = x1*cos(theta_rot)+z1*sin(theta_rot);
      // y2 = y1;
      // z2 = -x1*sin(theta_rot)+z1*cos(theta_rot);

      // //for rotation about x-axis (out of reaction-plane)
      // x2 = x1;
      // y2 = y1*cos(theta_rot) - z1*sin(theta_rot);
      // z2 = y1*sin(theta_rot) + z1*cos(theta_rot);

      
      //decay.theta = acos(z2);



      
      hist_cosPsi_P.Fill(cos(decay.theta));


      
      cosPsi_Chi_P.Fill(cos(decay.theta),decay.deltaPhi*180./3.14159);
      if(theta*180/3.1415927 < midGate_lo)hist_cosPsi_lowAngle_P.Fill(cos(decay.theta));
      else if (theta*180/3.1415927 < midGate_hi)hist_cosPsi_midAngle_P.Fill(cos(decay.theta));
      else hist_cosPsi_hiAngle_P.Fill(cos(decay.theta));

      
      hist_chi_P.Fill(decay.deltaPhi*180./3.14159);
      if (theta > .0349)
	{
	if (cos(decay.theta) > 0.) 
              hist_chi_plus_P.Fill(decay.deltaPhi*180./3.14159);
	else hist_chi_minus_P.Fill(decay.deltaPhi*180./3.14159);
	
	if(theta*180/3.14159 <5)
	  {
	    cosPsi_Chi_smallAngle_P.Fill(cos(decay.theta),decay.deltaPhi*180./3.14159);
	  }
	if (theta*180/3.14151927 < midGate_lo && theta*180/3.14151927<midGate_hi) 
	  {
	    hist_chi_midAngle_P.Fill(decay.deltaPhi*180./3.14159);
	    cosPsi_Chi_midAngle_P.Fill(cos(decay.theta),decay.deltaPhi*180./3.14159);
	  }
	else 
	  {
	    hist_chi_hiAngle_P.Fill(decay.deltaPhi*180./3.14159);
	    cosPsi_Chi_hiAngle_P.Fill(cos(decay.theta),decay.deltaPhi*180./3.14159);
	  }
	}
      
    
      //cout << plf.frame->theta*180/3.1415927 << endl;


      // angle of Q vector
      //rotate deutro com vector to reaction plane
      float vz_He3 = frag1.real->v[2];
      float vx_He3 = frag1.real->v[0]*cos(plf.frame->phi) 
             + frag1.real->v[1]*sin(plf.frame->phi);
      float vy_He3 = -frag1.real->v[0]*sin(plf.frame->phi) 
             + frag1.real->v[1]*cos(plf.frame->phi); 


      float angle =acos(-ppara_t/pp_t);

      //rotate so that Q is the z axis
      float vvy_He3 = vy_He3;
      float vvz_He3 = vz_He3*cos(angle) + vx_He3*sin(angle);
      float vvx_He3 = vx_He3*cos(angle) - vz_He3*sin(angle);

      float chiQ  = atan2(vvy_He3,vvx_He3);
      if (chiQ < 0.) chiQ += 2.*acos(-1.);

      if (plf.frame->theta > .0349)hist_chiQ_P.Fill(chiQ*180./acos(-1.));

      angle = plf.QP->theta;
      vvy_He3 = vy_He3;
      vvz_He3 = vz_He3*cos(angle) + vx_He3*sin(angle);
      vvx_He3 = vx_He3*cos(angle) - vz_He3*sin(angle);

      float chiQP  = atan2(vvy_He3,vvx_He3);
      if (chiQP < 0.) chiQP += 2.*acos(-1.);

      if (plf.frame->theta > .0349)hist_chiQP_P.Fill(chiQP*180./acos(-1.));



      float dotQ = 0.;
      float dotQP = 0.;
      pp = 0.;
      float tt = 0.;
      for (int i=0;i<3;i++)
	{
	  dotQ += decay.partCM[0]->v[i]*plf.Q->v[i];
	  dotQP += decay.partCM[0]->v[i]*plf.QP->v[i];
          tt += pow(decay.partCM[0]->v[i],2);
	}
      dotQ /= plf.Q->velocity*sqrt(tt);
      dotQP /= plf.QP->velocity*sqrt(tt);
      hist_cosPsiQ_P.Fill(dotQ);
      hist_cosPsiQP_P.Fill(dotQP);


      //interaction in target, continue if stopped
       if (frag1.targetInteraction(dthick,thickness)) continue;
       if (frag2.targetInteraction(dthick,thickness)) continue;

      //removing above lines for resolution determination
     
      //frag1.real->theta = decay.ran.Rndm()*3.14159/4; //just for test
      // frag2.real->theta = decay.ran.Rndm()*3.14159/4;
      

      //int nhit = frag1.hit2(xTarget,yTarget) + frag2.hit2(xTarget,yTarget);


      //use_me->triton in Rus, use_2->triton in S2, use_me3->alpha in Rus, use_me4->alpha in S2

      is_triton=1;
      int nhit1 = frag1.hit6(is_triton);
      int tRus = frag1.is_hit;
      int tS2 = frag1.is_hit2;


      CsIresTest.Fill(frag1.real->energy + sqrt(frag1.real->energy)*CsiRes*frag1.ran.Gaus(0.,1.));
      CsIresTest2.Fill(frag1.ran.Gaus(frag1.real->energy,CsiRes));
      
      is_triton =0;
      int nhit2 = frag2.hit6(is_triton);
      int aRus = frag2.is_hit;
      int aS2 = frag2.is_hit2;

      



      float ErelR = decay.getErelRecon();

      //   cout << aRus << " " << aS2 << endl;  

      // need to put scenarios here and decompose "holes" in detection

      float vv = decay.partCM[0]->v[2]/decay.partCM[0]->velocity;


      float phiDecay = atan2(decay.partCM[0]->v[1],decay.partCM[0]->v[0]);
      if (phiDecay < 0.) phiDecay += 2.*acos(-1.);
      if (phiDecay > 2.*3.14159) phiDecay -= 2.*3.14151927;
      float deltaPhi = phiDecay -  decay.plfRecon->phi;
      if(deltaPhi < 0.) deltaPhi += 2.*acos(-1.);
      if(deltaPhi > 2.*3.14159) deltaPhi -= 2.*3.1415927;

      float vrel[3];
      float velocity=0;
      for (int ll=0;ll<3;ll++)
	{
	  vrel[ll] = decay.partCM[0]->v[ll]-decay.partCM[1]->v[ll];
	  velocity  += pow(vrel[ll],2.);
	}

      velocity = sqrt(velocity);
      vv = vrel[2]/velocity;
      vv = decay.partCM[0]->v[2]/decay.partCM[0]->velocity;

      
      // deltaPhi = atan2(vrel[1],vrel[2])- decay.plfRecon->phi;
      // if(deltaPhi < 0.) deltaPhi += 2.*acos(-1.);
      // if(deltaPhi > 2.*3.14159) deltaPhi -= 2.*3.1415927;


      x =frag1.recon->theta*180./plf.pi*cos(frag1.recon->phi);//remember to change these back
      y =frag1.recon->theta*180./plf.pi*sin(frag1.recon->phi);

      if((deltaPhi*180/3.1415927 > 55 && deltaPhi*180/3.1415927 <80)||(deltaPhi*180/3.1415927 > 280 && deltaPhi*180/3.1415927 <305))
	{
	  x =frag1.real->theta*180./plf.pi*cos(frag1.recon->phi);
	  y =frag1.real->theta*180./plf.pi*sin(frag1.recon->phi);
	  xy_test.Fill(x,y);

	  x =frag2.real->theta*180./plf.pi*cos(frag1.recon->phi);
	  y =frag2.real->theta*180./plf.pi*sin(frag1.recon->phi);
	  xy_test.Fill(x,y);
	  
	}

      bool exclude = decay.OnTopOf6();
      float theta_min = atan2(frag1.Ring2->rmax,frag1.Ring2->dist);
      float theta_max = atan2(frag1.Ring->rmin,frag1.Ring->dist);


      //diagnosing "holes" in the detector setup

      if((nhit1+nhit2)==2)
	{
	  //cout << "nhit = " << nhit1 + nhit2 << endl;
	  // counts++;
	  if(tRus)
	    {
	      counts++;
	      if(aRus&&(frag1.Ring->CsIstrip==frag2.Ring->CsIstrip||frag1.Ring->hitPie==frag2.Ring->hitPie)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_exclude.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(frag1.Ring->CsIstrip!=frag1.RingCsI->CsIstrip){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeT.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(aRus&&(frag2.Ring->CsIstrip!=frag2.RingCsI->CsIstrip)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(aRus&&(frag1.Ring->hitRing==frag2.Ring->hitRing)){hist_chi_Ring_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_Ring_exclude.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_Ring_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(aS2&&(frag2.Ring2->CsIstrip!=frag2.Ring2CsI->CsIstrip)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(exclude)excludedCounts++;//cout << "Missed something with a d in Rus" << endl;
	      //if(frag1.Ring->hitPie==frag2.Ring->hitPie){hist_chi_wide.Fill(deltaPhi*180/3.1415927);cosPsi_wide.Fill(vv);}
	      // if(frag2.Ring2->hitRing>=44){cosPsi_wideA.Fill(cos(decay.theta));throwawayCounts++;}
	    }
	  else if(tS2)
	    {
	      counts++;
	      if(aS2&&(frag1.Ring2->CsIstrip==frag2.Ring2->CsIstrip||frag1.Ring2->hitPie==frag2.Ring2->hitPie)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_exclude.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(frag1.Ring2->CsIstrip!=frag1.Ring2CsI->CsIstrip){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeT.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(aS2&&(frag2.Ring2->CsIstrip!=frag2.Ring2CsI->CsIstrip)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      else if(aS2&&(frag1.Ring2->hitRing==frag2.Ring2->hitRing)){hist_chi_Ring_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_Ring_exclude.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_Ring_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(aRus&&(frag2.Ring->CsIstrip!=frag2.RingCsI->CsIstrip)){hist_chi_CsI_exclude.Fill(decay.deltaPhi*180/3.1415927);cosPsi_CsI_excludeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_CsI_exclude.Fill(cos(decay.theta),decay.deltaPhi*180/3.14159);}
	      //if(frag1.Ring2->hitPie==frag2.Ring2->hitPie){hist_chi_wide.Fill(deltaPhi*180/3.1415927);cosPsi_wide.Fill(vv);}
	      //if(frag1.Ring2->hitRing>=44){cosPsi_wideT.Fill(cos(decay.theta));throwawayCounts++;}
	      //if(aS2&&frag2.Ring2->hitRing>=44){cosPsi_wideA.Fill(cos(decay.theta));throwawayCounts++;}
	      else if(exclude)excludedCounts++;//cout << "Missed something with a d in S2" << endl;
	    }
	}
      else if((nhit1+nhit2)==1)
	{
	  // counts++;
	  if(nhit1==1)
	    {
	      counts++;
	      if(frag2.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)){hist_chi_beamhole.Fill(decay.deltaPhi*180/3.1415927);cosPsi_beamholeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(frag2.real->theta > atan2(frag1.Ring->rmax,frag1.Ring2->dist)){hist_chi_wide.Fill(decay.deltaPhi*180/3.1415927);cosPsi_wideA.Fill(cos(decay.theta));throwawayCounts++;}
	      else if(frag2.real->theta > theta_min  && frag2.real->theta < theta_max && frag2.Ring2->hitRing<44){hist_chi_lucky.Fill(decay.deltaPhi*180/3.1415927);cosPsi_luckyA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_lucky.Fill(cos(decay.theta),decay.deltaPhi*180./3.1415927);}
	      else if(frag2.Ring2->hitRing>=44){cosPsi_beamholeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(frag2.Ring2->above_lineS2){cosPsi_beamholeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      // else if(tRus&&(frag1.Ring->CsIstrip!=frag1.RingCsI->CsIstrip)) cosPsi_CsI_excludeT.Fill(cos(decay.theta));
	      // else if(tS2&&(frag1.Ring2->CsIstrip!=frag1.Ring2CsI->CsIstrip)) cosPsi_CsI_excludeT.Fill(cos(decay.theta));
	      else {excludedCounts++;}  // cout << "theta_min = " << theta_min*180/3.14159 << " theta_max = " << theta_max*180/3.14159 <<endl;excludedCounts++;cout << frag1.real->theta*180/3.14159 <<  " " << frag2.real->theta*180/3.14159 << endl << frag1.Ring2->hitRing << " " << frag2.Ring2->hitRing <<  endl;cout << frag1.Ring->hitRing << " " << frag2.Ring->hitRing << endl;}
	    }
	  else if (nhit2==1)
	    {
	      counts++;
	      if(frag1.real->theta > atan2(frag1.Ring->rmax,frag1.Ring2->dist)){hist_chi_wide.Fill(decay.deltaPhi*180/3.1415927);cosPsi_wideT.Fill(cos(decay.theta));throwawayCounts++;}
	      else if(frag1.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)){hist_chi_beamhole.Fill(decay.deltaPhi*180/3.1415927);cosPsi_beamholeT.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(frag1.real->theta > theta_min && frag1.real->theta < theta_max && frag1.Ring2->hitRing <44){hist_chi_lucky.Fill(decay.deltaPhi*180/3.1415927);cosPsi_luckyT.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_lucky.Fill(cos(decay.theta),decay.deltaPhi*180./3.1415927);}
	      else if(frag1.Ring2->hitRing>=44){cosPsi_beamholeT.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      else if(frag1.Ring2->above_lineS2){cosPsi_beamholeA.Fill(cos(decay.theta));throwawayCounts++;cosPsi_Chi_beamhole.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);}
	      // else if(aRus&&(frag2.Ring->CsIstrip!=frag2.RingCsI->CsIstrip)) cosPsi_CsI_excludeT.Fill(cos(decay.theta));
	      //else if(aS2&&(frag2.Ring2->CsIstrip!=frag2.Ring2CsI->CsIstrip)) cosPsi_CsI_excludeT.Fill(cos(decay.theta));
	      else {excludedCounts++;} //cout << "theta_min = " << theta_min*180/3.14159 << " theta_max = " << theta_max*180/3.14159 <<endl;excludedCounts++;cout << frag1.real->theta*180/3.14159 <<  " " << frag2.real->theta*180/3.14159 << endl << frag1.Ring2->hitRing << " " << frag2.Ring2->hitRing << endl;cout << frag1.Ring->hitRing << " " << frag2.Ring->hitRing << endl;}
	    }
	}
      else if((nhit1+nhit2)==0)
	{
	  hist_chi_brick.Fill(deltaPhi*180/3.1415927);
	  cosPsi_Chi_brick.Fill(cos(decay.theta),decay.deltaPhi*180/3.1415927);
	  if(frag1.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)&&frag2.real->theta > atan2(frag1.Ring->rmax,frag1.Ring2->dist)){cosPsi_brickTinAout.Fill(cos(decay.theta));counts++;throwawayCounts++;}
	  else if(frag2.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)&&frag1.real->theta > atan2(frag1.Ring->rmax,frag1.Ring2->dist)){cosPsi_brickToutAin.Fill(cos(decay.theta));counts++;throwawayCounts++;}
	  else if(frag1.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)&&frag2.real->theta < atan2(frag1.Ring2->rmin,frag1.Ring2->dist)){cosPsi_brickTinAin.Fill(cos(decay.theta));counts++;throwawayCounts++;}
	  else if(frag1.real->theta  > atan2(frag1.Ring->rmax,frag1.Ring2->dist)&&frag2.real->theta > atan2(frag1.Ring->rmax,frag1.Ring2->dist)){cosPsi_brickToutAout.Fill(cos(decay.theta));counts++;throwawayCounts++;}
	  else {cosPsi_brick.Fill(cos(decay.theta));counts++;throwawayCounts++;}


	}

      
      int nhit = nhit1+ nhit2;
      if (nhit != 2) continue;
      // if (nhit ==2 && (frag1.Ring2->hitRing>=44||frag2.Ring2->hitRing>=44))continue;
      // float sup = plf.ran.Rndm();
      // if(deltaPhi*180/3.1415927 > 170 && deltaPhi*180/3.1415927 < 190 && sup <.67)continue;



      /*
      // reactions in CsI
      float elab = frag1.real->energy/2.;
      float prob = 1. -(-3.+.125*elab/3.)/200.;
      if (decay.ran.Rndm() > prob) continue;



      elab = frag2.real->energy/4.;
      prob = 1. -(-4.33+.175*elab)/200.;
      if (decay.ran.Rndm() > prob) continue;
      */

      
      Ndet++;
      if (exclude)
	{
	  continue;
	}    
      // if (frag1.recon->getVelocity()>7.07) continue; // for 1cm thick CsI
      //if (frag1.recon->getVelocity()>8.56) continue; //velocity (cm/ns) where the triton will punch through 2cm thick CsI



      //determining coincidences in RusS2;



      NdetClean++;
      
      if(tRus) //use_me->triton in Rus, use_2->triton in S2, aRus->alpha in Rus, aS2->alpha in S2
	{
	  if(tS2){counting++;}//cout << "wtf mate" << endl;}
	  if(aRus) RusC++;
	  if(aS2) RusS2C++;
	}
      else if(tS2)
	{
	  if(aRus) RusS2C++;
	  if(aS2) S2C++;
	}

      hist_theta_check.Fill(frag2.real->theta*180./3.14159,frag2.recon->theta*180./3.14159);
      hist_phi_check.Fill(frag2.real->phi*180./3.14159,frag2.recon->phi*180./3.14159);
      hist_theta_check.Fill(frag1.real->theta*180./3.14159,frag1.recon->theta*180./3.14159);
      hist_phi_check.Fill(frag1.real->phi*180./3.14159,frag1.recon->phi*180./3.14159);
      

      //correct for energy loss in half of target and get velocity
      //taking out two frag lines below for resolution determination
      
      frag1.Egain(thickness/2.);
      frag2.Egain(thickness/2.);
      

      double thetaR = decay.plfRecon->theta;
      double phiR = decay.plfRecon->phi;
      Vector3d r2(sin(thetaR)*cos(phiR),sin(thetaR)*cos(phiR),cos(thetaR));
     //  r2[0]=sin(thetaR)*cos(plf.frame->phi);
     //  r2[1]=sin(thetaR)*sin(plf.frame->phi);
     //  r2[2]=cos(thetaR);
     //  //r2 = R.vect_mult(r2);
     // //cout << r[0] << endl;
     //  theta = atan2(sqrt(pow(r2[0],2.)+pow(r2[1],2.)),r2[2]);
     //  plf.frame->phi = atan2(r2[1],r2[0]);
      //delete r2;     

     //  R = I - sin(theta_off)*K+(1-cos(theta_off))*K2;

     //  r2 = R*r2;
      
      
     // B[0]=-B[0];
     // B[1]=-B[1];
     
     // double dot2 =r2[0]*B[0]+r2[1]*B[1]+r2[2]*B[2];
     // dot2 = dot2/sqrt((r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2])*(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));

     //thetaR = atan2(sqrt(pow(r[0],2.)+pow(r[1],2.)),r[2]);
     //phiR = atan2(r[1],r[0]);

      
     //thetaR = acos(dot);

      

      
      //float ErelR = decay.getErelRecon();
      hist_theta_reactionCM_R.Fill(decay.theta_reactionCM*180./acos(-1.));
      //hist_theta_reactionCM_R.Fill(thetaR*180./acos(-1.));

      if(theta > 0.0084)
	hist_theta_reactionCM_R_xsection.Fill(decay.theta_reactionCM*180./3.1415927,1./sin(decay.theta_reactionCM));

      hist_theta_reactionCM_R_xsection.SetBinContent(1,0);
      hist_theta_reactionCM_R_xsection.SetBinContent(125,0);
      //cout << decay.theta_reactionCM*180./acos(-1.) << endl;

      float ErelT = frag1.recon->energy; //stuff i added
      // cout << ErelT << endl;
      //reconstruct target momentum
      float ptr[3];
      ptr[0] = -decay.plfRecon->pc[0];
      ptr[1] = -decay.plfRecon->pc[1];
      ptr[2] = Pbeam2_mean -decay.plfRecon->pc[2];      


      float pptr = sqrt(pow(ptr[0],2)+pow(ptr[1],2)+pow(ptr[2],2));
      float etr = sqrt(pow(pptr,2)+pow(mass_target*931.478,2));
      float ekr = etr - mass_target*931.478;
     


      //transferred momentum Q is negative of ptr
      float dot = 0.;
      for (int i=0;i<3;i++)
	{
          dot += -ptr[i]*decay.partCM[0]->v[i];
	}
      dot = dot/pptr/decay.partCM[0]->velocity;
      hist_cosPsiQ_R.Fill(dot);
      


      float Ex_target = epa2_mean*7. - ekr - frag1.recon->energy - 
	frag2.recon->energy - 9.9754;//2.468;//1.4743; 
      hist_E_triton_R.Fill(ErelT);// stuff i added
      ErelT= frag2.recon->energy;
      hist_E_alpha_R.Fill(ErelT);
      hist_ExTarget.Fill(Ex_target);

      hist_ET_R.Fill(ErelR);//I changed R->T
      //hist_Ex_R.Fill(ErelR + 2.468);//1.4743);
      hist_Ex_R.Fill(ErelR + 9.9754);//1.4743);

      if(tRus && aRus){
	hist_Ex_R_RUS.Fill(ErelR+2.468);
      }

      if(tS2 && aS2){
      hist_Ex_R_S2.Fill(ErelR+2.468);
      }
      
      hist_Ppara.Fill(decay.plfRecon->pc[2]);
      hist_Px.Fill(decay.plfRecon->pc[0]);
      hist_Py.Fill(decay.plfRecon->pc[1]);
      hist_theta_R.Fill(decay.plfRecon->theta*180./3.14159);
      hist_phi_R.Fill(decay.plfRecon->phi*180./3.14159);
      hist_vel_R.Fill(decay.plfRecon->velocity);

      //will need to add rotation here for effeciency correction of rotated quantization axis: 8/16/2017
      //theta_rot = 0.; //in degrees
      //theta_rot = theta_rot*3.1415927/180.;
      //double x1,y1,z1;
      // psi = acos(vv);
      // x1 = cos(deltaPhi)*sin(psi);
      // y1 = sin(deltaPhi)*sin(psi);
      // z1 = cos(psi);

      // //cout << i << " " << deltaPhi << endl;
      // //cout << i << " " << psi << endl;
      
      // //double x2,y2,z2;

      // //for rotation about y-axis (in reaction-plane)
      // x2 = x1*cos(theta_rot)+z1*sin(theta_rot);
      // y2 = y1;
      // z2 = -x1*sin(theta_rot)+z1*cos(theta_rot);

      // //for rotation about x-axis (out of reaction-plane)

      // x2 = x1;
      // y2 = y1*cos(theta_rot) - z1*sin(theta_rot);
      // z2 = y1*sin(theta_rot) + z1*cos(theta_rot);
      
      //vv = z2;
     

     
      hist_cosPsi_R.Fill(vv);

      //need to fill histograms with each component of the detector array

 

     if (decay.theta_reactionCM*180/3.1415927 < midGate_lo)hist_cosPsi_lowAngle_R.Fill(vv);
     else if (decay.theta_reactionCM*180/3.1415927 < midGate_hi)hist_cosPsi_midAngle_R.Fill(vv);
     else hist_cosPsi_hiAngle_R.Fill(vv);

    

     cosPsi_Chi_R.Fill(vv,deltaPhi*180./acos(-1.));

     
     hist_chi_check.Fill(decay.deltaPhi*180/3.14159,deltaPhi*180/3.14159);
     hist_cosPsi_check.Fill(cos(decay.theta),vv);
     
     if(tRus)
       {
	 if(aRus) 
	   {
	     if(tS2||aS2) cout << "bughh" << endl;
	     cosPsi_Chi_Rus.Fill(vv,deltaPhi*180./acos(-1.));
	     /*
	     if(deltaPhi*180./acos(-1.)>150 && deltaPhi*180./acos(-1.) <210)
	       {
		 cout << " event #" << i << endl;
		 cout << "frag1.real->theta = " << frag1.real->theta << endl;
		 cout << "frag1.recon->theta = " << frag1.recon->theta << endl;
		 cout << "frag2.recon->theta = " << frag2.real->theta << endl;
		 cout << "frag2.real->theta = " << frag2.recon->theta << endl;
		 cout << "frag1.real->phi = " << frag1.real->phi << endl;
		 cout << "frag1.recon->phi = " << frag1.recon->phi << endl;
		 cout << "frag2.real->phi = " << frag2.real->phi << endl;
		 cout << "frag2.recon->phi = " << frag2.recon->phi << endl;
		 cout << "7Li actual theta = " << plf.frame->theta*180/3.1415927 << endl;
		 cout << "7Li reconstructed theta = " << decay.plfRecon->theta*180./3.14159 << endl;
		 cout << "7Li actual phi = " << plf.frame->phi*180/3.14159 << endl;
		 cout << "7Li reconstructed phi = " << decay.plfRecon->phi*180/3.14159<< endl;
		 cout << "Real cosPsi = " << cos(decay.theta) << endl;
		 cout << "Reconstructed cosPsi = " << vv << endl;
		 cout << "Real Chi = " << decay.deltaPhi*180/3.14159 << endl;
		 cout << "Reconstructed Chi = " << deltaPhi*180/3.14159 << endl;
		 hist_theta_test.Fill(decay.plfRecon->theta*180/acos(-1.));
		 hist_theta_test2.Fill(plf.frame->theta*180/3.14159);
		 if(decay.deltaPhi < deltaPhi+0.5 || decay.deltaPhi > deltaPhi-0.5 )
		   {
		     cout << "BAD!" << endl;
		   }
	       }
	     */

	   }
	 if(aS2) cosPsi_Chi_RusS2.Fill(vv,deltaPhi*180./acos(-1.));   
       }
     else
       {
	 if(aRus)cosPsi_Chi_RusS2.Fill(vv,deltaPhi*180./acos(-1.));
	 if(aS2)cosPsi_Chi_S2.Fill(vv,deltaPhi*180./acos(-1.));
	 
       }
     
     
     
     if (decay.plfRecon->theta > .0349)
       {
         hist_chi_R.Fill(deltaPhi*180./acos(-1.));
         if(vv > 0.) hist_chi_plus_R.Fill(deltaPhi*180./acos(-1.));
	 else hist_chi_minus_R.Fill(deltaPhi*180./acos(-1.));
	 
	 if(decay.theta_reactionCM*180/3.1415927 <midGate_lo)
	   {
	     cosPsi_Chi_smallAngle_R.Fill(vv,deltaPhi*180./3.1415927);
	   }
         else if (decay.theta_reactionCM*180/3.1415927 < midGate_hi)
	   {
	     hist_chi_midAngle_R.Fill(deltaPhi*180./acos(-1.));
	     cosPsi_Chi_midAngle_R.Fill(vv,deltaPhi*180./acos(-1.));
	   }
	 else 
	   {
	     hist_chi_hiAngle_R.Fill(deltaPhi*180./acos(-1.));
	     cosPsi_Chi_hiAngle_R.Fill(vv,deltaPhi*180./acos(-1.));
	   }
	 
         float cosDP = cos(deltaPhi);
         if (cosDP < -0.5)      hist_cosPsi_R_center.Fill(vv);
         else if (cosDP < 0.5)      hist_cosPsi_R_side.Fill(vv);
         else hist_cosPsi_R_away.Fill(vv);
       }
     
     
     
     float x_t,y_t,x_a,y_a;
     
       x =frag1.recon->theta*180./plf.pi*cos(frag1.recon->phi);//remember to change these back
       y =frag1.recon->theta*180./plf.pi*sin(frag1.recon->phi);
       xy.Fill(x,y);
       x_t = x;
       y_t = y;
       xy_t.Fill(x_t,y_t);
       
       
       x = frag2.recon->theta*180./plf.pi*cos(frag2.recon->phi);
       y = frag2.recon->theta*180./plf.pi*sin(frag2.recon->phi);
       xy.Fill(x,y);
       x_a = x;
       y_a = y;
       xy_a.Fill(x_a,y_a);
     
     
     
     if (decay.plfRecon->theta*180./plf.pi > 2. && decay.plfRecon->theta*180./plf.pi < 6.)
       {
	 hist_cosPsi_R_peak1.Fill(vv);
         hist_chi_R_peak1.Fill(deltaPhi*180./acos(-1.));
       }
     else if (decay.plfRecon->theta*180./plf.pi > 6. && decay.plfRecon->theta*180./plf.pi < 10.)
       {
	 hist_cosPsi_R_peak2.Fill(vv);
         hist_chi_R_peak2.Fill(deltaPhi*180./acos(-1.));
       }
     
     if((deltaPhi*180./acos(-1.)>150&&deltaPhi*180./acos(-1.)<210) || (deltaPhi*180./acos(-1.)>150&&deltaPhi*180./acos(-1.)<210))
       {
	 if(vv>-0.6 && vv < 0.7)
	   {
	     if(tRus&&aRus)
	     {
	       //xy_hole.Fill(x_t,y_t);
	     xy_hole.Fill(x_a,y_a);
		 cosPsi_Chi_hole.Fill(vv,deltaPhi*180/3.1415927);
		 // if(tS2)cout << "It says the triton also hit the S2?!" << endl;
		 // if(aS2)cout << "It says the alpha also hit the S2?!" << endl;


		 /*
		 if(fabs(decay.deltaPhi-deltaPhi)>0.3 || fabs(cos(decay.theta)-vv) > 0.3)
		   {
		 cout << "Real decay.phi = " << decay.deltaPhi*180./3.1415927 << endl;
		 cout << "Recon decay.phi = " << deltaPhi*180./3.1415927 << endl << endl;

		 cout << "Real cosPsi = " << cos(decay.theta) << endl;
		 cout << "Recon cosPsi = " << vv << endl << endl;

		 cout << "frag1 real thetaCM = " << frag1.real->theta*180./3.1415927 << endl;
		 cout << "frag1 recon thetaCM = " << frag1.recon->theta*180./3.1415927 << endl << endl;
		 
		 cout << "frag1 real phi = " << frag1.real->phi*180./3.1415927 << endl;
		 cout << "frag1 recon phi = " << frag1.recon->phi*180./3.1415927 << endl << endl;

		 cout << "frag1 real v = " << frag1.real->velocity << endl;
		 cout << "frag1 recon v = " << frag1.recon->velocity << endl << endl;

		 cout << "frag2 real theta = " << frag2.real->theta*180./3.1415927 << endl;
		 cout << "frag2 recon theta = " << frag2.recon->theta*180./3.1415927 << endl << endl;
		 
		 cout << "frag2 real phi = " << frag2.real->phi*180./3.1415927 << endl;
		 cout << "frag2 recon phi = " << frag2.recon->phi*180./3.1415927 << endl << endl;
		 
		 cout << "frag2 real v = " << frag2.real->velocity << endl;
		 cout << "frag2 recon v = " << frag2.recon->velocity << endl << endl;

		   }
		 */
	     }
	   }
       }
     
     
     
     
     





      vz_He3 = frag1.recon->v[2];
      vx_He3 = frag1.recon->v[0]*cos(decay.plfRecon->phi) + frag1.recon->v[1]*sin(decay.plfRecon->phi);
      vy_He3 = -frag1.recon->v[0]*sin(decay.plfRecon->phi) + frag1.recon->v[1]*cos(decay.plfRecon->phi); 
      
      float vz_He4 = frag2.recon->v[2];
      float vx_He4 = frag2.recon->v[0]*cos(decay.plfRecon->phi) + frag2.recon->v[1]*sin(decay.plfRecon->phi);
      float vy_He4 = -frag2.recon->v[0]*sin(decay.plfRecon->phi) + frag2.recon->v[1]*cos(decay.plfRecon->phi); 

      hist_vzvx_He4.Fill(vz_He4,vx_He4);
      hist_vzvx_He3.Fill(vz_He3,vx_He3);
      hist_vxvy_He4.Fill(vx_He4,vy_He4);
      hist_vxvy_He3.Fill(vx_He3,vy_He3);

      // angle of Q vector
      //rotate deutro com vector to reaction plane
      vz_He3 = decay.partCM[0]->v[2];
      vx_He3 = decay.partCM[0]->v[0]*cos(decay.plfRecon->phi) 
             + decay.partCM[0]->v[1]*sin(decay.plfRecon->phi);
      vy_He3 = -decay.partCM[0]->v[0]*sin(decay.plfRecon->phi) 
             + decay.partCM[0]->v[1]*cos(decay.plfRecon->phi); 


      angle = acos(-ptr[2]/pptr);

      //rotate so that Q is the z axis
      vvy_He3 = vy_He3;
      vvz_He3 = vz_He3*cos(angle) + vx_He3*sin(angle);
      vvx_He3 = vx_He3*cos(angle) - vz_He3*sin(angle);

      chiQ  = atan2(vvy_He3,vvx_He3);
      if (chiQ < 0.) chiQ += 2.*acos(-1.);

      if (decay.plfRecon->theta > .0349)hist_chiQ_R.Fill(chiQ*180./acos(-1.));


      float QP[3];  //perpendiculat to Q
      float angleQP = angle +0.*.0175;
      QP[0] = -cos(angle)*cos(decay.plfRecon->phi);
      QP[1] = -cos(angle)*sin(decay.plfRecon->phi);
      QP[2] = sin(angle);

      dotQP = 0.;
      for (int i=0;i<3;i++)
	  {
	    dotQP += QP[i]*decay.partCM[0]->v[i];
	  } 
     float cosPsiQP = dotQP/sqrt(tt);  
     hist_cosPsiQP_R.Fill(cosPsiQP);




      //rotate so that Q is the z axis
     angle = acos(QP[2]);
      vvy_He3 = vy_He3;
      vvz_He3 = vz_He3*cos(angle) + vx_He3*sin(angle);
      vvx_He3 = vx_He3*cos(angle) - vz_He3*sin(angle);

      chiQP  = atan2(vvy_He3,vvx_He3);
      if (chiQP < 0.) chiQP += 2.*acos(-1.);

      if (decay.plfRecon->theta > .0349)hist_chiQP_R.Fill(chiQP*180./acos(-1.));



      //cout << dot << " " << vvz_He3/sqrt(pow(vvx_He3,2)+pow(vvy_He3,2)+pow(vvz_He3,2)) << endl;

    }

  //cout << "detection efficiency = " << (float)Ndet/(float)Nevents << endl;
  cout << "clean detection efficiency (separate stips and Csi) = " 
       << (float)NdetClean/(float)Nevents << endl;
  cout << "Fraction of coincidences in both  = " << (float)RusS2C/(float)NdetClean << endl;
  cout << "Fraction in S2 = " <<  (float)S2C/(float)NdetClean << endl;
  cout << "Fraction in Rus = " <<  (float)RusC/(float)NdetClean << endl;
  cout << "Total Check2(Rus+S2+RusS2)  " << RusC+S2C+RusS2C<< " NdetClean = " << NdetClean << endl;
  cout << "Total Check3(counts) " << counts << endl;
  cout << "Weirdness (counting)" << counting << endl;
  cout << "Throwaway Counts " << throwawayCounts << endl;
  cout << "Throwaway Check: Missing " << Nevents - throwawayCounts - NdetClean << endl;
  cout << "Excluded Counts (Should equal Line above) " << excludedCounts << endl;
  mean = mean/(float)Ndet;
  sig = sig/(float)(Ndet-1) - (float)(Ndet)/(float)(Ndet-1)*pow(mean,2);
  sig = sqrt(sig);
  cout << "mean Erel= " << mean << " mean Ex " << mean 
       << " sig= " << sig << " FWHM= " 
  << sig*2.35 << endl;
  


  TFile f ("sim.root","RECREATE");

  hist_vel_P.Write();
  hist_vel_R.Write(); 
  hist_theta_P.Write();
  hist_theta_R.Write();
  hist_phi_P.Write();
  hist_phi_R.Write();
  hist_ET_R.Write();
  hist_ET_P.Write();


  hist_theta_reactionCM_P.Write();
  hist_theta_reactionCM_R.Write();
  hist_theta_reactionCM_Detector_axis_P.Write();

  
  hist_theta_reactionCM_P_xsection.Write();
  hist_theta_reactionCM_R_xsection.Write();
  hist_thetaInput.Write();

  hist_cosPsi_P.Write();
  hist_cosPsi_R.Write();

  hist_chi_P.Write();
  hist_chi_R.Write();


  cosPsi_Chi_P.Write();
  cosPsi_Chi_R.Write();

  cosPsi_Chi_Rus.Write();
  cosPsi_Chi_S2.Write();
  cosPsi_Chi_RusS2.Write();
  cosPsi_Chi_hole.Write();

  hist_chi_CsI_exclude.Write();
  hist_chi_Ring_exclude.Write();
  hist_chi_beamhole.Write();
  hist_chi_wide.Write();
  hist_chi_lucky.Write();
  hist_chi_brick.Write();

  cosPsi_Chi_CsI_exclude.Write();
  cosPsi_Chi_Ring_exclude.Write();
  cosPsi_Chi_beamhole.Write();
  cosPsi_Chi_lucky.Write();
  cosPsi_Chi_brick.Write();


  cosPsi_CsI_excludeT.Write();
  cosPsi_Ring_exclude.Write();
  cosPsi_beamholeT.Write();
  cosPsi_wideT.Write();
  cosPsi_luckyT.Write();
  //cosPsi_brickT.Write();

  cosPsi_CsI_excludeA.Write();
  //cosPsi_Ring_exclude.Write();
  cosPsi_beamholeA.Write();
  cosPsi_wideA.Write();
  cosPsi_luckyA.Write();
  cosPsi_brick.Write();
  cosPsi_brickTinAout.Write();
  cosPsi_brickToutAin.Write();
  cosPsi_brickTinAin.Write();
  cosPsi_brickToutAout.Write();


  hist_Ex_R.Write();
  hist_Ex_R_RUS.Write();
  hist_Ex_R_S2.Write();
  hist_ExTarget.Write(); 
  hist_Ppara.Write();
  hist_Px.Write();
  hist_Py.Write();
  hist_E_triton_R.Write();
  hist_E_alpha_R.Write();

  hist_theta_check.Write();
  hist_phi_check.Write();
  hist_chi_check.Write();
  hist_cosPsi_check.Write();

  hist_theta_test.Write();
  hist_theta_test2.Write();

  hist_cosPsi_lowAngle_P.Write();
  hist_cosPsi_midAngle_P.Write();
  hist_cosPsi_hiAngle_P.Write();

  hist_cosPsi_lowAngle_R.Write();
  hist_cosPsi_midAngle_R.Write();
  hist_cosPsi_hiAngle_R.Write();
  hist_cosPsi_R_peak1.Write();
  hist_cosPsi_R_peak2.Write();
  hist_cosPsi_R_center.Write();
  hist_cosPsi_R_side.Write();
  hist_cosPsi_R_away.Write();

  hist_chi_midAngle_R.Write();
  hist_chi_hiAngle_R.Write();

  hist_chi_R_peak1.Write();
  hist_chi_R_peak2.Write();
  xy.Write();
  xy_t.Write();
  xy_a.Write();
  xy_hole.Write();
  xy_test.Write();
  beamAtRus.Write();
  beamAtTarget.Write();
  
  hist_vzvx_He4.Write();
  hist_vzvx_He3.Write();
  hist_vzvx_T.Write();
  hist_vxvy_He4.Write();
  hist_vxvy_He3.Write();

  hist_chi_midAngle_P.Write();
  hist_chi_hiAngle_P.Write();

  hist_EPA.Write();
  hist_cosPsiQ_R.Write();
  hist_cosPsiQ_P.Write();
  hist_cosPsiQP_R.Write();
  hist_cosPsiQP_P.Write();
  hist_chiQ_R.Write();
  hist_chiQ_P.Write();
  hist_chiQP_R.Write();
  hist_chiQP_P.Write();

  cosPsi_Chi_smallAngle_P.Write();
  cosPsi_Chi_smallAngle_R.Write();
  cosPsi_Chi_midAngle_P.Write();
  cosPsi_Chi_midAngle_R.Write();
  cosPsi_Chi_hiAngle_P.Write();
  cosPsi_Chi_hiAngle_R.Write();
  hist_chi_plus_R.Write();
  hist_chi_plus_P.Write();
  hist_chi_minus_R.Write();
  hist_chi_minus_P.Write();
  AngDistConvoluted->Write();
  probDist->Write();
  AngDist->Write();
  inelFRESCO->Write();
  CsIresTest.Write();
  CsIresTest2.Write();
  
  f.Write();


  cout << "here2!" << endl;


}
