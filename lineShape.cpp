#include "Profile.h"
#include "coul.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
int main()
{

  ofstream file("line.dat");
  coul Coul;
  Profile * profile;
  double r0 = 1.45;

  //5Li decay parameters
  double Er2;  //Q1p
  double gamma2;
  double B2;
  int l2=1;
  int z1_2 = 1;
  int z2_2 = 2;

  double a2 = 5.5;//r0*(pow(3.,1./3.)+pow(4.,1./3.));
  double mu2 = 1.*4./5.;
  gamma2 = .952;
  Er2 = 1.861;
  B2 = profile->Shift(Er2,a2,l2,mu2,z1_2,z2_2);
  cout << "B2 = " << B2 << endl;

  //find normalization
  double dE = .1;
  double sumY = 0.;
  double Y2[1000];
  double UU[1000];
  for (int i= 0;i<200;i++)
    {
      double U = ((double)i+.5)*dE;
      UU[i] = U;
      double P2 = profile->P_l(U,a2,l2,mu2,z1_2,z2_2);
      double S2 = -gamma2*(profile->Shift(U,a2,l2,mu2,z1_2,z2_2)-B2);
      double G2 = 2.*gamma2*P2;
      Y2[i] = G2/(pow(U-Er2-S2,2)+pow(G2/2.,2));
      file << U -1.9642<< " " << Y2[i] << endl;
    }

}
