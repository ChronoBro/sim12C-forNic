#ifndef SEQ_PROFILE_H
#define SEQ_PROFILE_H

#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <iostream>

using namespace std;

class Profile
{
 public:
  double Et;
  Profile(double Et, double Er, int z1, int z2, int z3, double mu123, 
double ac123, int l1, double mu23, double ac23, int l2, double rwidth2_23, 
double B23);
  Profile(double Er, int z1,int z2,double mu, double ac, int l, double rwidth2,
	  double B);

  Profile(){};

  ~Profile();

  double rand(double xran);

  double Gamma(double E, double Et, double r, double l, double mu, int z1, int z2, double rwidth2);
  double Shift(double E, double r, double l, double mu, int z1, int z2);


  int N;
  double gamma;
  double Er;
  int l;
  double *profile;
  double dE;
  double P_l(double E, double r, double l, double mu, int z1, int z2);

 private:
  double Sommerfeld(double E, double mu, int z1, int z2);


  static const double hbarc;
  static const double amu;
  static const double pi;




};



#endif //SEQ_PROFILE_BE9_H
