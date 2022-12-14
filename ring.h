#ifndef _ring
#define _ring
#include "random.h"
#include "matrix.h"
#include <iostream>
using namespace std;

class ring
{
 public:
  static CRandom ran;
  ring(float dist,float rmin, float rmax, int Nring, int Npi);
  int hit(float theta, float phi);

  int hitRing; // ring which was hit
  int hitPie;  // pie sector which was hit

  float thetaHit; //zenith angle in radians assigned to hit
  float phiHit; // azimuth angle in radians assigned to hit

  float rmin;  //minium radius of ring counter
  float rmax;  //maximum radius of ring counter
  float dist;  //distance of ring counter from target
  int Nring;   // number of annullar ring sectors
  int Npie;    //number of pie sectors
  float theta_min;
  float theta_max;
  float dtheta;
  float dphi;
  float CsIstrip;
  float dr;
  int useful;
  int please;
  float * vect;
  int S2=0;
  int above_lineS2;

};
#endif
