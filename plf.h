#include "random.h"
#include <iostream>
#include "frame.h"
#include "mScat.h"

using namespace std;

/**
 * !\brief parent fragment stuff
 *
 * deals with the velocity and angular distribution of the parent
 * fragment before its decays
 */

class CPlf
{
 public:
  float vel;  //!< velocity of parent fragment in cm/ns
  float theta; //!< theta angle of parent fragment in radians
  float phi; //!< phi angle
  float v[3]; //!< velcoity vector of parent fragment in cm/ns

  float theta_breakup;
  float phi_breakup;
  float unit_breakup[3];

  CPlf(float);
  ~CPlf();
  void GetPlf(float);
  void isotropic();
  void inplane();
  float temp1,temp2;
  static float const pi;
  static CRandom ran;
  CFrame *frame;
  CFrame *Q;
  CFrame *QP;
  CMultScat *multScat;
  void MultiScat(float);
  void BeamDivergence(float,float);

  int npara;
  double Ppara[200];
  double ypara[200];
  float mean;


  double angArray[200];


  int N;
  double xx[200];
  double yy[200];

};
