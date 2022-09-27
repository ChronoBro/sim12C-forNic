#include "frag.h"
#include "fragment.h"
#include <valarray>
#include "random.h"


struct prop
{
  float Erel;
  float Dvelocity;
  float plfTheta;
  float plfPhi;
  float plfVel;
  float thetaEmission;
};


/**
 *!\brief selects the veloity vectors of the secondary fragments
 */

class CDecay
{
 public:


  bool einstein;
  static CRandom ran;
  CFrame *real[5]; 
  CFrame *recon[5];
  CFrame *plfRecon;
  CFrame *partCM[2];
  CFrag *frag[2];  //!< information about the decay fragments
  static float const pi; //!< 3.14159....

  CDecay(CFrag*,CFrag*,bool einstein0, float vCM);
  ~CDecay();
  float getErelReal();
  float getErelRecon();
  float getErel(CFrame**);
  float getErelNewton(CFrame**);
  float getErelRel(CFrame**);

  void Mode(float Et);
  void Mode(float Et, float phi);
  void Mode(float Et, float phi, float alpha);
  void ModeMicroCanonical();
  void micro(int,CFrame**,float,float);
  bool OnTopOf();
  bool OnTopOf2();
  bool OnTopOf3();
  bool OnTopOf4();
  bool OnTopOf5();
  bool OnTopOf6();
  bool OnTopOf7();

  prop Recon;

  float  sumA;
  float ErelRecon; //!<reconstructed relative kinetic energy




  int Nsolution; //!< number of p-p pairs with correct 6Be energy
  int Isolution; //!< solution #

  float angle[180];
  float phiAngle[360];
  float theta;
  float deltaPhi;
  float theta_reactionCM; // projectile angle in reaction com
  float vCM; // reaction center of mass velocity
};

