#include <fstream>
#include <iostream>
#include <cmath>
#include "sle.h"

using namespace std;

struct location
{
  double x,y,z;
  double theta,phi;
};


struct teles
{
  float r_front[3];
  float r_back[3];
  float r_center[3];
};


struct teleP
{
  location Location[32][32];
};

class pixels
{
 public:
  pixels();
  teleP TeleP[14];
  location getCenter(int itele);
  float getCsiCenter(int itele, int iCsi);
  float phi;
  float getAngle(int itele,int ifront, int iback);
  void prepareSim();
  teles Tele[14];
  bool sim(float,float,float,float, float,float);
  
  int ixStrip;
  int iyStrip;
  int ICsI;
  float thetaRecon;
  float phiRecon;
  int hitTele;
};
