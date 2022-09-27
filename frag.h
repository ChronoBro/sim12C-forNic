#ifndef _frag
#define _frag

#include "array.h"
#include "loss.h"
#include <string>
#include "mScat.h"
#include "random.h"
#include "frame.h"
#include "fragment.h"
#include "pixels.h"
#include "ring.h"
#include "ttt.h"
using namespace std;

/**
 *!\brief information about a single fragment and its interection with the detector
 *
 */
class CFrag : public fragment
{
 public:
  static CRandom ran;
  static float const pi; //!< 3.14159....

  float CsI_res;   

  CFrag(float,float,string,float,float);
  CFrag(float,float,string,float,float,string);
  ~CFrag();
  CArray * Array;
  CArray * shadowArray;
  int hit();
  int hit(float,float);
  int hit2(float,float);
  int hitShadow(float,float);
  int hit3();
  int hit4();
  int hit5();
  int hit6(int);
  int hit7();
  int  is_hit;
  int  is_hit2;

  void SetVelocity(float*,float*);
  void AddVelocity(float*);
  void AddVelocityRel(float*);
  float Eloss(float); // energy loss in target
  float Egain(float); // corrects for energy loss in target
  void MultiScat(float);
  bool targetInteraction(float,float);


    CLoss *loss;
    CMultScat *multScat;

    CFrame *real;      //<!real particles energy, velocity, etc
    CFrame *recon;    //<!reconstructed properties

    bool protonHole(int itower, int itele, int ifront, int iback);
    bool alphaHole(int itower, int itele, int ifront, int iback);

    bool protonHole2(int itele, int ifront, int iback);
    bool alphaHole2(int itele, int ifront, int iback);

    pixels Pixels;

    //Russian ring detector
    ring *  Ring;
    // S2
    ring * Ring2;
    ttt * TTT;
    ttt * TTT2;
    ttt * TTT3;
    ttt * TTT4;
    ring * RingCsI;
    ring * Ring2CsI;

    float thres[56];

};

#endif
