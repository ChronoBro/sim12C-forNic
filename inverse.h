#include <iostream>
using namespace std;

class inverse
{
 public: 
  inverse(int n);
  ~inverse();
  int n;
  double **A;
  double **A0;
  double *V;
  double *B;
  double **AI;
  double ** AP;
  double *X;
  double D,S,XA,XB,XG,XL;


  double Sign(double);
  void S1000();
  void S2000();
  void S3000();
  void S4000();
  void solve();

};
