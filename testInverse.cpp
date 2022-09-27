#include "inverse.h"
#include <iostream>
using namespace std;

int main()
{
  inverse Inverse(4);
  Inverse.A[0][0] = 1.;
  Inverse.A[0][1] = 0.5;
  Inverse.A[0][2] = .333333;
  Inverse.A[0][3] = 0.25;


  Inverse.A[1][0] = .5;
  Inverse.A[1][1] = 0.333333;
  Inverse.A[1][2] = .25;
  Inverse.A[1][3] = 0.2;


  Inverse.A[2][0] = .333333;
  Inverse.A[2][1] = 0.25;
  Inverse.A[2][2] = .2;
  Inverse.A[2][3] = 0.166667;


  Inverse.A[3][0] = .25;
  Inverse.A[3][1] = 0.2;
  Inverse.A[3][2] = .166667;
  Inverse.A[3][3] = 0.142857;

  Inverse.solve();

    for (int i=0;i<4;i++)
    {
      for (int j=0;j<4;j++) cout << Inverse.AI[i][j] << " " ;
       cout << endl;
    }

    Inverse.S3000();

}
