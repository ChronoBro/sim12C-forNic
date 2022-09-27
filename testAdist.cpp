#include "adist.h"
#include <iostream>
using namespace std;


int main()
{
  adist Adist;

  for (int i=0;i<100;i++)
    {
      double theta = ((double)i+0.5)*20./180.*acos(-1.)/100.;

      cout << Adist.getDist(theta,50.,10.) << endl;
    }
}

