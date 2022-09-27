#include "adist.h"


adist::adist()
{
  // Leg = new legendre(100);
}
//*****************************************
adist::~adist()
{
  //delete Leg;
}

double adist::getDist(double theta, double L0, double sigL)
{
  int lmin = (int)(L0 - 3.*sigL);
  if (lmin < 0) lmin = 0;
  int lmax = (int)(L0 + 3.*sigL);
  if (lmax > 100)  
    {
      cout << "lmax=  " << lmax << " resetting to 100" << endl; 
    }

  out0 = 0.;
  out1 = 0.;
  double sum = 0.;
  lmin = 0;
  double p0[101];
  double p1[101];
  double p2[101];
  double x = cos(theta);

  for (int l= 0;l<= lmax;l++)
    {

      if (l == 0) 
	{
         p0[0] = 1.;
	 p1[0] = 0.;
	}
      else if (l==1) 
	{
         p0[1] = x;
	 p1[1] = -sqrt(1.-x*x);
	}
      else 
	{
         p0[l] = ((double)(2*l-1)*x*p0[l-1] - (double)(l-1)*p0[l-2])/(double)l; 
	 p1[l] = (double)l*(x*p0[l]-p0[l-1])/sqrt(1.-x*x);
         p2[l] = ((double)(l-1)*x*p1[l] - (double)(l+1)*p1[l-1])/sqrt(1.-x*x);
	}
      double fact = exp(-pow(((double)l-L0)/sigL,2)/2.);

       double y0 = p0[l];
       double y1 = p1[l];
       double y2 = p2[l];


       //cout << l << " " << y0 <<  " " << y1 << " " << theta << endl;
      sum += fact;
      out0 += fact*y0;
      out1 += fact*y1;
      out2 += fact*y2;
    }

  out0 /= sum;
  out1 /= sum;
  out2 /= sum;

  return out0;

}
