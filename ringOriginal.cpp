#include "ring.h"

CRandom ring::ran;  //random number generator

//constructor
ring::ring(float dist0, float rmin0, float rmax0, int Nring0, int Npie0)
{
  dist = dist0;  // distance of ring counter to target
  rmin = rmin0;  // minimum active radius of ring counter
  rmax = rmax0;  // maximum active radius of ring counter
  Nring = Nring0; // number of annular ring sectors
  Npie  = Npie0;  //number of pie sectors
  theta_min = atan2(rmin,dist);
  theta_max = atan2(rmax,dist);
  dr = (rmax-rmin)/(float)Nring;
  dtheta = (theta_max - theta_min)/(float)Nring;
  dphi = 2*acos(-1)/(float)Npie;
  vect = new float[3];

  /*
  K = new matrix(3,3);
  K2 = new matrix(3,3);
  R = new matrix(3,3);

  float x1=4.;
  float y1=-2.;
  
  float L[3] ={0,0,358.}; //center of S2
  float B[3] = {x1,y1,358.}; // beam vector
  
  float theta_off = acos(pow(358.,2.)/(pow(358.,2.)*(sqrt(pow(B[0],2.)+pow(B[1],2.)+pow(B[2],2.)))));
  
  float cross[3] = {B[1]*L[2]-L[1]*B[2],B[0]*L[2]-L[0]*B[2],B[0]*L[1]-L[0]*B[1]};
  float cross_mag = sqrt(pow(cross[0],2.)+pow(cross[1],2.)+pow(cross[2],2.));
  float k[3] = {cross[0]/cross_mag,cross[1]/cross_mag,cross[2]/cross_mag};
  
  K->C[0][1]=-k[2];
  K->C[0][2]=k[1];
  K->C[1][2]=-k[0];
  K->C[1][0]=-1.*K->C[0][1];
  K->C[2][0] = -1.*K->C[0][2];
  K->C[2][1]=-1.*K->C[1][2];
  
  K2->C = K->C;
  K2->C = K2->mult(K2->C);
  
  R->C = R->I;
  K->C = K->scale(sin(theta_off));
  K2->C = K2->scale(1-cos(theta_off));
  R->C = R->add(K->C);
  cout << "here1" << endl;
  R->C = R->add(K2->C);
  */

}

// returns 1 if particle hits detector, else returns zero
// inputs  theta = emission zenith angle in radians of scattered particle
//         phi =   emission azimuth angle in radians of scattered particle
//         phi goes from -pi to pi
int ring::hit(float theta, float phi)
{
  float x = ran.Rndm();
  float y = ran.Rndm();
  float dead = 1.0/(2*acos(-1)*rmax);
  if (y<dead) return (0); // getting rid of section between pies, just input "angular difference"/360 //integer division gives int back
  if (y>(1-dead)) return (0);
  //specify ring and pie section which were hit

  if(S2)
    {
      float d = dist*tan(theta);
      float r = sqrt(pow(dist,2.)+pow(d,2.));
      float y1 = r*sin(theta)*sin(phi);
      if(y1 > 29.4){please=0;return 0;}
    }
      
  if (theta_min<theta && theta<theta_max)
    {
      if(theta<0.001) { hitRing = int((theta-theta_min)/dtheta);}
      else
	{
	  for (useful=1;useful < (Nring+1); useful++)
	    {
	      if(atan2(useful*dr+rmin,dist)> theta && theta > atan2((useful-1)*dr+rmin,dist))
		{ hitRing=useful-1;break;}
	    }	    
	}
      /*
      cout << "initial theta = " << theta << endl;
      cout << "initial phi = " << phi << endl;
      float radius = sqrt(pow(dist,2.)+pow(dr*(float)hitRing+rmin,2.));
				      

      vect[0] = radius*sin(theta)*cos(phi);
      vect[1] = radius*sin(theta)*sin(phi);
      vect[2] = radius*cos(theta);
				      
      vect = R->vect_mult(vect);

      radius = sqrt(pow(vect[0],2.)+pow(vect[1],2.));
      theta = atan2(radius,dist);
      phi = atan2(vect[1],vect[0]);
      */
 
      if(phi>=0)
	hitPie = int(phi/dphi);
      else
	hitPie = int(phi/dphi)-1;
      thetaHit = atan2((x+(float)hitRing)*dr+rmin,dist);
      //thetaHit = dtheta*(x+(float)hitRing)+theta_min;
      phiHit = dphi*(y+(float)hitPie);
     


      if(phi>=0)
	CsIstrip = (int) hitPie*16/Npie;
      else
	CsIstrip = (int) hitPie*16/Npie;

      
      // if(S2&&phiHit < (.168456+3.1415927/2.)&& phiHit >(3.1415927/2.-.168456)&&ran.Rndm()>0.75){please=0;return 0;}
	

      //cout << "theta = " << theta<< endl; //for testing
      //cout << "thetaHit = " << thetaHit << endl;
      //cout << "# of Rings = " << Nring << endl;
      //cout << "hitRing = " << hitRing << endl;
      please=1;
      return 1;
    }
  else
    please = 0;
    return (0);
}
