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

  
}

// returns 1 if particle hits detector, else returns zero
// inputs  theta = emission zenith angle in radians of scattered particle
//         phi =   emission azimuth angle in radians of scattered particle
//         phi goes from -pi to pi
int ring::hit(float theta, float phi)
{
  hitRing=-1;
  hitPie=-1;
  CsIstrip=17;
  above_lineS2=0;

  float x = ran.Rndm();
  float y = ran.Rndm();
  float dead;
  /*
  if(!S2)
    {
      dead = 1.0/(2*acos(-1)*rmax);
      
      if (y<dead) return (0); // getting rid of section between pies, just input "angular difference"/360 //integer division gives int back
      if (y>(1-dead)) return (0);
    }
  */


  //specify ring and pie section which were here

  if(S2)
    {
      float d = dist*tan(theta);
      float r = sqrt(pow(dist,2.)+pow(d,2.));
      float y1 = r*sin(theta)*sin(phi);
      if(y1 > 29.4){above_lineS2=1;please=0;return 0;}
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

      if(S2){if(hitRing>=44){return 0;}}//cout << "hitRing >=44 theta = " << theta*180/3.14159 << endl;return 0;}}
      // if(!S2){if(hitRing<=1)return 0;}

 
      if(phi>=0)
	hitPie = int(phi/dphi);
      else
	hitPie = int(phi/dphi)-1;
      if(S2){if(hitPie==3||hitPie==4)return 0;} 
      thetaHit = atan2((x+(float)hitRing)*dr+rmin,dist);
      //thetaHit = theta; //testing resolution determination
      //thetaHit = dtheta*(x+(float)hitRing)+theta_min;
      //if(S2)
      //	phiHit = (101.25+22.5*(hitPie))*acos(-1)/180.;
      //else
      phiHit =dphi*(y+(float)hitPie);
      //phiHit = phi; //testing resolution determination


      if(phi>=0)
	CsIstrip = (int) hitPie*16/Npie;
      else
	CsIstrip = (int) hitPie*16/Npie;

      
	

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
