#include "pixels.h"


pixels::pixels()
{
  ifstream file("pixels.dat");


  double a,b,c,x,y,z;
  double theta1 = -2.82725; 
  double theta2 = 1.57277;
  int e,f;
  for (int itele=0;itele<14;itele++)
    for (int ifront=0;ifront<32;ifront++)
     for  (int iback=0;iback<32;iback++)
      {
         file >> x >> y >> z >> e >> f;


         //x = -18.3083;
         //y = -29.1029;
         //z = 8.5688;

         //take away "target position" in inches
         x -= 15.2146;
         y -= -18.2034;
         z -= 8.63831;


         x *= 2.54;
         y *= 2.54;
         z *= 2.54;

          
         double xx = x*cos(theta1) + y*sin(theta1);
         double yy = y*cos(theta1) - x*sin(theta1);


         //cout << xx << " " << yy << " " << z << endl;

         double zz = z*cos(theta2) + xx*sin(theta2);
         double xxx = xx*cos(theta2) - z*sin(theta2);

	 //cout << xxx << " " << yy << " " << zz << endl;

         TeleP[itele].Location[ifront][iback].x = -yy;
         TeleP[itele].Location[ifront][iback].y = -xxx;
         TeleP[itele].Location[ifront][iback].z = zz-.79-.05;

      }
  file.close();
  file.clear();

  location center1 = getCenter(6);
  cout << "det6" << endl;
  cout << center1.x << " " << center1.y << " " << center1.z << endl;
  cout << sqrt(pow(center1.x,2)+pow(center1.y,2)+pow(center1.z,2)) << endl;
  location center2 = getCenter(7);

  //float xx = (center1.x + center2.x)/2.;
  //float yy = (center1.y + center2.y)/2.;
  //float zz = (center1.z + center2.z)/2.;

  for (int itele=0;itele<14;itele++)
    for (int ifront=0;ifront<32;ifront++)
     for  (int iback=0;iback<32;iback++)
      {

	//TeleP[itele].Location[ifront][iback].z -= zz;
	//TeleP[itele].Location[ifront][iback].x -= xx;
	//TeleP[itele].Location[ifront][iback].y -= yy;
         float r = pow(TeleP[itele].Location[ifront][iback].x,2)
           + pow(TeleP[itele].Location[ifront][iback].y,2)
	   + pow(TeleP[itele].Location[ifront][iback].z,2);
	 r = sqrt(r);
         TeleP[itele].Location[ifront][iback].theta =
	   acos(TeleP[itele].Location[ifront][iback].z/r);
         TeleP[itele].Location[ifront][iback].phi =
	   atan2(TeleP[itele].Location[ifront][iback].y,
		 TeleP[itele].Location[ifront][iback].x);
      }

  /*
  cout << "det #6 " << endl;
  cout << TeleP[6].Location[0][0].x << " " <<TeleP[6].Location[0][0].y
       << " " << TeleP[6].Location[0][0].z << endl;
  cout << TeleP[6].Location[0][32].x << " " <<TeleP[6].Location[0][32].y
       << " " << TeleP[6].Location[0][32].z << endl;
  cout << TeleP[6].Location[32][0].x << " " <<TeleP[6].Location[32][0].y
       << " " << TeleP[6].Location[32][0].z << endl;
  cout << TeleP[6].Location[32][32].x << " " <<TeleP[6].Location[32][32].y
       << " " << TeleP[6].Location[32][32].z << endl;
  */

}
//*********************************
location pixels::getCenter(int itele)
{
  float x = (TeleP[itele].Location[15][15].x 
          + TeleP[itele].Location[15][16].x
          + TeleP[itele].Location[16][15].x
          + TeleP[itele].Location[16][16].x)/4.;

  float y = (TeleP[itele].Location[15][15].y
          + TeleP[itele].Location[15][16].y
          + TeleP[itele].Location[16][15].y
          + TeleP[itele].Location[16][16].y)/4.;

  float z = (TeleP[itele].Location[15][15].z 
          + TeleP[itele].Location[15][16].z
          + TeleP[itele].Location[16][15].z
          + TeleP[itele].Location[16][16].z)/4.;
 
  location out;
 out.x = x;
 out.y = y;
 out.z = z;

 return out;
}
//***************************************
float pixels::getAngle(int itele, int ifront, int iback)
{
  phi = TeleP[itele].Location[ifront][iback].phi;
  return TeleP[itele].Location[ifront][iback].theta;
  
}
//*********************************
float pixels::getCsiCenter(int itele, int iCsi)
{
  float x,y,z;
  if (iCsi == 0)
    {
     x = (TeleP[itele].Location[7][7].x 
          + TeleP[itele].Location[7][8].x
          + TeleP[itele].Location[8][7].x
          + TeleP[itele].Location[8][8].x)/4.;

     y = (TeleP[itele].Location[7][7].y
          + TeleP[itele].Location[7][8].y
          + TeleP[itele].Location[8][7].y
          + TeleP[itele].Location[8][8].y)/4.;

     z = (TeleP[itele].Location[7][7].z 
          + TeleP[itele].Location[7][8].z
          + TeleP[itele].Location[8][7].z
          + TeleP[itele].Location[8][8].z)/4.;
    }
  if (iCsi == 1)
    {
     x = (TeleP[itele].Location[7][23].x 
          + TeleP[itele].Location[7][24].x
          + TeleP[itele].Location[8][23].x
          + TeleP[itele].Location[8][24].x)/4.;

     y = (TeleP[itele].Location[7][23].y
          + TeleP[itele].Location[7][24].y
          + TeleP[itele].Location[8][23].y
          + TeleP[itele].Location[8][24].y)/4.;

     z = (TeleP[itele].Location[7][23].z 
          + TeleP[itele].Location[7][24].z
          + TeleP[itele].Location[8][23].z
          + TeleP[itele].Location[8][24].z)/4.;
    }
  if (iCsi == 2)
    {
     x = (TeleP[itele].Location[23][23].x 
          + TeleP[itele].Location[23][24].x
          + TeleP[itele].Location[24][23].x
          + TeleP[itele].Location[24][24].x)/4.;

     y = (TeleP[itele].Location[23][23].y
          + TeleP[itele].Location[23][24].y
          + TeleP[itele].Location[24][23].y
          + TeleP[itele].Location[24][24].y)/4.;

     z = (TeleP[itele].Location[23][23].z 
          + TeleP[itele].Location[23][24].z
          + TeleP[itele].Location[24][23].z
          + TeleP[itele].Location[24][24].z)/4.;
    }
  if (iCsi == 3)
    {
     x = (TeleP[itele].Location[23][7].x 
          + TeleP[itele].Location[23][8].x
          + TeleP[itele].Location[24][7].x
          + TeleP[itele].Location[24][8].x)/4.;

     y = (TeleP[itele].Location[23][7].y
          + TeleP[itele].Location[23][8].y
          + TeleP[itele].Location[24][7].y
          + TeleP[itele].Location[24][8].y)/4.;

     z = (TeleP[itele].Location[23][7].z 
          + TeleP[itele].Location[23][8].z
          + TeleP[itele].Location[24][7].z
          + TeleP[itele].Location[24][8].z)/4.;
    }

  float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return acos(z/r)*180./3.14159;
}
//*****************************************************

void pixels::prepareSim()
{
  for (int itele = 0;itele < 14;itele++)
    {

      location out = getCenter(itele);
      Tele[itele].r_center[0] = out.x;
      Tele[itele].r_center[1] = out.y;
      Tele[itele].r_center[2] = out.z;


      Tele[itele].r_front[0] = TeleP[itele].Location[31][0].x - TeleP[itele].Location[0][0].x;
      Tele[itele].r_front[1] = TeleP[itele].Location[31][0].y - TeleP[itele].Location[0][0].y;
      Tele[itele].r_front[2] = TeleP[itele].Location[31][0].z - TeleP[itele].Location[0][0].z;

      Tele[itele].r_back[0] = TeleP[itele].Location[0][31].x - TeleP[itele].Location[0][0].x;      
      Tele[itele].r_back[1] = TeleP[itele].Location[0][31].y - TeleP[itele].Location[0][0].y;      
      Tele[itele].r_back[2] = TeleP[itele].Location[0][31].z - TeleP[itele].Location[0][0].z;      

      for (int i=0;i<3;i++)
	{
	  Tele[itele].r_front[i] *=32./31./2.;
	  Tele[itele].r_back[i] *=32./31./2.;
	}

    }

}

//************************************************
bool pixels::sim(float theta, float phi, float xtarget, float ytarget, float rand1, float rand2)
{
  //cout << theta << endl;
  bool hitit = false;
  sle Sle(3);
  for (int itele=0;itele<14;itele++)
    {
     Sle.M[0][0] = sin(theta)*cos(phi);
     Sle.M[0][1] = -Tele[itele].r_front[0];
     Sle.M[0][2] = -Tele[itele].r_back[0];

     Sle.M[1][0] = sin(theta)*sin(phi);
     Sle.M[1][1] = -Tele[itele].r_front[1];
     Sle.M[1][2] = -Tele[itele].r_back[1];

     Sle.M[2][0] = cos(theta);
     Sle.M[2][1] = -Tele[itele].r_front[2];
     Sle.M[2][2] = -Tele[itele].r_back[2];

     Sle.Y[0] = Tele[itele].r_center[0] - xtarget;
     Sle.Y[1] = Tele[itele].r_center[1] - ytarget;
     Sle.Y[2] = Tele[itele].r_center[2];

     Sle.solve();




     float dist_front = Sle.Y[1];
     float dist_back = Sle.Y[2];


      if (fabs(dist_front) < 1. && fabs(dist_back) < 1.) 
      hitit = true;

      if (hitit == false) continue;;


      hitTele = itele;

      //find strip hit
      ixStrip = (int)((dist_front + 1.)*32./2.);
      iyStrip = (int)((dist_back + 1.)*32./2.);

      float xRecon = ((float)ixStrip + rand1)/32.*2. - 1.;
       float yRecon = ((float)iyStrip +rand2)/32.*2. - 1.;

      //for 2x2 array of Csi
       if (ixStrip < 16)
          {
           if (iyStrip < 16) ICsI = 0;
           else ICsI = 1;
           }
       else 
          {
           if (iyStrip < 16) ICsI = 3;
           else ICsI = 2;
           }



  float rRecon[3];
  float rr = 0.;
  for (int i=0;i<3;i++) 
    {
      rRecon[i] = Tele[itele].r_center[i] + xRecon*Tele[itele].r_front[i] + yRecon*Tele[itele].r_back[i];
      rr += pow(rRecon[i],2);
    }




  float r = sqrt(rr);
  thetaRecon = acos(rRecon[2]/r);
  phiRecon = atan2(rRecon[1],rRecon[0]);



  return hitit;


        
    }
  return hitit;

}
