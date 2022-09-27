#include "inverse.h"
#include <cmath>

// ----------------------------------------------------------- 
// Ref.: "Methodes de calcul numerique - Tome 2 - By Claude    
//        Nowakowski, PSI Editions, 1984".                     
//                                                             
//                  Visual C++ Version By J-P Moreau, Paris.   
//                             (www.jpmoreau.fr)               
//**************************************************************
inverse::inverse(int n0)
{
  n = n0;
  V = new double [n];
  B = new double [n];
  X = new double [n];

  A = new double * [n];
  A0 = new double * [n];
  AI = new double * [n];
  AP = new double * [n];


  for (int i=0;i<n;i++)
    {
      A[i] = new double [2*n];
      A0[i] = new double [n];
      AI[i] = new double [n];
      AP[i] = new double [n];
    }

}

//**************************************************
inverse::~inverse()
{
  delete V;
  delete B;
  delete X;
  for (int i=0;i<n;i++)
    {
      delete A[i];
      delete A0[i];
      delete AI[i];
      delete AP[i];
    }
  delete A;
  delete A0;
  delete AI;
  delete AP;
}
//*************************************
double inverse::Sign(double x) 
        {
	   double s;
	   if (x<0.0) s = -1.0;
	   else if (x>0.0) s = 1.0;
	   else s = 0.0;
	   return s;
	}



//************************************************
void inverse::S1000() 
{
   int i,j,k;

   for (k = 0; k<n-1; k++)
     {
       S = 0.0;
       for (i = k; i<n; i++)S += (A[i][k] * A[i][k]);
       XL = -Sign(A[k][k]) * sqrt(S);
       XA = XL * (XL - A[k][k]);
       V[k] = A[k][k] - XL;
       for (i = k+1 ; i<n; i++) V[i] = A[i][k];
       for (j = k+1 ; j<2*n; j++)
          {
            XB = 0.0;
            for (i = k; i<n; i++) XB += (V[i] * A[i][j]);
            XG = XB / XA;
            for (i = k; i<n; i++) A[i][j] -= (XG * V[i]);
	  }
       A[k][k] = XL;
       for (i = k+1 ; i<n; i++) A[i][k] = 0.0;

       //       cout << "transformation " << k << endl;
       //for (int i=0;i<n;i++)
       //  {
       //    for (int j=0;j<n;j++) cout << A[i][j] << " " ;
       //    cout << endl;
       //  }
       //cout << endl;
     }




} // end of S1000

//************************************
void inverse::S2000() 
{
    int i, j;
    for (i = n-1; i>=0; i--) 
      {
            j = n; S = 0.0;
e2020:      if (i == j) goto e2040;
            S -= A[i][j] * X[j];
            j--;
            goto e2020;
e2040:      X[i] = (B[i] + S) / A[i][i];
      }
 }
//*******************************************
  //check A*A-1 = identity

void inverse::S3000() 
{
    int i,j,k;
    for (i = 0; i<n; i++)
	for (j = 0; j<n; j++) 
          {
            S = 0.0;
            for (k = 0; k<n; k++)
            S += A0[i][k] * AI[k][j];
            AP[i][j] = S;
	  }

       cout << "check A*A-1 " << endl;
       for (int i=0;i<n;i++)
         {
           for (int j=0;j<n;j++) 
	     {
	       if (fabs(AP[i][j]) < 1.e-7) AP[i][j] = 0.;
                cout << AP[i][j] << " " ;
	     }
           cout << endl;
         }
       cout << endl;

}

//***********************************
//* 
void inverse::S4000() 
{
   int i;
   D = 1.0;
   for (i = 0; i<n; i++)D *= A[i][i];
}

//************************************************
  // matrix to be inverse must be in A
void inverse::solve()
{


  //   cout << "input" << endl;
  // for (int i=0;i<n;i++)
  //   {
  //     for (int j=0;j<n;j++) cout << A[i][j] << " " ;
  //     cout << endl;
  //   }
  // cout << endl;
    
    for (int i = 0; i<n; i++)
       {
        for (int j = 0; j<n; j++) 
           {
            A[i][j + n] = 0.0;
            A0[i][j] = A[i][j];
	   }
        A[i][i + n] = 1.0;
       }

        // Transform A into triangular matrix
        S1000();
        // N linear systems to solve
        for (int k = 0; k<n; k++)
          {
           for (int i = 0; i<n; i++) B[i] = A[i][k + n];
            // Solve triangular system
            S2000();
	    for (int i = 0; i<n; i++) AI[i][k] = X[i];
	  }

        // Calculate determinant
        //S4000();
        //fprintf(fp1, " AP=A*A^-1 matrix:\n");
        //S3000();
} //end of main program
