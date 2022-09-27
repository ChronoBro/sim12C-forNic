#include "stdio.h"
#include "math.h"
#include <iostream>

using namespace std;

#define SIZE 25
   
int n;    //size of matrix A
double A[SIZE][2*SIZE], A0[SIZE][SIZE], V[SIZE];
double D, S, XA, XB, XG, XL;
double B[SIZE], AI[SIZE][SIZE], AP[SIZE][SIZE], X[SIZE];
FILE *fp, *fp1; 

    double Sign(double x) {
		double s;
		if (x<0.0) s = -1.0;
		else if (x>0.0) s = 1.0;
		else s = 0.0;
		return s;
	}

    void Print_A(int n) {
      int i, j;
        for (i = 1; i<=n; i++)
          for (j = 1; j<=n; j++)
            if (j == n)
              fprintf(fp1, "  %f\n", A[i][j]);
            else
			  fprintf(fp1, "  %f", A[i][j]);
    }

    void S1000() {
		int i,j,k;
		cout << endl;
        for (k = 1; k<n; k++) {
            S = 0.0;
            for (i = k; i<=n; i++)
                S += (A[i][k] * A[i][k]);
            XL = -Sign(A[k][k]) * sqrt(S);
            XA = XL * (XL - A[k][k]);
            V[k] = A[k][k] - XL;
            for (i = k + 1; i<=n; i++)
                V[i] = A[i][k];
            for (j = k + 1; j<=2*n; j++) {
                XB = 0.0;
                for (i = k; i<=n; i++)
                    XB += (V[i] * A[i][j]);
                XG = XB / XA;
                for (i = k; i<=n; i++)
                    A[i][j] -= (XG * V[i]);
			}
            A[k][k] = XL;
            for (i = k + 1; i<=n; i++)
                A[i][k] = 0.0;
            cout << " Transformation #" << k << endl;
       for (int i=1;i<=n;i++)
         {
           for (int j=1;j<=n;j++) cout << A[i][j] << " " ;
           cout << endl;
         }
       cout << endl;



 		} // end of k loop
	} // end of S1000

    void S2000() {
        int i, j;
        for (i = n; i>0; i--) {
            j = n; S = 0.0;
e2020:      if (i == j) goto e2040;
            S -= A[i][j] * X[j];
            j--;
            goto e2020;
e2040:      X[i] = (B[i] + S) / A[i][i];
		}
	}

    void S3000() {
		int i,j,k;
        for (i = 1; i<=n; i++)
			for (j = 1; j<=n; j++) {
                S = 0.0;
                for (k = 1; k<=n; k++)
                    S += A0[i][k] * AI[k][j];
                AP[i][j] = S;
			}
	}

    void S4000() {
		int i;
        D = 1.0;
        for (i = 1; i<=n; i++)
            D *= A[i][i];
	}


int main() {
	    int i,j,k;
        // read matrix to be inverted
        // Open input/output text files
        //fp = fopen("householder.dat", "r");
        //fp1 = fopen("householder.txt", "w");
        //fprintf(fp1,"\n");
        //fprintf(fp1, " Input square matrix:\n");
        //fscanf(fp,"%d", &n);
        //for (i = 1; i<=n; i++) 
	//		for (j = 1; j<=n; j++)
        //        fscanf(fp,"%lf", &A[i][j]);
        //fclose(fp);
	//	Print_A(n);
	    n = 4;
  A[1][1] = 1.;
  A[1][2] = 0.5;
  A[1][3] = .333333;
  A[1][4] = 0.25;


  A[2][1] = .5;
  A[2][2] = 0.333333;
  A[2][3] = .25;
  A[2][4] = 0.2;


  A[3][1] = .333333;
  A[3][2] = 0.25;
  A[3][3] = .2;
  A[3][4] = 0.166667;


  A[4][1] = .25;
  A[4][2] = 0.2;
  A[4][3] = .166667;
  A[4][4] = 0.142857;

		for (i = 1; i<=n; i++) {
            for (j = 1; j<=n; j++) {
                A[i][j + n] = 0.0;
                A0[i][j] = A[i][j];
			}
            A[i][i + n] = 1.0;
		}
        // Transform A into triangular matrix
        S1000();
        // N linear systems to solve
        for (k = 1; k<=n; k++) {
            for (i = 1; i<=n; i++)
                B[i] = A[i][k + n];
            // Solve triangular system
            S2000();
			for (i = 1; i<=n; i++)
                AI[i][k] = X[i];
		}
        fprintf(fp1,"\n");
        fprintf(fp1, " Inverted matrix:\n");
		 for (i = 1; i<=n; i++)
            for (j = 1; j<=n; j++) {
                if (j == n)
                    fprintf(fp1, "  %f\n", AI[i][j]);
                else
					fprintf(fp1, "  %f", AI[i][j]);
			}
        // Calculate determinant
        S4000();
        fprintf(fp1,"\n");
        fprintf(fp1, " Determinant: %e\n", D);
        // Check AP=A*A^-1=Unity Matrix
		fprintf(fp1,"\n");
        fprintf(fp1, " AP=A*A^-1 matrix:\n");
        S3000();
        for (i = 1; i<=n; i++)
            for (j = 1; j<=n; j++) {
                if (j == n)
                    fprintf(fp1, "  %f\n", AP[i][j]);
                else
					fprintf(fp1, "  %f", AP[i][j]);
			}
        fclose(fp1);
		printf("\n See results in file householder.txt...\n\n");

		return 0;
} //end of main program
