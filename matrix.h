#ifndef matrix_
#define matrix_
#include <iostream>
#include <math.h>
using namespace std;

class matrix
{
 public:
  matrix(int rows0,int columns0);
  ~matrix();
  float ** mult(float ** A);
  float ** add(float ** A);
  float ** sub(float **A);
  float * vect_mult(float * v);
  float ** scale(float num);
  float ** C;
  float ** I;
  float ** M;

 private:
  float *row_vect;
  float *col_vect;
  int i,j,k,l,m;
  int size;
  float use;
  int rows;
  int columns;


};
#endif
