/* rint.c -- round to integral value in floating point format

   replacement for missing math routine under VMS

   H. J. Bernstein, yaya@dowling.edu, 28 Apr 02                */

#include <math.h>
double rint(double);
float rintf(float);

double rint(double x) {
  return ( (x<0.)? -floor(-x+.5):floor(x+.5) );
}

float rintf(float x) {
  return ( (x<(float)0.)? -(float)floor(-x+.5):(float)floor(x+.5) );
}
