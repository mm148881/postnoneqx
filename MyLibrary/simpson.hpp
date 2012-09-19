/*
 * simpson.hpp
 *
 *  Created on: Jul 19, 2011
 *      Author: marchi
 */

#ifndef SIMPSON_HPP_
#define SIMPSON_HPP_
namespace simpson{
static double * func0=NULL;
static int ntot;
static double dxx;
static double cut;
void func(double f[], double cut0,int ntot0){
	cut=cut0;
	ntot=ntot0;
	dxx=cut/double (ntot);
	if(func0) delete [] func0;
	func0=new double[ntot];
	for(int i=0;i<ntot;i++) func0[i]=f[i];
}
double func(double x){
	int h=(x<cut)? int (x/dxx):int (cut/dxx);
	double g0=(x/dxx)- double (h);
	double g1=1.0-g0;
	double f=(h==int (cut/dxx))?func0[h]:g0*func0[h]+g1*func0[h+1];
	return f;
}
double simpson(const int N, const double A, const double B)
{
  double X, h, Iapp0, Iapp1, Iapp2, Iapp;
  int NN, i;

  // Etape 1
  h = (B - A) / N;

  // Etape 2
  Iapp0 = func(A) + func(B);
  Iapp1 = 0.0;
  Iapp2 = 0.0;

  // Etape 3
  NN = N -1;
  for (i=1; i<=NN; i++)
    {
      // Etape 4
      X = A + i*h;
      // Etape 5
      if ((i%2) == 0)
        Iapp2 = Iapp2 + func(X);
      else
        Iapp1 = Iapp1 + func(X);
    }

  // Etape 6
  Iapp = (Iapp0 + 2.0 * Iapp2 + 4.0 * Iapp1) * h / 3.0;

  // Etape 7
  return (Iapp);

}

double simpson(const double a, double b, const int n0) {
	double s, dx, x;
// if n is odd - add +1 interval to make it even
	int n=n0;
	if(n%2) n++;
	s = 0.0;
	dx = (b-a)/static_cast<float>(n);
	for ( int i=1; i<n; i++){
		x = a+static_cast<float>(i)*dx;
		s+=(i%2)?2.0*func(x):4.0*func(x);
	}
	s = (s + func(a)+func(b))*dx/3.0;

return s;
}

}

#endif /* SIMPSON_HPP_ */
