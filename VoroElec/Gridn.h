/*
 * Gridn.h
 *
 *  Created on: Jun 6, 2011
 *      Author: marchi
 */

#ifndef GRIDN_H_
#define GRIDN_H_
//#include "vec.h"
#include "Array.h"
#include "fftw++.h"
#include <cstdlib>
#include <cmath>
#include "MyVec.h"
#include <vector>
#include <string>
#include <sstream>
#include "DiffCoeffs.h"
typedef float Real;
using namespace Array;
using namespace fftwpp;
using std::vector;
using std::string;
using std::cout;
using std::endl;
#define VERTEX 8

template <unsigned int Ndim>
class Gridn: public array4<double> {
protected:
	static array3<Complex> * filter;
	static array3<Complex> * filterrc;
	std::string name;
	static matrix co,oc;
	static double Volume;
	static unsigned int nnx,nny,nnz;
	Real qq[VERTEX][Ndim];
	static const rvec cube[VERTEX];
	static void setMetric(const matrix &);
	double Alpha();
	double Beta();
	double Gamma();
	Real TriDiff(Real g, Real cubi){return((g-cubi>=0.0)? cubi-g+1 : g-cubi+1);};
	void TriLinear(Real, Real, Real, Real [Ndim]);
public:
	static void setNN(unsigned int nx0, unsigned int ny0,unsigned int nz0){nnx=nx0;nny=ny0;nnz=nz0;};
	static int SetNo;
	static const rvec * Getcube(){return cube;};
	Gridn();
	Gridn(Gridn &);
	Gridn(const Gridn &);
	Gridn(array3<Complex> &);
	Gridn(array4<Complex> &);
	Gridn(array4<double> &);
	Gridn(array3<double> &);
	Gridn(const double &);
	Gridn & operator=(const double);
	Gridn & operator=(const array3<double> &);
	Gridn & operator=(const array3<Complex> &);
	Gridn & operator=(const array4<Complex> &);
	Gridn & operator=(const array4<double> &);
	Gridn operator*(const double &);
	Gridn operator/(const double &);
	Gridn operator/(Gridn &);
	Gridn operator-(Gridn &);
	Gridn operator+(Gridn &);
	static void set(const matrix & CO,unsigned int nx0=0,unsigned int ny0=0,unsigned int nz0=0){
		if(nx0+ny0+nz0) setNN(nx0,ny0,nz0);setMetric(CO);
	};

	void Allocate(){
		if(!this->allocated) array4<double>::Allocate(nnx,nny,nnz,Ndim);
		else {
			this->Deallocate();
			array4<double>::Allocate(nnx,nny,nnz,Ndim);
		}
	};
	double getVol(){return Volume;};
	static matrix & getCO(){return co;};
	static matrix & getOC(){return oc;};
	static double getDV(){return Volume/static_cast<double>(nnx*nny*nnz);};
	static unsigned int getnnx(){return nnx;};
	static unsigned int getnny(){return nny;};
	static unsigned int getnnz(){return nnz;};
	unsigned int getn0(){return Ndim;};
	unsigned int getnx(){return ny;};
	unsigned int getny(){return nz;};
	unsigned int getnz(){return nw;};
	void Xplor(FILE *,double=1.0);
	vector<double> Rdf(FILE * fp,const double [DIM], std::string, const double & = 2.0, const double & = 0.02);
	void RdfK(FILE * fp, std::string, const double & = 200.0, const double & = 0.02);

	Gridn<1> Derivative(const unsigned int);
	Gridn<DIM> Der1();
	Gridn<DIM> Der2();
	Gridn<1> Div();
	Gridn<DIM> Field(){
		Gridn<DIM> temp=this->Der1();
		return -temp;
	};
	Gridn<DIM> Field2(){
		Gridn<DIM> temp=this->Der2();
		return -temp;
	};
	void setname(std::string & nme){name=nme;};
	static void MakeFilter(int, int, int);
	static void MakeFilter(const double);
	void Filter();
	template <unsigned int Ndima>
	friend Gridn<Ndima>  operator-(const Gridn<Ndima> &);
	template <class T,unsigned int Ndima>
	friend Gridn<Ndima>  operator*(const T &,const Gridn<Ndima> &);
	virtual ~Gridn();
};

#include "Gridn.cpp"
#endif /* GRIDN_H_ */
