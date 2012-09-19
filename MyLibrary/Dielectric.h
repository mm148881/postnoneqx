/*
 * Dielectric.h
 *
 *  Created on: Jun 17, 2011
 *      Author: marchi
 */

#ifndef DIELECTRIC_H_
#define DIELECTRIC_H_
#include "Array.h"
#include "typedefs.h"
#include "fftw++.h"
#include "Grid.h"
#include <cmath>
#include <vector>
#include <string>
#include "Rho.h"

using namespace Array;
using namespace fftwpp;
using std::cout;
using std::endl;
using std::string;

class Dielectric {
	const double eps;
	static unsigned int nnx,nny,nnz;
	Grid<DIM> epsilon;
	void Allocate(unsigned int nx,unsigned int ny,unsigned int nz){
		nnx=nx;nny=ny;nnz=nz;
		epsilon[XX].Allocate(nnx,nny,nnz);
		epsilon[YY].Allocate(nnx,nny,nnz);
		epsilon[ZZ].Allocate(nnx,nny,nnz);
	};
public:
	Dielectric(): eps(1.0e-5){};
	Dielectric(unsigned int nx,unsigned int ny,unsigned int nz): eps(1.0e-5){
		nnx=nx;nny=ny;nnz=nz;
		epsilon[XX].Allocate(nnx,nny,nnz);
		epsilon[YY].Allocate(nnx,nny,nnz);
		epsilon[ZZ].Allocate(nnx,nny,nnz);
	}

	bool Allocated(){bool a=nnx+nny+nnz != 0; return a;};
	void GetEps(FILE *, string, Grid<DIM> &, Grid<DIM> &);
	void GetEps(FILE *, string, Rho  &, Rho &);
	void GetEps(FILE *, string, Rho &);
	bool isSmall(Complex &);
	void Print();
	virtual ~Dielectric();

};

#endif /* DIELECTRIC_H_ */
