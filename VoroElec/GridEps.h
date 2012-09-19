/*
 * GridEps.h
 *
 *  Created on: Feb 27, 2012
 *      Author: marchi
 */

#ifndef GRIDEPS_H_
#define GRIDEPS_H_
#include "MyUtilClass.h"
#include "Array.h"
#include "Polarization.h"
#include "Field.h"
#include "Euler.h"
#include "Communicator.hpp"

using namespace Array;
using namespace DVECT;
using namespace MATRIX;
using namespace PolarizationNS;

namespace GridEpsNS {
struct IndexedMatrix{
	int i;
	Matrix dm;
	Matrix gm;
	Dvect p;
	Dvect e;
};

class GridEps {
protected:
	bool is_COset,is_xRefset;
	Matrix co,oc;
	double Volume;
	Dvect xRef;
	array3<Matrix> D;
	array3<Matrix> G;
	array3<Dvect> E;
	array3<Dvect> M;
	Dvect M_avg;
	array3<Matrix> eps;
	static unsigned int nnx,nny,nnz;
	int TotalCount;
	static void setNN(int nx,int ny, int nz){
		nnx=nx,nny=ny;nnz=nz;
	}
public:
	virtual double getDV() const {return Volume/static_cast<double> (nnx*nny*nnz);}
	void Allocate(){
		try{if(!(nnx || nny || nnz)) throw " Trying to allocate, but at least "
				"one of the dimensions is thero ";}
		catch(const char * s){
			std::cout << s << std::endl;
			exit(1);
		}

		D.Allocate(nnx,nny,nnz);
		G.Allocate(nnx,nny,nnz);
		M.Allocate(nnx,nny,nnz);
		E.Allocate(nnx,nny,nnz);
		eps.Allocate(nnx,nny,nnz);
		D=0.0;
		G=0.0;
		M=0.0;
		E=0.0;
		eps=0.0;
	}
	GridEps():is_COset(false),is_xRefset(false),co(0),oc(0),xRef(0)
		,TotalCount(0), M_avg(0.0){
		try{
			if(nnx && nny && nnz) Allocate();
			else throw " Warning GridEps Allocation is differred !!!";
		}
		catch(const char * s) {
			std::cout << s << std::endl;
		}
	}
	GridEps(int nx,int ny,int nz):is_COset(false),is_xRefset(false),co(0)
		,oc(0),xRef(0),TotalCount(0), M_avg(0.0){
		setNN(nx,ny,nz);
		Allocate();
	}

	bool operator()(Polarization & , Field & , bool=false);
	void setMetric(Matrix &);
	void setMetric(matrix & co_in0){
		Matrix co_in=co_in0;
		setMetric(co_in);
	}
	void Avg();
	void Avg2();

	static void Setup(int nx,int ny, int nz){setNN(nx,ny,nz);};

	void setxRef(Dvect & x){is_xRefset=true;xRef=x;};
	void setxRef(Dvect x){is_xRefset=true;xRef=x;};
	Matrix getCO(){return co;};
	Matrix getOC(){return oc;};
	unsigned int GetNx(){return nnx;};
	unsigned int GetNy(){return nny;};
	unsigned int GetNz(){return nnz;};
	double getD(){return D[3][2][5][0][0];};
	void Deallocate(){
		D.Deallocate();
		G.Deallocate();
		M.Deallocate();
		E.Deallocate();
		eps.Deallocate();
	}
	virtual ~GridEps(){
		D.Deallocate();
		G.Deallocate();
		M.Deallocate();
		E.Deallocate();
		eps.Deallocate();
	}
	friend std::ofstream & operator<<(std::ofstream &, const GridEps & );
	friend std::ifstream & operator>>(std::ifstream &, GridEps &);
};

} /* namespace GridEpsNS */
#endif /* GRIDEPS_H_ */
