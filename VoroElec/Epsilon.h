/*
 * Epsilon.h
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#ifndef EPSILON_H_
#define EPSILON_H_
#include "MyUtilClass.h"
#include "Array.h"
#include "Polarization.h"
#include "Field.h"
#include "Euler.h"
#include "Communicator.hpp"
#include "Voronoi.h"
#include "Residue.h"

using namespace Array;
using namespace DVECT;
using namespace MATRIX;
using namespace PolarizationNS;

namespace EpsilonNS {
struct HistData{
	Matrix dm;
	Matrix gm;
	Dvect  e;
	Dvect  p;
	int idx;
	HistData(){
		dm=0.0;
		gm=0.0;
		e=0.0;
		p=0.0;
		idx=0;
	}
};
struct IndexedMatrix{
	int i;
	Matrix dm;
	Matrix gm;
	Dvect p;
	Dvect e;
};

class Epsilon {
protected:
	bool is_COset,is_xRefset;
	Matrix co,oc;
	double Volume;
	Dvect M_avg;
	int TotalCount;
	static unsigned int nnx,nny,nnz;
	static void setNN(int nx,int ny, int nz){
		nnx=nx,nny=ny;nnz=nz;
	}
	virtual void WriteIt(std::ofstream &) const=0;
	virtual void ReadIt(std::ifstream &)=0;
public:
	Epsilon():is_COset(false),is_xRefset(false),co(0),oc(0),M_avg(0.0),TotalCount(0){};
	Epsilon(int nx,int ny,int nz):is_COset(false),is_xRefset(false),co(0)
		,oc(0), M_avg(0),TotalCount(0){
		setNN(nx,ny,nz);
	}
	virtual int GetNresid(){return 0;};
	virtual void Rdf(const double x, const double y){};
	virtual void NERdf(const double x, const double y){};
	virtual void setxRef(Dvect & x){is_xRefset=true;};
	virtual void setxRef(Dvect x){is_xRefset=true;};
	virtual void Statistics(){};
	virtual void setdir(const int D){};


	virtual void Allocate()=0;
	virtual void Deallocate()=0;
	virtual void filter(){};
	virtual Epsilon & operator-=(const Epsilon &);
	virtual Epsilon & operator=(const Epsilon &);

	double getDV() const {return Volume/static_cast<double> (nnx*nny*nnz);}
	virtual bool operator()(Polarization & x, Field & y, Matrix & M, Voronoi * z){return false;};

	void setMetric(Matrix &);
	void setMetric(matrix & co_in0){
		Matrix co_in=co_in0;
		setMetric(co_in);
	}

	static void Setup(int nx,int ny, int nz){setNN(nx,ny,nz);};

	Matrix getCO(){return co;};
	Matrix getOC(){return oc;};
	unsigned int GetNx(){return nnx;};
	unsigned int GetNy(){return nny;};
	unsigned int GetNz(){return nnz;};
	virtual ~Epsilon(){};
	friend std::ofstream & operator<<(std::ofstream & ofout, const Epsilon & y){
		y.WriteIt(ofout);
		return ofout;
	};
	friend std::ifstream & operator>>(std::ifstream & ifin, Epsilon & y){
		y.ReadIt(ifin);
		return ifin;
	}

};

} /* namespace EpsilonNS */
#endif /* EPSILON_H_ */
