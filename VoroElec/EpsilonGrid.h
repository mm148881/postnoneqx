/*
 * EpsilonGrid.h
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#ifndef EPSILONGRID_H_
#define EPSILONGRID_H_
#include "Epsilon.h"
#include "Polarization.h"
namespace EpsilonNS {

class EpsilonGrid: public Epsilon {
protected:
	static Matrix zero;
	Dvect xRef;
	array3<Matrix> D;
	array3<Matrix> G;
	array3<Dvect> E;
	array3<Dvect> M;
	array3<Matrix> eps;
	virtual void WriteIt(std::ofstream &) const;
	virtual void ReadIt(std::ifstream &);
public:
	EpsilonGrid();
	virtual ~EpsilonGrid();
	virtual bool operator()(Polarization &, Field &, Matrix & M=zero, Voronoi * b=NULL);
	virtual EpsilonGrid & operator-=(const Epsilon &);
	virtual EpsilonGrid & operator=(const Epsilon &);


	virtual void Allocate(){
		try{if(!(nnx || nny || nnz)) throw " Trying to allocate, but at least "
				"one of the dimensions is zero ";}
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
	virtual void Deallocate(){
		D.Deallocate();
		G.Deallocate();
		M.Deallocate();
		E.Deallocate();
		eps.Deallocate();
	}
	virtual void filter(){
		Polarization::setNN(nnx,nny,nnz);
		Polarization M0,E0;
		M0.Allocate();E0.Allocate();
		double * Mt0=&M0[0][0][0][0];
		double * Mt1=&M[0][0][0][0];
		double * Et0=&E0[0][0][0][0];
		double * Et1=&E[0][0][0][0];
		for(size_t i=0;i<nnx*nny*nnz;i++){
			Mt0[i]=Mt1[i];
			Et0[i]=Et1[i];
		}
		M0.Filter();
		E0.Filter();
		for(size_t i=0;i<nnx*nny*nnz;i++){
			Mt1[i]=Mt0[i];
			Et1[i]=Et0[i];
		}

	}
	void setxRef(Dvect & x){is_xRefset=true;xRef=x;};
	void setxRef(Dvect x){is_xRefset=true;xRef=x;};
};

} /* namespace EpsilonNS */
#endif /* EPSILONGRID_H_ */
