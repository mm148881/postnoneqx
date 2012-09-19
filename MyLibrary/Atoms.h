/*
 * Atoms.h
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */

#ifndef ATOMS_H_
#define ATOMS_H_
#include "Metric.h"
#include "AtomIndex.h"
#include <cmath>

class Atoms {
protected:
	int nr;
	int status;
	rvec * x;
	rvec * xa;
	Metric Mt;
public:
	Atoms():nr(0), x(NULL), xa(NULL){};
	Atoms(const int);
	Atoms(const AtomIndex &);
	Atoms(const Atoms &);
	virtual ~Atoms();
	struct plane{
		real xc[4];
		int n;
		plane & operator=(rvec & xa){for(int o=0;o<DIM;o++) xc[o]=xa[o];return *this;};
		plane & operator=(real dd){xc[3]=dd;return *this;};
		plane & operator=(double dd){xc[3]=dd;return *this;};
		plane & operator=(int i){n=i;return *this;};
	} ax;
	void setDim(const int n){nr=n;if(x) delete [] x;x=new rvec [nr];if(xa) delete [] xa;xa=new rvec [nr];}
	void setCoord(const Metric &, const rvec *, const AtomIndex & );
	void setCoord(const Metric &, const rvec *);
	void setMT(const Metric &);
	Atoms &  operator=(const Atoms &);
	Atoms & operator=(const real);
	rvec & operator[](const int i)const {return x[i];};
	Atoms & operator()(const int);

	void Rot(const matrix);
	void Tra(const rvec );
	int getNR()const{return nr;};
	const Metric & getMt() const {return Mt;};
	Atoms COtoOC();
	Atoms OCtoCO();
	void doCOtoOC();
	Atoms Shift(const rvec);
	real Dist(const int,const int);
	plane VectorDist(const int,const int);
	rvec * getX(){return x;};
	rvec * getXA(){return xa;};
	virtual void pbc();
};

#endif /* ATOMS_H_ */

