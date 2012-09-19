/*
 * Dipoles.h
 *
 *  Created on: Jun 21, 2011
 *      Author: marchi
 */

#ifndef DIPOLES_H_
#define DIPOLES_H_

#include "Atoms.h"
#include "AtomIndex.h"
#include <iostream>
#include <cmath>
using std::cout;
using std::endl;
class Dipoles: public Atoms {
	const double unitd;
	int ndip;
	int natoms;
	rvec * dip;
	real * chg;
	real * mass;
	AtomIndex * cidx;
	int * at_nr;
public:
	Dipoles();
	virtual ~Dipoles();
	void set(const int, const int,const bool *,const int *,const real *, const real * );
	void setCoord(const Metric &, const rvec *, AtomIndex &);
	rvec & gdip(int i){return dip[i];}
	rvec & operator()(int i){return dip[i];}
	int getCDIX(int i, int j) const {
		if(j < cidx[i].getN())
			return cidx[i].getI(j);
		else return -1;
	};
	bool test(){if(dip==NULL) return FALSE;else return TRUE;}
};

#endif /* DIPOLES_H_ */
