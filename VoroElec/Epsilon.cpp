/*
 * Epsilon.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#include "Epsilon.h"

namespace EpsilonNS {
	unsigned int Epsilon::nnx=0,Epsilon::nny=0,Epsilon::nnz=0;

	Epsilon & Epsilon::operator-=(const Epsilon & z) {
		return *this;
	}
	Epsilon & Epsilon::operator=(const Epsilon & z) {
		return *this;
	}


	void Epsilon::setMetric(Matrix & co_in){
/* Save some time by assuming lower right part is zero */
		is_COset=true;
		for(int i=0;i<DIM;i++)
			for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];

		double tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
		oc[XX][XX]=co[YY][YY]*co[ZZ][ZZ]*tmp;
		oc[YY][XX]=0;
		oc[ZZ][XX]=0;
		oc[XX][YY]=-co[XX][YY]*co[ZZ][ZZ]*tmp;
		oc[YY][YY]=co[XX][XX]*co[ZZ][ZZ]*tmp;
		oc[ZZ][YY]=0;
		oc[XX][ZZ]=(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ])*tmp;
		oc[YY][ZZ]=-co[YY][ZZ]*co[XX][XX]*tmp;
		oc[ZZ][ZZ]=co[XX][XX]*co[YY][YY]*tmp;
		Volume=(co[XX][XX]*(co[YY][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[YY][ZZ])
		  -co[YY][XX]*(co[XX][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[XX][ZZ])
		  +co[ZZ][XX]*(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ]));
	}


} /* namespace EpsilonNS */
