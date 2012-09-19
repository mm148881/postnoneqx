/*
 * EpsofK.h
 *
 *  Created on: Jul 12, 2011
 *      Author: marchi
 */

#ifndef EPSOFK_H_
#define EPSOFK_H_
#include "Array.h"
#include "typedefs.h"
#include "fftw++.h"
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include "Rho.h"
#include "RhoDip.h"
#include "Phi.h"


const double TEMP=300;
const double AVOGAD=6.0221367e23;
const double EPS0=8.854e-12;
const double BOLTZ0=1.380658e-23;
const double UNITL=1.0e-9; // nanometers!!
const double ELECHG=1.602e-19;
const double UNITC=4.0*M_PI*EPS0*UNITL/(ELECHG*ELECHG);
const double EFACT=AVOGAD/UNITC/1000.0;
const double KT=BOLTZ0*AVOGAD*TEMP/1000.0;

using namespace Array;

class EpsofK: public Rho {
protected:
	virtual void SofK(Rho &);
	virtual void SofK(Rho &, Rho &);
	virtual void SofK(Phi &, Phi &);
	virtual void SofK(RhoDip &);
public:
	EpsofK(){(*this).Allocate();};
	EpsofK(Rho & x){(*this).Allocate();SofK(x);};
	EpsofK(Rho & x, Rho & y){(*this).Allocate();SofK(x,y);};
	EpsofK(Phi & x, Phi & y){(*this).Allocate();SofK(x,y);};
	EpsofK(RhoDip & x){(*this).Allocate();SofK(x);};

	void operator()(Rho & x){SofK(x);};
	void operator()(Rho & x, Rho & y){SofK(x,y);};
	void operator()(Phi & x, Phi & y){SofK(x,y);};
	void operator()(RhoDip & x){SofK(x);};

	virtual ~EpsofK();
};

#endif /* EPSOFK_H_ */
