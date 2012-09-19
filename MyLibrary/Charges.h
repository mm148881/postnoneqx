/*
 * Charges.h
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */

#include "types/simple.h"
#ifndef CHARGES_H_
#define CHARGES_H_

#include "Atoms.h"

typedef float Real;

class Charges: public Atoms {
protected:
	Real * q;
	void copy(const Charges &);
public:
	Charges(): Atoms(){q=NULL;};
	Charges(const int nint):
		Atoms(nint){
		q=NULL;
		if(nint)
			q=new Real [nr];
	};
	Charges(const AtomIndex & id):
		Atoms(id){
		nr=id.getN();
		q=NULL;
		if(nr)
			q=new Real [nr];
	};
	Charges(const Charges &);
	virtual ~Charges();
	virtual Charges & operator=(const Charges & y){copy(y);return *this;};
	virtual Charges & operator()(const Charges & y){copy(y);return *this;};
	virtual void setQ(const Real *, const int );
	virtual void setQ(const Real *, const AtomIndex &);
	virtual Real getQ(const int i) const {return q[i];};
	bool test(){if(q==NULL) return FALSE;else return TRUE;}
};

#endif /* CHARGES_H_ */
