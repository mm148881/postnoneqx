/*
 * RhoDipAverage.h
 *
 *  Created on: Jul 13, 2011
 *      Author: marchi
 */

#ifndef RHODIPAVERAGE_H_
#define RHODIPAVERAGE_H_
#include "RhoDip.h"

class RhoDipAverage: public RhoDip {
protected:
	int n;
public:
	RhoDipAverage(){n=0;Allocate();for(int m=0;m<DIM;m++) (*this)[m]=0.0;};
	void Accumulate(){
		n++;
	}
	RhoDip getavg(){
		RhoDip temp(*this);
		temp/=static_cast<double>(n);
		return temp;
	}
	RhoDipAverage operator+(RhoDipAverage & a){RhoDipAverage temp(*this); temp+=a; return temp;};

	virtual ~RhoDipAverage();
};

#endif /* RHODIPAVERAGE_H_ */
