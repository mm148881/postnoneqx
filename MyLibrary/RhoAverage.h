/*
 * RhoAverage.h
 *
 *  Created on: Jul 13, 2011
 *      Author: marchi
 */

#ifndef RHOAVERAGE_H_
#define RHOAVERAGE_H_
#include "Rho.h"

class RhoAverage: public Rho {
protected:
	int n;
public:
	RhoAverage(){n=0;Allocate();(*this)[0]=0.0;

	}
	void Accumulate(){
		n++;
	}
	Rho getavg(){
		Rho temp(*this);
		temp/=static_cast<double>(n);
		return temp;
	}
	void getavgthis(){
		try{
			if(!n) throw "Cannot take average if count is zero ";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
		*this/=static_cast<double>(n);
	}
	virtual ~RhoAverage();

	RhoAverage & operator+=(const Rho & a){
		n++;
		for(unsigned int i=0;i<this->size;i++)
			this->v[i]+=a(i);
		return *this;
	}
};

#endif /* RHOAVERAGE_H_ */
