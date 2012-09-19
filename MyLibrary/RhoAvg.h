/*
 * RhoAvg.h
 *
 *  Created on: Jun 8, 2011
 *      Author: marchi
 */

#ifndef RHOAVG_H_
#define RHOAVG_H_
#include "Rho.h"
#include "RhoDip.h"
#include <iostream>
#include <typeinfo>
#include <cstring>
using std::cout;
using std::cerr;
using std::endl;

template <class T>
class RhoAvg {
protected:
	int n;
	T * avg;
public:
	RhoAvg(){
		try{
			if(strstr(typeid(T).name(),"Rho") == NULL && strstr(typeid(T).name(),"RhoDip") == NULL)
				throw " Can use RhoAvg class template only with classes Rho and RhoDip";
		}
		catch(const char * s ){
			cerr << s << endl;
			exit(1);
		}
		n=0;
		avg=NULL;
	}
	void Accumulate(T & y){
		if(avg) (*avg)+=y;
		else avg=new T(y);
		n++;
	}
	T getavg(){
		T temp(*avg);
		temp=temp/static_cast<double>(n);
		return temp;
	}
	int getN(){return n;}
	virtual ~RhoAvg(){avg->T::~T();delete avg;};
};
#include "RhoAvg.cpp"
#endif /* RHOAVG_H_ */
