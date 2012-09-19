/*
 * Accumulate.h
 *
 *  Created on: Apr 18, 2012
 *      Author: marchi
 */

#ifndef ACCUMULATE_H_
#define ACCUMULATE_H_
#include <vector>
using std::vector;

template <class T>
class Accumulate {
	vector<T> Ac;
	int count;
public:
	Accumulate();
	Accumulate(size_t n):count(0){Allocate(n);};
	void Allocate(size_t n){Ac=vector<T>(n,T(0));};
	void zero(){count=0;Ac.assign(Ac.size(),T(0));}
	T & operator[](size_t n){return Ac[n];}
	Accumulate<T> & operator+=(vector<T> & y){
		for(std::size_t n=0;n<y.size();n++)
			(*this)[n]+=y[n];
		count++;
		return *this;
	}
	Accumulate<T> & operator()(vector<T> & y){
		for(std::size_t n=0;n<y.size();n++)
			(*this)[n]+=y[n];
		count++;
		return *this;
	}
	Accumulate<T> & operator()(T * begin, std::size_t length){
		T * it=begin;
		for(std::size_t n=0;n<length;n++){
			(*this)[n]+=*(it+n);
		}
		count++;
		return *this;
	}
	inline T  at(int n){
		return Ac[n]/static_cast<double> (count);
	}
	vector<T> & Avg(){
		for(size_t n=0;n<Ac.size();n++)
			Ac[n]=at(n);
		return Ac;
	}
	virtual ~Accumulate();
};
#include "Accumulate.cpp"
#endif /* ACCUMULATE_H_ */
