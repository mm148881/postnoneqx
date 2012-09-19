/*
 * Field.h
 *
 *  Created on: Dec 22, 2011
 *      Author: marchi
 */

#ifndef FIELD_H_
#define FIELD_H_

#include "Gridn.h"
#include "Grid.h"

class Field: public Gridn<DIM> {
	static array3<double> * ro0;
	static array3<Complex> * rok0;
	static array3<Complex> * dgi;
	int times;
public:
	static void setup(){
		try{
			if(ro0) throw "Cannot call setup twice. Exit here.";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
		unsigned int nzp=nnz/2+1;
		size_t align=sizeof(Complex);
		ro0=new array3<double>;
		ro0->Allocate(nnx,nny,nnz,align);
		rok0=new array3<Complex>;
		rok0->Allocate(nnx,nny,nzp,align);
		dgi=new array3<Complex>[DIM];
		for(int i=0;i<DIM;i++) dgi[i].Allocate(nnx,nny,nzp,align);
	}
	Field():Gridn<DIM>::Gridn(){
		std::string gh("Fieldoutput.bin");
		name=gh;
	};
	Field & operator()(Grid<1> &);
	Field & operator()(const double y){
		for(unsigned int i=0;i<this->size;i++)
			this->v[i]=y;

		return *this;
	}
	Field & operator=(const double y){
		for(unsigned int i=0;i<this->size;i++)
			this->v[i]=y;

		return *this;
	}

	Field & operator+=(Field &);
	Field getavg(){
		try{
			if(!times) throw "Cannot take average if count is zero ";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}

		Field temp(*this);
		temp/=static_cast<double>(times);
		return temp;
	}
	void getavgthis(){
		try{
			if(!times) throw "Cannot take average if count is zero ";
		}
		catch (const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
		*this/=static_cast<double>(times);
	}
	virtual ~Field();
};

#endif /* FIELD_H_ */
