/*
 * Charges.cpp
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */

#include "Charges.h"
void Charges::copy(const Charges & y){
	try{
		if(nr != y.nr) throw "Cannot copy Charges of different size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++) {
		q[i]=y.q[i];
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
	}
	Mt(y.Mt);
}

Charges::Charges(const Charges & y): Atoms(y) {
	q=new real [y.nr];
	for(int i=0;i<nr;i++) q[i]=y.q[i];
}
Charges::~Charges() {
	if(q != NULL) {
		delete [] q;
		delete [] x;
		nr=0;
	}
}
void Charges::setQ(const real * c, const int nind){
	this->setDim(nind);
	if(q) delete [] q;
	q=new real [nr];
	for(int i=0;i<nr;i++) q[i]=c[i];
}
void Charges::setQ(const real * c, const AtomIndex & id){
	this->setDim(id.getN());
	if(q) delete [] q;
	q=new real [nr];
	for(int i=0;i<nr;i++) q[i]=c[id[i]];
}
