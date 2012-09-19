/*
 * ChangeInput.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: marchi
 */

#include "ChangeInput.h"

namespace Parallel {

ChangeInput::ChangeInput(int argci, char ** argvi, int proc, int nproc) {
	vector<string> a(argvi,argvi+argci);
	vector<string>::iterator its,its3,its2;
	its=find(a.begin(),a.end(),"-b");
	its3=find(a.begin(),a.end(),"-e");
	its2=find(a.begin(),a.end(),"-dt");
	try{
		if(its2 == a.end()) throw " Ouch!! can't find -dt option in input. How can I do my job in parallel??";
		if(its == a.end())  throw " Ouch!! can't find -b option in input. How can I do my job in parallel??";
		if(its3 == a.end())  throw " Ouch!! can't find -e option in input. How can I do my job in parallel??";
	}
	catch(const char * s){
		cout << s << endl;
		exit(1);
	}
	argc=argci;
	argv=new char * [argci];
	++its;
	++its3;
	++its2;
	double Beg0=Conv<double>()(*its);
	double Beg;
	double End=Conv<double>()(*its3);
	double dtt=Conv<double>()(*its2);
	double tmp=(End-Beg0)/nproc;
	Beg=Beg0+proc*tmp;
	End=Beg0+(proc+1)*tmp-dtt;
	*its=Conv<double>()(Beg);
	*its3=Conv<double>()(End);

	for(int i=0;i<argc;i++){
		argv[i]=new char [a[i].size()];
		strcpy(argv[i], a[i].c_str());
	}
}

ChangeInput::~ChangeInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace Parallel */
