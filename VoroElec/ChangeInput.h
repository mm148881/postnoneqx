/*
 * ChangeInput.h
 *
 *  Created on: Mar 13, 2012
 *      Author: marchi
 */

#ifndef CHANGEINPUT_H_
#define CHANGEINPUT_H_
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>

namespace Parallel {
using namespace std;
template <class T>
class Conv {
public:
	T operator()(string a){
		T b;
		stringstream ss;
		ss<< a;
		ss >> b;
		return b;
	}
	string operator()(T b){
		stringstream ss;
		ss<< b;
		string h=ss.str();
		return h;
	}
};


class ChangeInput {
	int argc;
	char ** argv;
public:
	ChangeInput(int , char * [], int=0, int=1);
	char ** getargv(){return argv;};
	int * getargc(){return &argc;};
	virtual ~ChangeInput();
};

} /* namespace Parallel */
#endif /* CHANGEINPUT_H_ */
