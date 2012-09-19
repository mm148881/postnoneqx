/*
 * MyStrings.h
 *
 *  Created on: Dec 26, 2011
 *      Author: marchi
 */

#ifndef MYSTRINGS_H_
#define MYSTRINGS_H_
#include <string>
#include <iostream>
#include <sstream>
using std::string;

class MyStrings {
	string sep;
	string ext;
public:
	MyStrings(const string e,const string s="."): ext(e),sep(s){};
	string  operator()(const string);
	string  operator()(const string, const int );
	virtual ~MyStrings();
};

#endif /* MYSTRINGS_H_ */
