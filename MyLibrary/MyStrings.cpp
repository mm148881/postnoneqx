/*
 * MyStrings.cpp
 *
 *  Created on: Dec 26, 2011
 *      Author: marchi
 */

#include "MyStrings.h"

string MyStrings::operator()(const string filename){
	string basename;
	string output;
	string::size_type idx=filename.find(sep);
	if (idx == string::npos) {
		output=filename + sep + ext;
		return output;
    } else {
		basename = filename.substr(0, idx);
		output = basename+sep+ext;
	}
	return output;
}
string MyStrings::operator()(const string filename, const int n){
	string basename, extname, tmpname;
	string output;
	std::ostringstream Mymsg;
	Mymsg << n;
	string::size_type idx=filename.find(sep);
	if (idx == string::npos) {
		output=filename + Mymsg.str() + sep + ext;
		return output;
    } else {
		basename = filename.substr(0, idx);
		output = basename + Mymsg.str() + sep + ext;
	}
	return output;
}
MyStrings::~MyStrings() {
	// TODO Auto-generated destructor stub
}

