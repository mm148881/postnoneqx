/*
 * BSpline.cpp
 *
 *  Created on: Jul 21, 2011
 *      Author: marchi
 */

#include "BSpline.h"
int BSpline::order=4;

BSpline::BSpline():nders(2) {
		spline=new vector<double> [nders];
		for(int o=0;o<nders;o++) spline[0].resize(order);
	};

void BSpline::Fill(const double & w){
	Init(w);
	for(int i=3;i<order;i++) OnePass(w,i);
	Diff();
	OnePass(w,order);
}

BSpline::~BSpline() {
	for(int i=0;i<nders;i++) spline[i].clear();
	delete [] spline;
}

