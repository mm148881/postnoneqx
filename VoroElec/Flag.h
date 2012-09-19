/*
 * Flag.h
 *
 *  Created on: Apr 18, 2012
 *      Author: marchi
 */

#ifndef FLAG_H_
#define FLAG_H_
#include <iostream>

class Flag {
	static int freq;
public:
	Flag();
	static void set(int y){freq=y;};
	bool operator()(int x){return freq!=0?x?!(x%freq):false:false;};
	int Freq(){return freq;}
	double fFreq(){return static_cast<double> (freq);}
	virtual ~Flag();
};

#endif /* FLAG_H_ */
