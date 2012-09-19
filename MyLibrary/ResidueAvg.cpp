/*
 * ResidueAvg.cpp
 *
 *  Created on: Sep 21, 2011
 *      Author: marchi
 */

#include "ResidueAvg.h"
int ResidueAvg::counter=0;

ResidueAvg::ResidueAvg(int n) {
	nr=n;
	nr0=static_cast<int> (nr);
	Avg.Allocate(nr,DIM);
	Coord.Allocate(nr,DIM);
	Avg2.Allocate(nr,nr,DIM,DIM);
	Rms.Allocate(nr,nr,DIM,DIM);
	Avg=0.0; Avg2=0.0; Rms=0.0;
};
ResidueAvg::ResidueAvg() {
	// TODO Auto-generated constructor stub
}

ResidueAvg & ResidueAvg::operator+=(ResidueCM & y){
	try{
		if(y.Size() != nr0) throw "Sizes do not match !! ";
	}
	catch(const char * s) {
		cout <<" \n" << s << endl;
		cout << y.Size() << " " << nr0  << endl;
		exit(1);
	}
	for(int i=0;i<nr0; i++){
		for(int o=0;o<DIM;o++) Avg[i][o]+=y[i][o];
		for(int j=i;j<nr0; j++){
			for(int o=0;o<DIM;o++)Avg2[i][j][o][o]+=y[i][o]*y[j][o];
		}
	}
	counter++;
	return *this;
}
void ResidueAvg::Spit(){
	cout << nr0 << endl;
}
void ResidueAvg::Average(){
	double A[DIM];
	double A2;

	for(int i=0;i<nr0; i++){
		for(int o=0;o<DIM;o++) {
			Coord[i][o]=Avg[i][o]/static_cast<double> (counter);
			A[o]=Coord[i][o];
		}
		for(int j=i;j<nr0; j++){
			for(int o=0;o<DIM;o++) {
				A2=Avg2[i][j][o][o]/static_cast<double> (counter);
				Rms[i][j][o][o]=A2-A[o]*A[o];
			}
		}
	}
}
ResidueAvg::~ResidueAvg() {
	// TODO Auto-generated destructor stub
}

