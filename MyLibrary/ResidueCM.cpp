/*
 * ResidueCM.cpp
 *
 *  Created on: Sep 21, 2011
 *      Author: marchi
 */

#include "ResidueCM.h"
#include <iostream>
using std::cout;
using std::endl;
vector<int> * ResidueCM::ind=NULL;
vector<int> * ResidueCM::iD=NULL;
int ResidueCM::nnr=0;
string * ResidueCM::label=NULL;


ResidueCM::ResidueCM(const rvec * y) {
	SetCM(y);
}
ResidueCM::ResidueCM(const int natoms,const int n, const vector<int> * y, const string * l): Atoms::Atoms(natoms) {
	nnr=n;
	ind=new vector<int> [nnr];
	label=new string [nnr];
	for(int i=0; i< nnr;i++) {ind[i]=y[i]; label[i]=l[i];}
}

void ResidueCM::Allocate(const int n){
	nnr=n;
	if(!ind) delete [] ind;
	if(!xcm) delete [] xcm;
	if(!xcma) delete [] xcma;
	if(!label) delete [] label;
	ind=new vector<int> [nnr];
	xcm=new rvec [nnr];
	xcma=new rvec [nnr];
	label=new string [nnr];
}
void ResidueCM::Allocate(){
	if(!xcm) delete [] xcm;
	if(!xcma) delete [] xcma;
	xcm=new rvec [nnr];
	xcma=new rvec [nnr];
}
ResidueCM & ResidueCM::operator ()(const int n){
	Allocate(n);
	return *this;
}
ResidueCM & ResidueCM::operator ()(const int natoms, const int n, const vector<int> * y, const string * l){
	Atoms::operator()(natoms);
	Allocate(n);
	for(int i=0; i< nnr;i++) {ind[i]=y[i];label[i]=l[i];}
	return *this;
}
ResidueCM & ResidueCM::operator ()(const ResidueCM & y){
	Allocate(y.nnr);
	for(int i=0; i< nnr;i++) ind[i]=y.ind[i];
	for(int i=0; i< nnr;i++) label[i]=y.label[i];
	for(int i=0; i< nnr;i++)
		for(int o=0;o<DIM;o++) {
			xcm[i][o]=y.xcm[i][o];
			xcma[i][o]=y.xcma[i][o];
		}
	return *this;
}
ResidueCM & ResidueCM::operator ()(const rvec * y){
	SetCM(y);
	return *this;
}
ResidueCM & ResidueCM::operator=(const float z){
	float * y=&x[0][0];
	for(int i=0;i<DIM*nr;i++) *y=z;
	return *this;
}
ResidueCM & ResidueCM::operator=(const ResidueCM & y){
	Allocate(y.nnr);
	for(int i=0; i< nnr;i++) ind[i]=y.ind[i];
	for(int i=0; i< nnr;i++) label[i]=y.label[i];
	for(int i=0; i< nnr;i++)
		for(int o=0;o<DIM;o++) {
			xcm[i][o]=y.xcm[i][o];
			xcma[i][o]=y.xcma[i][o];
		}
	return *this;
}
ResidueCM & ResidueCM::operator+=(const ResidueCM & y){
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++){
			xcm[i][o]+=y.xcm[i][o];
			xcma[i][o]+=y.xcma[i][o];
		}
	return *this;
}
ResidueCM ResidueCM::operator+(const ResidueCM & y){
	ResidueCM temp(y);
	temp+=*this;
	return temp;
}

void ResidueCM::doCOtoOC(){
	try{
	if(!xcm) throw "ResidueCM class coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	if(!xcma) xcma=new rvec [nnr];
	for(int i=0;i<nnr;i++)
		for(int o=0;o<DIM;o++) xcma[i][o]=Mt.getOC()[o][0]*xcm[i][0]+Mt.getOC()[o][1]*xcm[i][1]+
			Mt.getOC()[o][2]*xcm[i][2];
}
void ResidueCM::SetCM(const rvec * y){
	if(!xcm) xcm=new rvec [nnr];
	for(int o=0; o<nnr;o++){
		int n=ind[o].size();
		rvec xcmm={0.0,0.0,0.0};
		for(int i=0;i<n;i++) {
			int p=ind[o][i];
			xcmm[XX]+=y[p][XX];
			xcmm[YY]+=y[p][YY];
			xcmm[ZZ]+=y[p][ZZ];
		}
		for(int i=0;i<DIM;i++) xcmm[i]/=static_cast<float> (n);
		for(int i=0;i<DIM;i++) xcm[o][i]=xcmm[i];
	}

}
void ResidueCM::SetCM(){
	if(!x) x=new rvec [nnr];
	for(int o=0; o<nnr;o++){
		int n=ind[o].size();
		rvec xcmm={0.0,0.0,0.0};
		for(int i=0;i<n;i++) {
			int p=ind[o][i];
			xcmm[XX]+=xcm[p][XX];
			xcmm[YY]+=xcm[p][YY];
			xcmm[ZZ]+=xcm[p][ZZ];
		}
		for(int i=0;i<DIM;i++) xcmm[i]/=static_cast<float> (n);
		for(int i=0;i<DIM;i++) xcm[o][i]=xcmm[i];
	}

}
void ResidueCM::ResidueSelect(const int ngrps,const AtomIndex * cidx){
	if(iD) iD->clear();
	else iD=new vector<int>;
    array1<int> SelIndex;
    SelIndex.Allocate(nr);
    SelIndex=-1;
    try{
    	int sum=0;
    	for(int i=0;i<ngrps;i++) sum+=cidx[i].getN();
    	if(sum > nr) throw "The atoms of the selected groups overlap, cannot compute ";
    	if(sum < nr) throw "The selected groups do not encompass the entire system ";
    }
    catch(const char * s) {
    	cout << s << endl;
    	exit(1);
    }

    for(int i=0;i<ngrps;i++)
    	for(int o=0;o<cidx[i].getN();o++){
    		int ia=cidx[i][o];
    		SelIndex[ia]=i;
    	}


    int nr0=0;
    vector<vector<int> > ind0;
	vector<string> label0;

    for(int i=0;i<nnr;i++){
    	vector<int> * Sel=new vector<int> [ngrps];
    	for(unsigned int o=0;o<ind[i].size();o++){
    		int ia=ind[i][o];
    		Sel[SelIndex[ia]].push_back(ia);
    	}

    	for(int p=0;p<ngrps;p++){
    		if(Sel[p].empty()) continue;
    		std::stringstream io;
    		io << label[i]<<p;
    		label0.push_back(io.str());
    		ind0.push_back(Sel[p]);
    		iD->push_back(i);
    		nr0++;
    	}
    	delete [] Sel;
    }
    Allocate(nr0);
    for(int p=0;p<nnr;p++){
    	ind[p]=ind0[p];
    	label[p]=label0[p];
    }
}
void ResidueCM::pbc(){
	for(int i=0;i<nnr;i++)
		for(int o=0;o<DIM;o++)
			xcma[i][o]=xcma[i][o]-rint(xcma[i][o]);
	for(int i=0;i<nnr;i++)
		for(int o=0;o<DIM;o++) xcm[i][o]=Mt.getCO()[o][0]*xcma[i][0]+Mt.getCO()[o][1]*xcma[i][1]+
			Mt.getCO()[o][2]*xcma[i][2];

}

ResidueCM::~ResidueCM() {
	if(xcm) delete [] xcm;
	if(xcma) delete [] xcma;
	// TODO Auto-generated destructor stub
}
std::ostream & operator<<(std::ostream & out, const ResidueCM  & y){
	for(int i=0; i< y.nnr;i++) {
		for(unsigned int j=0;j<y.ind[i].size(); j++)
		out << " "  << i+1 << " " << y.ind[i][j]+1 << " " << y.label[i]<< endl;
	}
	return out;
}
