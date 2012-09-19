/*
 * Residue.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: marchi
 */

#include "Residue.h"
	int Residue::nnr=0;
	int Residue::nr=0;
	vector<int> * Residue::ind=NULL;
	int * Residue::resid=NULL;
	int * Residue::resid0=NULL;
	string * Residue::label=NULL;

Residue::Residue() {
	// TODO Auto-generated constructor stub

}
void Residue::Select(const AtomIndex & y){
	int * nvec=new int [nr];
	for(int i=0;i<nr;i++) nvec[i]=-1;
	for(int n=0;n<y.getN();n++) nvec[y[n]]=resid[y[n]];
	for(int i=0;i<nr;i++) resid[i]=nvec[i];
	delete [] nvec;
}
void Residue::SetUp(const int natoms, const int nres, const vector<int> * y
		,const vector<int> * y0, const string * l){
	nr=natoms;
	nnr=nres;
	ind=new vector<int> [nnr];
	label=new string [nnr];
	resid=new int [nr];
	resid0=new int [nr];
	for(int i=0;i<nr;i++) resid[i]=-1;
	for(int i=0; i< nnr;i++) {ind[i]=y[i]; label[i]=l[i];}

	for(int n=0; n< nnr; n++)
		for(int i=0;i<y0[n].size();i++){
			resid0[y0[n][i]]=n;
		}
	for(int n=0; n< nnr; n++)
		for(int i=0;i<ind[n].size();i++){
			resid[ind[n][i]]=n;
		}
}

Residue::~Residue() {
	// TODO Auto-generated destructor stub
}

