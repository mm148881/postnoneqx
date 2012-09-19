/*
 * Dipoles.cpp
 *
 *  Created on: Jun 21, 2011
 *      Author: marchi
 */

#include "Dipoles.h"

Dipoles::Dipoles():unitd(1.0){ndip=0; cidx=NULL; chg=NULL; mass=NULL;};

void Dipoles::set(const int natom0,const int ndip0, const bool * Solv, const int * molidx,
		const real * q, const real * mm){
	ndip=0;
	for(int i=0;i<ndip0;i++) if(Solv[i]) ndip++;
	natoms=natom0;
	if(!dip) dip=new rvec [ndip];
	if(!at_nr) at_nr=new int [ndip];
	if(!cidx) cidx=new AtomIndex [ndip];
	if(!chg) chg=new real [natoms];
	if(!mass) mass=new real [natoms];
	for(int i=0;i<natoms;i++) chg[i]=q[i];
	for(int i=0;i<natoms;i++) mass[i]=mm[i];

	int nn=0;
	for(int i=0;i<ndip0;i++) {
		if(Solv[i]) {
			int ncindex0;
			(i+1<ndip0)?
					ncindex0=molidx[i+1]-molidx[i]:
					ncindex0=natoms-molidx[i];
			int * cindex0=new int [ncindex0];
			int n=0;
			int m=molidx[i];
			for(int o=m;o<m+ncindex0;o++){
				cindex0[n++]=o;
			}
			cidx[nn].set(ncindex0,cindex0);
			at_nr[nn]=ncindex0;
			nn++;
		}
	}
	setDim(ndip);
}
void Dipoles::setCoord(const Metric & Mt_in, const rvec * xm, AtomIndex & ecidx){
	Mt(Mt_in);
	rvec * x0=new rvec [natoms];
	for(int i=0;i<natoms;i++) {
		x0[i][XX]=0.0;
		x0[i][YY]=0.0;
		x0[i][ZZ]=0.0;
	}
	int n=ecidx.getN();
	for(int m=0;m<n;m++) {
		int i=ecidx[m];
		x0[i][XX]=xm[i][XX];
		x0[i][YY]=xm[i][YY];
		x0[i][ZZ]=xm[i][ZZ];
	}


	for(int i=0;i<ndip;i++){
		for(int o=0;o<DIM;o++) {
			dip[i][o]=0.0;
			x[i][o]=0.0;
		}
		real masstot=0.0;
		n=cidx[i].getN();
		for(int o=0;o<n;o++){
			int m=cidx[i][o];
			real q=chg[m];
			real mm=mass[m];
			masstot+=mass[m];
			x[i][XX]+=mm*x0[m][XX];
			x[i][YY]+=mm*x0[m][YY];
			x[i][ZZ]+=mm*x0[m][ZZ];
			dip[i][XX]+=x0[m][XX]*q;
			dip[i][YY]+=x0[m][YY]*q;
			dip[i][ZZ]+=x0[m][ZZ]*q;
			}
		x[i][XX]/=masstot;
		x[i][YY]/=masstot;
		x[i][ZZ]/=masstot;
		dip[i][XX]*=unitd;
		dip[i][YY]*=unitd;
		dip[i][ZZ]*=unitd;
	}
	delete [] x0;
}

Dipoles::~Dipoles() {
	// TODO Auto-generated destructor stub
}
