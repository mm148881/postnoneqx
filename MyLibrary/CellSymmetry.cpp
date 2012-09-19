/*
 * CellSymmetry.cpp
 *
 *  Created on: Jun 8, 2011
 *      Author: marchi
 */

#include "CellSymmetry.h"
t_trxframe * CellSymmetry::fr=NULL;

void CellSymmetry::Center1(AtomIndex & cidx){
	Metric Met(fr->box);
	const matrix & co=Met.getCO();
	const matrix & oc=Met.getOC();
	int natoms=fr->natoms;
	rvec * x=new rvec [natoms];
	rvec * xa=new rvec [natoms];
	for(int i=0;i<natoms;i++)
		for(int m=0;m<DIM;m++)
			x[i][m]=fr->x[i][m];

	for(int i=0;i<natoms;i++){
		xa[i][XX]=oc[XX][XX]*x[i][XX]+oc[XX][YY]*x[i][YY]+oc[XX][ZZ]*x[i][ZZ];
		xa[i][YY]=oc[YY][XX]*x[i][XX]+oc[YY][YY]*x[i][YY]+oc[YY][ZZ]*x[i][ZZ];
		xa[i][ZZ]=oc[ZZ][XX]*x[i][XX]+oc[ZZ][YY]*x[i][YY]+oc[ZZ][ZZ]*x[i][ZZ];
	}

	rvec Mycma={0,0,0};
	real mtot=0.0;
	for(int p=0;p<cidx.getN();p++){
		int i=cidx[p];
		real mm=top->atoms.atom[i].m;
		Mycma[XX]+=mm*xa[i][XX];
		Mycma[YY]+=mm*xa[i][YY];
		Mycma[ZZ]+=mm*xa[i][ZZ];
		mtot+=mm;
	}
	for(int m=0;m<DIM;m++) Mycma[m]/=mtot;

	for(int i=0;i<natoms;i++){
		xa[i][XX]+=0.5-Mycma[XX];
		xa[i][YY]+=0.5-Mycma[YY];
		xa[i][ZZ]+=0.5-Mycma[ZZ];
	}
	int nmol=top->mols.nr;

	rvec * cmol=new rvec[nmol];
	for(int n=0;n<nmol;n++){
		cmol[n][XX]=0.0 ;cmol[n][YY]=0.0 ;cmol[n][ZZ]=0.0 ;
		real mtot=0.0;
		for(int i=top->mols.index[n];i< top->mols.index[n+1];i++){
			real mm=top->atoms.atom[i].m;
			cmol[n][XX]+=mm*xa[i][XX];
			cmol[n][YY]+=mm*xa[i][YY];
			cmol[n][ZZ]+=mm*xa[i][ZZ];
			mtot+=mm;
		}
		cmol[n][XX]/=mtot;
		cmol[n][YY]/=mtot;
		cmol[n][ZZ]/=mtot;
	}
	rvec xm;
	for(int i=0;i<nmol;i++){
		xm[XX]=-rint(cmol[i][XX]-0.5);
		xm[YY]=-rint(cmol[i][YY]-0.5);
		xm[ZZ]=-rint(cmol[i][ZZ]-0.5);
		for(int ia=top->mols.index[i];ia< top->mols.index[i+1];ia++){
			xa[ia][XX]=xa[ia][XX]+xm[XX];
			xa[ia][YY]=xa[ia][YY]+xm[YY];
			xa[ia][ZZ]=xa[ia][ZZ]+xm[ZZ];

		}
	}
	for(int i=0;i<natoms;i++){
		fr->x[i][XX]=co[XX][XX]*xa[i][XX]+co[XX][YY]*xa[i][YY]+co[XX][ZZ]*xa[i][ZZ];
		fr->x[i][YY]=co[YY][XX]*xa[i][XX]+co[YY][YY]*xa[i][YY]+co[YY][ZZ]*xa[i][ZZ];
		fr->x[i][ZZ]=co[ZZ][XX]*xa[i][XX]+co[ZZ][YY]*xa[i][YY]+co[ZZ][ZZ]*xa[i][ZZ];
	}
}
