/*
 * CellSymmetry.h
 *
 *  Created on: Jun 8, 2011
 *      Author: marchi
 */

#ifndef CELLSYMMETRY_H_
#define CELLSYMMETRY_H_
#include "gmx_polar.h"
#include "AtomIndex.h"
#include "princ.h"
#include "rmpbc.h"

using std::cout;
using std::endl;

class CellSymmetry {
	t_topology *top;
	static t_trxframe *fr;
	t_pbc *pbc;
	rvec * x_ref;
	gmx_rmpbc_t  gpbc;
	CellSymmetry(){};
	void Gpbc(){
		if (gpbc != NULL)
			gmx_rmpbc_done(gpbc);
	}
public:
	CellSymmetry(t_topology * its_top, t_trxframe * its_fr,t_pbc * its_pbc){
		top=its_top;
		fr=its_fr;
		pbc=its_pbc;
		int ePBC;
		ePBC=fr->ePBC;
		set_pbc(pbc,ePBC,fr->box);
		gpbc = gmx_rmpbc_init(&top->idef,ePBC,fr->natoms,fr->box);
		x_ref=NULL;
	}
	CellSymmetry & operator()(t_topology * its_top, t_trxframe * its_fr,t_pbc * its_pbc){
		top=its_top;
		fr=its_fr;
		pbc=its_pbc;
		int ePBC;
		ePBC=fr->ePBC;
		set_pbc(pbc,ePBC,fr->box);
		gpbc = gmx_rmpbc_init(&top->idef,ePBC,fr->natoms,fr->box);
		x_ref=NULL;
		return *this;
	}
	void Xref(int n, atom_id * c){
		Center0(n,c);
		if(!x_ref) x_ref=new rvec [Polar::AtomData::natoms];
		for(int i=0;i<Polar::AtomData::natoms;i++){
			x_ref[i][XX]=fr->x[i][XX];
			x_ref[i][YY]=fr->x[i][YY];
			x_ref[i][ZZ]=fr->x[i][ZZ];
			}
	}
	void Xref(int n, atom_id * c, rvec * x0, matrix M){
		gmx_rmpbc(gpbc,fr->natoms,M,x0);
		reset_x(n,c, Polar::AtomData::natoms, NULL,x0,Polar::AtomData::w_lst);
		if(!x_ref) x_ref=new rvec [Polar::AtomData::natoms];
		for(int i=0;i<Polar::AtomData::natoms;i++){
			x_ref[i][XX]=x0[i][XX];
			x_ref[i][YY]=x0[i][YY];
			x_ref[i][ZZ]=x0[i][ZZ];
			}
	}
	void Center(int nindex, atom_id * cindex, int natoms=Polar::AtomData::natoms){
		gmx_rmpbc_trxfr(gpbc,fr);
		rvec xcm;
		calc_xcm(fr->x,nindex,cindex,top->atoms.atom,xcm,FALSE);
		for(int i=0;i<DIM;i++) xcm[i]=-xcm[i];
		add_xcm(fr->x,natoms,NULL,xcm);
	}
	void Center(){
		gmx_rmpbc_trxfr(gpbc,fr);
		rvec xcm;
		calc_xcm(fr->x,Polar::AtomData::nindex,Polar::cindex,top->atoms.atom,xcm,FALSE);
		for(int i=0;i<DIM;i++) xcm[i]=-xcm[i];
		add_xcm(fr->x,Polar::AtomData::natoms,NULL,xcm);
	}
	void CenterCM(int nresid,rvec * x){
		gmx_rmpbc(gpbc,nresid,fr->box,x);
		Gpbc();
	}
	void Center1(AtomIndex &);

	void Center0(int nindex0, atom_id * cindex0){
		gmx_rmpbc(gpbc,fr->natoms,fr->box,fr->x);
		reset_x(nindex0,cindex0, Polar::AtomData::natoms, NULL,fr->x,Polar::AtomData::w_lst);
	}
	void Center0(int nindex0, atom_id * cindex0, rvec * x0, matrix box0){
		gmx_rmpbc(gpbc,fr->natoms,box0,x0);
		reset_x(nindex0,cindex0, Polar::AtomData::natoms, NULL,x0,Polar::AtomData::w_lst);
	}

	void Mypbc(){
		gmx_rmpbc_trxfr(gpbc,fr);
	}
	double Rms(float R[3][3]){
		double diff=0.0;
		int uo=0;
		double tm=0.0;
		for(int i=0;i<Polar::AtomData::natoms;i++){
			if(!Polar::AtomData::w_lst[i]) continue;
			uo++;
			rvec xo,xn,xp;
			for(int p=0;p<DIM;p++)xo[p]=x_ref[i][p];
			for(int p=0;p<DIM;p++)xn[p]=fr->x[i][p];
			for(int o=0;o<DIM;o++){
				xp[o]=0.0;
				for(int p=0;p<DIM;p++){
					xp[o]+=R[o][p]*xn[p];
				}
			}
			for(int o=0;o<DIM;o++) diff+=(xp[o]-xo[o])*(xp[o]-xo[o]);
		}
		diff/=static_cast<double> (Polar::AtomData::nindex);
		diff=sqrt(diff);
		return diff;
	}
	double gmxRms(float R[3][3]){
		double diff=0.0;
		int uo=0;
		double tm=0.0;
		for(int i=0;i<Polar::AtomData::natoms;i++){
			if(!Polar::AtomData::w_lst[i]) continue;
			double m=Polar::AtomData::w_lst[i];
			tm+=m;
			uo++;
			rvec xo,xn,xp;
			for(int p=0;p<DIM;p++)xo[p]=x_ref[i][p];
			for(int p=0;p<DIM;p++)xn[p]=fr->x[i][p];
			for(int o=0;o<DIM;o++){
				xp[o]=0.0;
				for(int p=0;p<DIM;p++){
					xp[o]+=R[o][p]*xn[p];
				}
			}
			for(int o=0;o<DIM;o++) diff+=m*(xp[o]-xo[o])*(xp[o]-xo[o]);
		}
		diff/=tm;
		return sqrt(diff);
	}
	void GetRot(float R[3][3]){
		calc_fit_R(DIM,Polar::AtomData::natoms,Polar::AtomData::w_lst,x_ref,fr->x,R);
	}
	void GetRot(float R[3][3],rvec * pt){
		calc_fit_R(DIM,Polar::AtomData::natoms,Polar::AtomData::w_lst,x_ref,pt,R);
	}
	void Fit(){
		do_fit(Polar::AtomData::natoms,Polar::AtomData::w_lst,x_ref,fr->x);
	}

	virtual ~CellSymmetry(){top=NULL;fr=NULL;pbc=NULL;x_ref=NULL; if(x_ref) delete [] x_ref;Gpbc();};
};

#endif /* CELLSYMMETRY_H_ */
