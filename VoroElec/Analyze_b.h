/*
 * Analyze_b.h
 *
 *  Created on: Mar 1, 2012
 *      Author: marchi
 */

#ifndef ANALYZE_B_H_
#define ANALYZE_B_H_


#include "Metric.h"
#include "CellSymmetry.h"
#include "Voronoi.h"
#include "gmx_voronoi.h"
#include "Field.h"
#include "Parameters.h"
#include "gmx_voroelec.h"
#include "Gridn_SpecTemplates.hpp"
#include "Grid_SpecTemplates.hpp"
#include "Communicator.hpp"

#include "DiffCoeffs_Spec.hpp"
#include "Residue.h"
#include "Accumulate.h"
#include "Flag.h"
#include <cmath>

using namespace gmx_Voronoi;
using namespace gmx_voroelec;
const bool TEST=true;
static int nframes=0,nstart,freq;
static int myframes=0;
static CellSymmetry * symmetry=NULL;
static bool bHydrogen=FALSE;
static double dtt;
static double Mytime=0.0;
static matrix M;
using std::vector;
static Rho Rho_f;
static Polarization P_f;
static Accumulate<Dvect> LocX;
static Dvect2rvec Tr0;

void getDipole(int nres, MyDipole * Dip, rvec * x0){
	int nresid=Residue::GetNNR();
	int natoms=Residue::GetNR();
	vector<Dvect> mdip(nresid,0.0);
	vector<int> ndip(nresid,0);
	vector<double> tchg(nresid,0.0);
	for(int m=0;m<natoms;m++){
		int res=Residue::ResiduePid0(m);
		rvec & x=x0[m];
		double qchg=Dip->getQ(m);
		tchg[res]+=qchg;
		for(int o=0;o<DIM;o++) mdip[res][o]+=qchg*x[o];
		ndip[res]++;
	}
}



static int Analyze(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{

    /* Here, you can do whatever analysis your program requires for a frame. */
    Parallel::comm->Barrier();

	Voronoi * vor0=static_cast<Voronoi *> (data);

	if(!symmetry){
		symmetry=new CellSymmetry(top,fr,pbc);
		symmetry->Xref(Polar::AtomData::nindex,Polar::cindex,Polar::AtomData::x,Polar::AtomData::box);
	}

	symmetry->Center0(Polar::AtomData::nindex,Polar::cindex);
	symmetry->GetRot(M);
	Matrix Rot(M);

	Metric Met;
    Met(fr->box);
    MyEps0->setMetric(fr->box);
    if(vor0) {
    	atm.setCoord(Met,fr->x);
    	atm.doCOtoOC();
    	vor0->Start(fr->time,atm);
    }

    Grid<1>::set(fr->box);
	Gridn<DIM>::set(fr->box);

    Dips->setCoord(Met, fr->x);

    P_t->Density(*Dips);
    P_t->Filter();
    Rho_t->Density(*Dips);
    Rho_t->Filter();


    (*E_t)(*Rho_t);
    if(!(*MyEps0)(*P_t, *E_t, Rot, vor0)) exit(-1);

    nframes++;
    return 0;
}
static int Analyze_Conv(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
	nframes++;

    /* Here, you can do whatever analysis your program requires for a frame. */

	Voronoi * vor0=static_cast<Voronoi *> (data);

	if(!symmetry){
		symmetry=new CellSymmetry(top,fr,pbc);
		symmetry->Xref(Polar::AtomData::nindex,Polar::cindex,Polar::AtomData::x,Polar::AtomData::box);
		LocX.Allocate(natoms);
		Rho_f.Allocate();
		Rho_f=0.0;
		P_f.Allocate();
		P_f=0.0;

	}
	symmetry->Center0(Polar::AtomData::nindex,Polar::cindex);
	Metric Met;
    Met(fr->box);
    vector<Dvect> Tr;
    Tr.assign(fr->x,fr->x+natoms);
    LocX+=Tr;

    Grid<1>::set(fr->box);
	Gridn<DIM>::set(fr->box);

    Dips->setCoord(Met, fr->x);

    P_t->Density(*Dips);
    P_f+=*P_t;
    Rho_t->Density(*Dips);
    Rho_f+=*Rho_t;

    if(Flag()(nframes)) {
    	double fact0=1.0/Flag().fFreq();
    	P_f*=fact0;
    	Rho_f*=fact0;
    	P_f.Filter();
    	Rho_f.Filter();
    	(*E_t)(Rho_f);

    	Tr0(LocX.Avg());
    	symmetry->GetRot(M,Tr0.pt());
    	Matrix Rot(M);
    	if(vor0) {
    		atm.setCoord(Met,Tr0.pt());
    		atm.doCOtoOC();
    		vor0->Start(fr->time,atm);
    	}
    	MyEps0->setMetric(fr->box);
    	if(!(*MyEps0)(P_f, *E_t, Rot, vor0)) exit(-1);
    	P_f=0.0;
    	Rho_f=0.0;
    	LocX.zero();
    	EpsilonNS::EpsilonGrid * opGrid=dynamic_cast<EpsilonNS::EpsilonGrid *> (MyEps0);
    	EpsilonNS::EpsilonVor  * opVor=dynamic_cast<EpsilonNS::EpsilonVor *> (MyEps0);
    	if(opVor) {
    		Dvect M=opVor->TestM(4);
    		Dvect E=opVor->TestE(4);
    		E*=Parameters::unitfield;
    		Matrix D=opVor->TestD(4);
    		double Dd=(D[XX][XX]+D[YY][YY]+D[ZZ][ZZ])/3.0;
    		Matrix G=opVor->TestG(4);
    		double Gg=(G[XX][XX]+G[YY][YY]+G[ZZ][ZZ])/3.0;
    		dataIn->push_back(M[XX]);
    		dataIn->push_back(M[YY]);
    		dataIn->push_back(M[ZZ]);
    		dataIn->push_back(E[XX]);
    		dataIn->push_back(E[YY]);
    		dataIn->push_back(E[ZZ]);
    		dataIn->push_back(Dd);
    		dataIn->push_back(Gg);
    	}

    }
    return 0;
}
static int Analyze0(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{

    /* Here, you can do whatever analysis your program requires for a frame. */

	Voronoi * vor0=static_cast<Voronoi *> (data);

	if(!symmetry){
		symmetry=new CellSymmetry(top,fr,pbc);
		symmetry->Xref(Polar::AtomData::nindex,Polar::cindex,Polar::AtomData::x,Polar::AtomData::box);
	}

	symmetry->Center0(Polar::AtomData::nindex,Polar::cindex);
	symmetry->GetRot(M);
	Matrix Rot(M);

	Metric Met;
    Met(fr->box);
    MyEps0->setMetric(fr->box);


    Grid<1>::set(fr->box);
	Gridn<DIM>::set(fr->box);

    Dips->setCoord(Met, fr->x);

    P_t->Density(*Dips);
    Dvect PP0=P_t->GetTotDip();
    Dvect PP=Rot*P_t->GetTotDip();
    nframes++;
    dataIn->push_back(PP0[XX]);
    dataIn->push_back(PP0[YY]);
    dataIn->push_back(PP0[ZZ]);
    dataIn->push_back(PP[XX]);
    dataIn->push_back(PP[YY]);
    dataIn->push_back(PP[ZZ]);
    return 0;
}
/*
static int Analyze_Conv(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{

     Here, you can do whatever analysis your program requires for a frame.
    nframes++;

	Voronoi * vor0=static_cast<Voronoi *> (data);

	if(!symmetry){
		symmetry=new CellSymmetry(top,fr,pbc);
		symmetry->Xref(Polar::AtomData::nindex,Polar::cindex,Polar::AtomData::x,Polar::AtomData::box);
		Rho_f.Allocate();
		Rho_f=0.0;
		P_f.Allocate();
		P_f=0.0;
	}

	symmetry->Center0(Polar::AtomData::nindex,Polar::cindex);
	symmetry->GetRot(M);
	Matrix Rot(M);

	Metric Met;
    Met(fr->box);
    MyEps0->setMetric(fr->box);
    vector<Dvect> Tr;
    Tr.assign(fr->x,fr->x+natoms);
    LocX+=Tr;
    if(Flag(nframes)) {
    	cout << " Here we are " << endl;
    	exit(-10);
    }
    if(vor0) {
    	atm.setCoord(Met,fr->x);
    	atm.doCOtoOC();
    	vor0->Start(fr->time,atm);
    }

    Grid<1>::set(fr->box);
	Gridn<DIM>::set(fr->box);

    Dips->setCoord(Met, fr->x);

    P_t->Density(*Dips);
    P_f+=*P_t;
    Rho_t->Density(*Dips);
    Rho_f+=*Rho_t;

    double dx=Rho_f.getCO()[ZZ][ZZ]/static_cast<double> (Rho_f.getnnz());

    cout <<  " Illa 0 " << endl;
    Phi Myphi(*Rho_t);

   // for(size_t o=0;o<Rho_t->getnnx();o++) cout << o*dx << " " << Rho_f[0][o][0][0]/static_cast<double>(nframes) << endl;
    cout << Rot << endl;

    (*E_t)(*Rho_t);
    if(!(*MyEps0)(*P_t, *E_t, Rot, vor0)) exit(-1);
    EpsilonNS::EpsilonGrid * opGrid=dynamic_cast<EpsilonNS::EpsilonGrid *> (MyEps0);
    EpsilonNS::EpsilonVor  * opVor=dynamic_cast<EpsilonNS::EpsilonVor *> (MyEps0);
    if(opVor) {
    	Dvect M=opVor->TestM(4);
    	Dvect E=opVor->TestE(4);
    	Matrix D=opVor->TestD(4);
    	double Dd=(D[XX][XX]+D[YY][YY]+D[ZZ][ZZ])/3.0;
    	Matrix G=opVor->TestG(4);
    	double Gg=(G[XX][XX]+G[YY][YY]+G[ZZ][ZZ])/3.0;
    	dataIn->push_back(M[XX]);
    	dataIn->push_back(M[YY]);
    	dataIn->push_back(M[ZZ]);
    	dataIn->push_back(E[XX]);
    	dataIn->push_back(E[YY]);
    	dataIn->push_back(E[ZZ]);
    	dataIn->push_back(Dd);
    	dataIn->push_back(Gg);
    }
    return 0;
}
*/
#endif /* ANALYZE_B_H_ */
