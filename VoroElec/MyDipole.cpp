/*
 * MyDipole.cpp
 *
 *  Created on: Dec 27, 2011
 *      Author: marchi
 */

#include "MyDipole.h"
#include <iomanip>
vector<int> MyDipole::resid;
vector<int> MyDipole::residDip;
int MyDipole::nresid=0;

void MyDipole::setCoord(const Metric & Mt_in, const rvec * x0){
	this->Atoms::setCoord(Mt_in,x0);
	this->doCOtoOC();

	vector<Dvect> dip0(nresid,0.0);
	for(int i=0;i<ndipoles;i++){
		dip[i]=0.0;
		for(int j=0;j<Conf[i].cidx.size();j++){
			int ia=Conf[i].cidx[j];
			Myreal mass=Conf[i].mass[j];
			for(int m=0;m<DIM;m++)  dip[i].x[m]+=static_cast<double> (mass*x[ia][m]);
			for(int m=0;m<DIM;m++) dip[i].xa[m]+=static_cast<double> (mass*xa[ia][m]);
		}

		rvec xcm;

		for(int m=0;m<DIM;m++) xcm[m]=dip[i].x[m];

		double sum=0.0;
		for(int j=0;j<Conf[i].cidx.size();j++){
			int ia=Conf[i].cidx[j];
			for(int m=0;m<DIM;m++) dip[i].d[m]+=static_cast<double> ((q[ia]-Conf[i].q)*(x[ia][m]-xcm[m]));
		}
	}
	for(int i=0;i<ndipoles;i++)
		for(int m=0;m<DIM;m++)
			dip0[residDip[i]][m]+=dip[i].d[m];
}
MyDipole::MyDipole(const MyDipole & y){
	natoms=y.natoms;
	ndipoles=y.ndipoles;
	dip=new Dip[y.ndipoles];
	Conf=new DipG[y.ndipoles];
	for(int i=0;i<ndipoles;i++) dip[i]=y.dip[i];
	for(int i=0;i<ndipoles;i++) Conf[i]=y.Conf[i];
}

MyDipole::MyDipole(const int nato, const int ndip, const int * cindex,
		const Myreal * chg, const Myreal * mass, const vector<int> * y,int nresid0){
	natoms=int(nato); ndipoles=int(ndip);
	dip=new Dip[ndip];
	Conf=new DipG[ndip];
	this->setQ(chg,natoms);

	nresid=nresid0;
	resid=vector<int>(natoms,0);
	for(int i=0;i<nresid;i++)
		for(std::size_t j=0;j<y[i].size();j++)
			resid[y[i][j]]=i;

	residDip=vector<int>(ndip,0);
	for(int i=0;i<ndip;i++){
		int end=ndip-i-1?cindex[i+1]:natoms;
		Myreal tmass=0.0;
		Myreal tchg=0.0;
		int m=0;
		for(int j=cindex[i];j<end;j++){
			tmass+=mass[j];
			tchg+=q[j];
			Conf[i].cidx.push_back(j);
			m++;
		}
		Conf[i].q=tchg/static_cast<Myreal> (m);
		Conf[i].tmass=tmass;
		for(int j=cindex[i];j< end;j++) Conf[i].mass.push_back(mass[j]/tmass);
		for(int j=cindex[i];j<end;j++) residDip[i]=resid[j];
	}
}
MyDipole::~MyDipole() {
	// TODO Auto-generated destructor stub
}

