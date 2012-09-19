/*
 * LocalEps.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: marchi
 */

#include "LocalEps.h"
void LocalEps::SetUp(){
	try {if(!nr) throw " Residue not initialized ";}
	catch(const char * s){std::cout << s << std::endl;exit(1);}
	if(P) delete [] P;
	if(E) delete [] E;
	if(P_avg) delete [] P_avg;
	if(E_avg) delete [] E_avg;
	if(PM_avg) delete [] PM_avg;
	if(EM_avg) delete [] EM_avg;
	if(eps) delete [] eps;
	if(G_fl) delete [] G_fl;
	if(Gm1_fl) delete [] Gm1_fl;
	if(D_fl) delete [] D_fl;
	if(count) delete [] count;
	P=new Dvect [nr];
	E=new Dvect [nr];
	count=new int [nr];

	P_avg=new Dvect [nr];
	E_avg=new Dvect [nr];

	eps=new Matrix [nnr];
	PM_avg=new Matrix [nr];
	EM_avg=new Matrix [nr];
	M_avg=0.0;
	D_fl=new Matrix [nnr];
	G_fl=new Matrix [nnr];
	Gm1_fl=new Matrix [nnr];
	for(int o=0;o<nr;o++) {
		P_avg[o]=0.0;
		E_avg[o]=0.0;
		PM_avg[o]=0.0;
		EM_avg[o]=0.0;
	}

	for(int o=0;o<nnr;o++) {
		D_fl[o]=0.0;
		G_fl[o]=0.0;
		Gm1_fl[o]=0.0;
	}
	P0_avg=0.0;
	E0_avg=0.0;
	PM0_avg=0.0;
	EM0_avg=0.0;
	D0_fl=0.0;
	G0_fl=0.0;
	Gm10_fl=0.0;
};

void LocalEps::Avg(){
	Matrix * PM_a=new Matrix [nnr];
	Matrix * EM_a=new Matrix [nnr];
	Dvect  * P_a =new Dvect  [nnr];
	Dvect  * E_a =new Dvect  [nnr];
	Dvect  M_a;
	M_a=M_avg/static_cast<double> (TotalCount);
	int totc=0;
	for(int o=0;o<nnr;o++){
		if(resid[ind[o][0]]<0) continue;
		int nop=0;
		for(int n=0;n<ind[o].size();n++){
			int ia=ind[o][n];
			if(count[ia]) {
				P_a[o]+=P_avg[ia]/static_cast<double> (count[ia]);
				E_a[o]+=E_avg[ia]/static_cast<double> (count[ia]);
				PM_a[o]+=PM_avg[ia]/static_cast<double> (count[ia]);
				EM_a[o]+=EM_avg[ia]/static_cast<double> (count[ia]);
				P0_avg+=P_avg[ia]/static_cast<double> (count[ia]);
				E0_avg+=E_avg[ia]/static_cast<double> (count[ia]);
				PM0_avg+=PM_avg[ia]/static_cast<double> (count[ia]);
				EM0_avg+=EM_avg[ia]/static_cast<double> (count[ia]);
				totc++;
				nop++;
			}
		}
		if(nop) P_a[o]/=static_cast<double> (nop);
		if(nop) E_a[o]/=static_cast<double> (nop);
		if(nop) PM_a[o]/=static_cast<double> (nop);
		if(nop) EM_a[o]/=static_cast<double> (nop);
		D_fl[o]=PM_a[o]-P_a[o] % M_a;
		G_fl[o]=EM_a[o]-E_a[o] % M_a;
		std::cout << " Inverting G matrix for residue "<< o << std::endl;
		Gm1_fl[o]=G_fl[o].Inversion();
		Matrix one;

		eps[o]=4.0*M_PI*D_fl[o]*Gm1_fl[o]+one.Unit();
	}
	D0_fl=PM0_avg-P0_avg%M_a;
	G0_fl=EM0_avg-E0_avg%M_a;
	std::cout << " Inverting G matrix for reference molecule " << std::endl;
	Gm10_fl=G0_fl.Inversion();
	Matrix one;

	eps0=4.0*M_PI*D0_fl*Gm10_fl+one.Unit();



	delete [] PM_a;
	delete [] EM_a;
	delete [] P_a;
	delete [] E_a;

};
void LocalEps::Zero(){
	try {if(!nr) throw "Goes above dimensions set";}
	catch(const char * s){std::cout << s << std::endl;exit(1);}
	double yy[DIM]={0,0,0};
	for(int i=0;i<nr;i++){
		P[i]=yy;
		E[i]=yy;
		count[i]=0;
		P_avg[i]=0.0;
		E_avg[i]=0.0;
		PM_avg[i]=0.0;
		EM_avg[i]=0.0;
		eps[i]=0.0;
	}
};
std::ostream & operator<<(std::ostream & ofout, LocalEps & y){
	ofout <<endl;
	for(int i=0;i<y.nnr;i++){
		if(y.Pid(i) < 0) continue;
		for(int p=0;p<DIM;p++){
			ofout << setw(4) << setprecision(0) << fixed<< i << " " << setw(4) << setprecision(0) << p ;
			ofout << fixed << " ---- D = " ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ----  G = " ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ---- Inv(G) =" ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ---- eps =" ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][ZZ] << ' ';
			ofout <<endl;
		}
	}
	for(int p=0;p<DIM;p++){
		ofout << setw(4) << setprecision(0) << fixed<< 0 << " " << setw(4) << setprecision(0) << p ;
		ofout << fixed << " ---- D = " ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][ZZ] << ' ';
		ofout << fixed << " ----  G = " ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][ZZ] << ' ';
		ofout << fixed << " ---- Inv(G) =" ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][ZZ] << ' ';
		ofout << fixed << " ---- eps =" ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][ZZ] << ' ';
		ofout <<endl;
	}
	return ofout;
}
std::ofstream & operator<<(std::ofstream & ofout, LocalEps & y){
	ofout <<endl;
	for(int i=0;i<y.nnr;i++){
		if(y.Pid(i) < 0) continue;
		for(int p=0;p<DIM;p++){
			ofout << setw(4) << setprecision(0) << fixed<< i << " " << setw(4) << setprecision(0) << p ;
			ofout << fixed << " ---- D = " ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.D_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ----  G = " ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.G_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ---- Inv(G) =" ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.Gm1_fl[i][p][ZZ] << ' ';
			ofout << fixed << " ---- eps =" ;
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][XX] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][YY] << ' ';
			ofout << setw(11) << setprecision(4) << fixed<<  y.eps[i][p][ZZ] << ' ';
			ofout <<endl;
		}
	}
	for(int p=0;p<DIM;p++){
		ofout << setw(4) << setprecision(0) << fixed<< 0 << " " << setw(4) << setprecision(0) << p ;
		ofout << fixed << " ---- D = " ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.D0_fl[p][ZZ] << ' ';
		ofout << fixed << " ----  G = " ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.G0_fl[p][ZZ] << ' ';
		ofout << fixed << " ---- Inv(G) =" ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.Gm10_fl[p][ZZ] << ' ';
		ofout << fixed << " ---- eps =" ;
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][XX] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][YY] << ' ';
		ofout << setw(11) << setprecision(4) << fixed<<  y.eps0[p][ZZ] << ' ';
		ofout <<endl;
	}
	return ofout;
}

void LocalEps::Print(){
	cout <<endl;
	for(int i=0;i<nnr;i++){
		if(Pid(i) < 0) continue;
		for(int p=0;p<DIM;p++){
			cout << setw(4) << setprecision(0) << fixed<< i << " " << setw(4) << setprecision(0) << p ;
			cout << fixed << " ---- D = " ;
			cout << setw(11) << setprecision(4) << fixed<<  D_fl[i][p][XX] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  D_fl[i][p][YY] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  D_fl[i][p][ZZ] << ' ';
			cout << fixed << " ----  G = " ;
			cout << setw(11) << setprecision(4) << fixed<<  G_fl[i][p][XX] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  G_fl[i][p][YY] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  G_fl[i][p][ZZ] << ' ';
			cout << fixed << " ---- Inv(G) =" ;
			cout << setw(11) << setprecision(4) << fixed<<  Gm1_fl[i][p][XX] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  Gm1_fl[i][p][YY] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  Gm1_fl[i][p][ZZ] << ' ';
			cout << fixed << " ---- eps =" ;
			cout << setw(11) << setprecision(4) << fixed<<  eps[i][p][XX] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  eps[i][p][YY] << ' ';
			cout << setw(11) << setprecision(4) << fixed<<  eps[i][p][ZZ] << ' ';
			cout <<endl;
		}
//			eps[o][p]=(4.0*M_PI*P_avg[o][p])/E_avg[o][p]+1.0;
//		for(int p=0;p<DIM;p++) std::cout << o << " " << P_avg[o][p] << " " << E_avg[o][p]  << " " << eps[o][p] << std::endl;
	}
};
