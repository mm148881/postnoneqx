/*
 * Polarization.cpp
 *
 *  Created on: Dec 27, 2011
 *      Author: marchi
 */

#include "Polarization.h"

namespace PolarizationNS {

Polarization & Polarization::operator=(const double x){
	Gridn<DIM> & My=*this;
	My=x;
	return *this;
}
void MyDensity(Rho & r,Polarization & d,MyDipole & x){
	try{
		if(!x.test()) throw "Cannot compute density if Dipoles not initialized.";
	}
	catch(const char *s){
		std::cout << s << std::endl;
	}
	int nx0=static_cast<int>(Gridn<DIM>::getnnx());
	int ny0=static_cast<int>(Gridn<DIM>::getnny());
	int nz0=static_cast<int>(Gridn<DIM>::getnnz());

	int ndip=x.getNdip();

	double xa,ya,za,xc,yc,zc,xb,yb,zb;
	Dvect c;
	int natoms=x.getNatom();
	d.SetDip()=0.0;
	double tchg=0.0;
	for(int i=0;i<natoms;i++){
		tchg+=x.getQ(i);
	}

	tchg/=static_cast<double> (natoms);
	for(int i=0;i<natoms;i++){
		c[XX]=x.Atoms::getX()[i][XX];
		c[YY]=x.Atoms::getX()[i][YY];
		c[ZZ]=x.Atoms::getX()[i][ZZ];
		double fact=(x.getQ(i)-tchg);
		d.SetDip()[XX]+=c[XX]*fact;
		d.SetDip()[YY]+=c[YY]*fact;
		d.SetDip()[ZZ]+=c[ZZ]*fact;
	}
	int mx,my,mz;
	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int ox,oy,oz,tx,ty,tz;
	double DV=Polarization::getDV();
	const rvec * cube=&Gridn<DIM>::Getcube()[0];

	d=0.0;
	for (int i = 0; i < ndip; i++) {
		x1 = x.getXa(i)[XX];
		y1 = x.getXa(i)[YY];
		z1 = x.getXa(i)[ZZ];
		r1 = static_cast<real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (real) mx;
		gy = s1 - (real) my;
		gz = t1 - (real) mz;
		rvec chg;
		float qchg=x.getQ(i);
		chg[XX]=x.getDip(i)[XX];
		chg[YY]=x.getDip(i)[YY];
		chg[ZZ]=x.getDip(i)[ZZ];
		vector<float> gg=MyTriLinear(gx, gy, gz);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			r[0][ox][oy][oz]+= static_cast<double>(gg[n]*qchg)/DV;
			for(int o=0;o<DIM; o++)
				d[ox][oy][oz][o] += static_cast<double>(gg[n]*chg[o])/DV;
		}
	}
}
void Polarization::Density(MyDipole & x){
	try{
		if(!x.test()) throw "Cannot compute density if Dipoles not initialized.";
	}
	catch(const char *s){
		std::cout << s << std::endl;
	}
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	int ndip=x.getNdip();

	double xa,ya,za,xc,yc,zc,xb,yb,zb;
	Dvect c;
	int natoms=x.getNatom();
	Mdip=0.0;
	double tchg=0.0;
	for(int i=0;i<natoms;i++){
		tchg+=x.getQ(i);
	}

	tchg/=static_cast<double> (natoms);
	for(int i=0;i<natoms;i++){
		c[XX]=x.Atoms::getX()[i][XX];
		c[YY]=x.Atoms::getX()[i][YY];
		c[ZZ]=x.Atoms::getX()[i][ZZ];
		double fact=(x.getQ(i)-tchg);
		Mdip[XX]+=c[XX]*fact;
		Mdip[YY]+=c[YY]*fact;
		Mdip[ZZ]+=c[ZZ]*fact;
	}
	int mx,my,mz;
	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int ox,oy,oz,tx,ty,tz;
	double DV=Polarization::getDV();
	*this=0.0;
	rvec summa={0.0,0.0,0.0};
	for (int i = 0; i < ndip; i++) {
		x1 = x.getXa(i)[XX];
		y1 = x.getXa(i)[YY];
		z1 = x.getXa(i)[ZZ];
		r1 = static_cast<real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (real) mx;
		gy = s1 - (real) my;
		gz = t1 - (real) mz;
		rvec chg;
		chg[XX]=x.getDip(i)[XX];
		chg[YY]=x.getDip(i)[YY];
		chg[ZZ]=x.getDip(i)[ZZ];
		summa[XX]+=chg[XX];
		summa[YY]+=chg[YY];
		summa[ZZ]+=chg[ZZ];
		TriLinear(gx, gy, gz, chg);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			for(int o=0;o<DIM; o++)
				(*this)[ox][oy][oz][o] += static_cast<double>(qq[n][o])/DV;
		}
	}
}

Polarization::~Polarization() {
	// TODO Auto-generated destructor stub
}

} /* namespace PolarizationNS */
