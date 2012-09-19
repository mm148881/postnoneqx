/*
 * Rho.cpp
 *
 *  Created on: May 24, 2011
 *      Author: marchi
 */

#include "Rho.h"
#include <cmath>
using std::ofstream;
using std::ifstream;
using std::ios_base;

Rho & Rho::operator=(const Grid<1> & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}
Rho & Rho::operator=(const double & x){
	Grid<1> & My=*this;
	My=x;
	return *this;
}
void Rho::Density(Charges & x){
// int get_charges(real*** Rho_slt, real*** Rho_sol,real* Charges,int nx,int ny,int nz,int nind, atom_id *cindex, t_trxframe* fr){

	try{
		if(!x.test()) throw "Cannot compute density if charges not initialized.";
	}
	catch(const char *s){
		std::cout << s << std::endl;
	}
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	int natoms=x.getNR();
	array2<real> xa(natoms,DIM);

	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int mx,my,mz;
	int ox,oy,oz,tx,ty,tz;
	Metric mtt=x.getMt();
	const matrix & co1=mtt.getCO();
	setMetric(co1);
	double DV=Rho::getDV();
	if(strstr(typeid(*this).name(),"RhoAverage") == NULL) *this=0.0;
	else this->Accumulate();

	for (int i = 0; i < natoms; i++) {
		xa[i][XX] = oc[XX][XX] * x[i][XX] + oc[XX][YY] * x[i][YY] + oc[XX][ZZ] * x[i][ZZ];
		xa[i][YY] = oc[YY][XX] * x[i][XX] + oc[YY][YY] * x[i][YY] + oc[YY][ZZ] * x[i][ZZ];
		xa[i][ZZ] = oc[ZZ][XX] * x[i][XX] + oc[ZZ][YY] * x[i][YY] + oc[ZZ][ZZ] * x[i][ZZ];
	}
	for (int i = 0; i < natoms; i++) {
		x1 = xa[i][XX];
		y1 = xa[i][YY];
		z1 = xa[i][ZZ];
		r1 = static_cast<real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (real) mx;
		gy = s1 - (real) my;
		gz = t1 - (real) mz;
		real chg[1];
		chg[0]=x.getQ(i);
		TriLinear(gx, gy, gz, chg);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			(*this)[0][ox][oy][oz] += static_cast<double>(qq[n][0])/DV;
		}
	}
	xa.Deallocate();
}
void Rho::Density1(Charges & x){
	try{
		if(!x.test()) throw "Cannot compute density if charges not initialized.";
	}
	catch(const char *s){
		std::cout << s << std::endl;
	}
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	int natoms=x.getNR();
	array2<real> xa(natoms,DIM);

	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int mx,my,mz;
	int ox,oy,oz,tx,ty,tz,oTot;
	Metric mtt=x.getMt();
	const matrix & co1=mtt.getCO();
	setMetric(co1);
	double DV=Rho::getDV();
	if(strstr(typeid(*this).name(),"RhoAverage") == NULL) *this=0.0;
	else this->Accumulate();

	multiset<DataDensity,DensityComp> mydens;
	DataDensity Val;
	for (int i = 0; i < natoms; i++) {
		xa[i][XX] = oc[XX][XX] * x[i][XX] + oc[XX][YY] * x[i][YY] + oc[XX][ZZ] * x[i][ZZ];
		xa[i][YY] = oc[YY][XX] * x[i][XX] + oc[YY][YY] * x[i][YY] + oc[YY][ZZ] * x[i][ZZ];
		xa[i][ZZ] = oc[ZZ][XX] * x[i][XX] + oc[ZZ][YY] * x[i][YY] + oc[ZZ][ZZ] * x[i][ZZ];
	}
	for (int i = 0; i < natoms; i++) {
		x1 = xa[i][XX];
		y1 = xa[i][YY];
		z1 = xa[i][ZZ];
		r1 = static_cast<real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (real) mx;
		gy = s1 - (real) my;
		gz = t1 - (real) mz;
		real chg[1];
		chg[0]=x.getQ(i);
		TriLinear(gx, gy, gz, chg);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			Val.n=oz+oy*nz0+ox*nz0*ny0;
			Val.value=static_cast<double>(qq[n][0])/DV;
			mydens.insert(Val);
		}
	}
	double * pt=&(*this)[0][0][0][0];
	for(multiset<DataDensity>::iterator it=mydens.begin();it!=mydens.end();++it){
		int n=it->n;
		pt[n]+=it->value;
	}
	xa.Deallocate();
}
void Rho::Linear(Rho Chi, Rho ro0){
	unsigned int nzp=nnz/2+1;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	Array3<double> ro_r,rosol_r;
	Array3<Complex> ro_k,rosol_k;

	ro_r.Allocate(nnx,nny,nnz);
	ro_k.Allocate(nnx,nny,nzp);
	rosol_k.Allocate(nnx,nny,nzp);
	rosol_r.Allocate(nnx,nny,nnz);

	(*this)[0].Allocate(nnx,nny,nnz);
	ro_r=ro0;
	rcfft3d Forward3(nnz,ro_r,ro_k);
	crfft3d Backward3(nnz,rosol_k,rosol_r);

	Forward3.fft(ro_r,ro_k);


	for(int i=0;i<nx0;i++)
		for(int j=0;j<ny0;j++)
			for(int k=0;k<nz0/2+1;k++){
				rosol_k[i][j][k].real()=-Chi[0][i][j][k]*ro_k[i][j][k].real();
				rosol_k[i][j][k].imag()=0.0;
			}

	Backward3.fftNormalized(rosol_k,rosol_r);
	(*this)[0]=rosol_r;
}
void Rho::Symmetrify(rvec & x){

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	double rx,ry,rz;
	rvec xa;
	real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int mx,my,mz;
	int ox,oy,oz,tx,ty,tz;
	double DV=Rho::getDV();
	int nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	int nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	int nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;
	xa[XX] = oc[XX][XX] * x[XX] + oc[XX][YY] * x[YY] + oc[XX][ZZ] * x[ZZ];
	xa[YY] = oc[YY][XX] * x[XX] + oc[YY][YY] * x[YY] + oc[YY][ZZ] * x[ZZ];
	xa[ZZ] = oc[ZZ][XX] * x[XX] + oc[ZZ][YY] * x[YY] + oc[ZZ][ZZ] * x[ZZ];
	double dx=0.02;
	double cut=0.9*co[0][0];
	int nrdf=int (cut/dx);

	double * val=new double[nrdf];
	int * nval=new int [nrdf];
	for(int i=0;i<nrdf;i++)  {
		*(val+i)=0.0;
		*(nval+i)=0;
	}

	double rnx=1.0/ (double)(nx0);
	double rny=1.0/ (double)(ny0);
	double rnz=1.0/ (double)(nz0);
	double xd[DIM],xc[DIM],dist;
	for(int i=0;i<nx0;i++){
		rx=(double) (i)*rnx;
		for(int j=0;j<ny0;j++){
			ry=(double) (j)*rny;
			for(int k=0;k<nz0;k++){
				rz=(double) (k)*rnz;
				xd[XX]=rx-xa[XX];
				xd[YY]=ry-xa[YY];
				xd[ZZ]=rz-xa[ZZ];
				xd[XX]=xd[XX]-rint(xd[XX]);
				xd[YY]=xd[YY]-rint(xd[YY]);
				xd[ZZ]=xd[ZZ]-rint(xd[ZZ]);

				xc[XX]=co[XX][XX]*xd[XX]+co[XX][YY]*xd[YY]+co[XX][ZZ]*xd[ZZ];
				xc[YY]=co[YY][XX]*xd[XX]+co[YY][YY]*xd[YY]+co[YY][ZZ]*xd[ZZ];
				xc[ZZ]=co[ZZ][XX]*xd[XX]+co[ZZ][YY]*xd[YY]+co[ZZ][ZZ]*xd[ZZ];
				dist=sqrt(xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ]);
				if(dist < cut){
					int h=(int) (dist/dx);
					val[h]+=(*this)[0][i][j][k];
					nval[h]++;
				}
			}
		}
	}
	for(int i=0;i<nrdf;i++) {
		val[i]=(nval[i])?val[i]/double (nval[i]):0.0;
	}

	for(int i=0;i<nx0;i++){
		rx=(double) (i)*rnx;
		for(int j=0;j<ny0;j++){
			ry=(double) (j)*rny;
			for(int k=0;k<nz0;k++){
				rz=(double) (k)*rnz;
				xd[XX]=rx-xa[XX];
				xd[YY]=ry-xa[YY];
				xd[ZZ]=rz-xa[ZZ];
				xd[XX]=xd[XX]-rint(xd[XX]);
				xd[YY]=xd[YY]-rint(xd[YY]);
				xd[ZZ]=xd[ZZ]-rint(xd[ZZ]);

				xc[XX]=co[XX][XX]*xd[XX]+co[XX][YY]*xd[YY]+co[XX][ZZ]*xd[ZZ];
				xc[YY]=co[YY][XX]*xd[XX]+co[YY][YY]*xd[YY]+co[YY][ZZ]*xd[ZZ];
				xc[ZZ]=co[ZZ][XX]*xd[XX]+co[ZZ][YY]*xd[YY]+co[ZZ][ZZ]*xd[ZZ];
				dist=sqrt(xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ]);
				int h=(int) (dist/dx);
				(*this)[0][i][j][k]=val[h];
			}
		}
	}
}

std::ofstream & operator<<(std::ofstream & fout, Rho & y){
	fout.write((char *) &y.nnx, sizeof y.nnx);
	fout.write((char *) &y.nny, sizeof y.nny);
	fout.write((char *) &y.nnz, sizeof y.nnz);
	fout.write((char *) &y.co, sizeof y.co);
	fout.write((char *) &y.co, sizeof y.oc);
	fout.write((char *) &y[0][0][0][0], (sizeof y[0][0][0][0])*y.nnx*y.nny*y.nnz);
	return fout;
}
std::ifstream & operator>>(std::ifstream & fin, Rho & y){
	unsigned int nnx,nny,nnz;
	matrix co,oc;
	fin.read((char *) &nnx, sizeof nnx);
	fin.read((char *) &nny, sizeof nny);
	fin.read((char *) &nnz, sizeof nnz);
	fin.read((char *) &co, sizeof co);
	fin.read((char *) &oc, sizeof oc);

	y.set(co,nnx,nny,nnz);
	y.Allocate();
	fin.read((char *) &y[0][0][0][0], (sizeof y[0][0][0][0])*nnx*nny*nnz);
	return fin;
}
