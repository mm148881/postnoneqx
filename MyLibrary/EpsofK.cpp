/*
 * EpsofK.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: marchi
 */

#include "EpsofK.h"


EpsofK::~EpsofK() {
	// TODO Auto-generated destructor stub
}
void EpsofK::SofK(Rho & ro){
	Array3<double> ro_r;
	Array3<Complex> ro_k;
	unsigned int nzp=nnz/2+1;
	try {
		if(ro.getnnx() != nnx || ro.getnny() != nny || ro.getnnz() != nnz ) throw "EpsofK: Dimensions of Rho do not match";
	}
	catch(const char * s){
		cout << s << endl;
		exit(1);
	}

	ro_r.Allocate(nnx,nny,nnz);
	ro_k.Allocate(nnx,nny,nzp);

	ro_r=ro;

	rcfft3d Forward3(nnz,ro_r,ro_k);

	Forward3.fft(ro_r,ro_k);

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	double factor=4.0*M_PI*(EFACT/KT)*sqr(EpsofK::getDV())/Volume;


	int ia,ja,ka;
	double mw1,mw2,mw3,mw;
	Complex zero(0.0,0.0);
	Complex vt,v0;
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;

	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*mw1;
				mw2=2.0*M_PI*mw2;
				mw3=2.0*M_PI*mw3;
				mw=mw1*mw1+mw2*mw2+mw3*mw3;
				vt=ro_k[i][j][k];
				v0=vt*conj(vt);
				if(mw != 0.0) (*this)[0][i][j][k]=factor*v0.real()/mw;
				else (*this)[0][i][j][k]=0.0;
			}
		}
	}

/*
	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*mw1;
				mw2=2.0*M_PI*mw2;
				mw3=2.0*M_PI*mw3;
				mw=mw1*mw1+mw2*mw2+mw3*mw3;
				vt=ro_k[i][j][k];
				v0=vt;
				if(mw != 0.0) (*this)[0][i][j][k]=v0.real();
				else (*this)[0][i][j][k]=0.0;
			}
		}
	}
*/
}

void EpsofK::SofK(RhoDip & dip){
	Array3<double> dip_r[DIM];
	Array3<Complex> dip_k[DIM];
	unsigned int nzp=nnz/2+1;
	double factor=4.0*M_PI*(EFACT/KT)*sqr(EpsofK::getDV())/Volume;
	try {
		if(dip.getnnx() != nnx || dip.getnny() != nny || dip.getnnz() != nnz ) throw "EpsofK: Dimensions of Rho do not match";
	}
	catch(const char * s){
		cout << s << endl;
		exit(1);
	}

	for(int m=0;m<DIM;m++){
		dip_r[m].Allocate(nnx,nny,nnz);
		dip_k[m].Allocate(nnx,nny,nzp);
		dip_r[m]=dip[m];
	}

	rcfft3d Forward3(nnz,dip_r[XX],dip_k[XX]);

	Forward3.fft(dip_r[XX],dip_k[XX]);
	Forward3.fft(dip_r[YY],dip_k[YY]);
	Forward3.fft(dip_r[ZZ],dip_k[ZZ]);

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;

	Complex vx,vy,vz,vt,v0;
	double mw1,mw2,mw3,mw;
	Complex zero(0.0,0.0);
	for(int i=0;i<nx0;i++){
		int ia=(i<nfx)?i : i-nx0;
		for(int j=0;j<ny0;j++){
			int ja=(j<nfy)?j : j-ny0;
			for(int k=0;k<nz0/2+1;k++){
				int ka=(k<nfz)?k : k-nz0;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw=mw1*mw1+mw2*mw2+mw3*mw3;
				vx=dip_k[XX][i][j][k];
				vy=dip_k[YY][i][j][k];
				vz=dip_k[ZZ][i][j][k];
				vt=mw1*vx+mw2*vy+mw3*vz;
				v0=vt*conj(vt);
				if(mw != zero) (*this)[0][i][j][k]=factor*v0.real()/mw;
				else (*this)[0][i][j][k]=0.0;
			}
		}
	}
}
void EpsofK::SofK(Rho & x_0, Rho & x_1){ // x_0=Rho_0 x_1 = Rho_Sol
	Array3<double> ro_r0;
	Array3<double> ro_r1;
	Array3<Complex> ro_k0;
	Array3<Complex> ro_k1;
	unsigned int nzp=nnz/2+1;
	try {
		if(x_0.getnnx() != nnx || x_0.getnny() != nny || x_0.getnnz() != nnz ) throw "EpsofK: Dimensions of Rho do not match";
	}
	catch(const char * s){
		cout << s << endl;
		exit(1);
	}
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	ro_r0.Allocate(nnx,nny,nnz);
	ro_r1.Allocate(nnx,nny,nnz);
	ro_k0.Allocate(nnx,nny,nzp);
	ro_k1.Allocate(nnx,nny,nzp);
	ro_r0=x_0;
	ro_r1=x_1;



	rcfft3d Forward3D(nnz,ro_r0,ro_k0);

	Forward3D.fft(ro_r0,ro_k0);
	Forward3D.fft(ro_r1,ro_k1);



	Complex zero(0.0,0.0);
	Complex v_0,v_1,v10;
	double v11;


	for(int i=0;i<nx0;i++){
		for(int j=0;j<ny0;j++){
			for(int k=0;k<nz0/2+1;k++){
				v_0=ro_k0[i][j][k];
				v_1=ro_k1[i][j][k];
				if(v_0 != zero) v10=v_1/v_0;
				else v10=zero;
				(*this)[0][i][j][k]=-v10.real();
			}
		}
	}
/*
	for(int i=0;i<nx0;i++){
		int ib=(nx0-i)%nx0;
		for(int j=0;j<ny0;j++){
			int jb=(ny0-j)%ny0;
			for(int k=0;k<nz0/2+1;k++){
				int kb=(nz0-k)%nz0;
				(*this)[0][ib][jb][kb]=(*this)[0][i][j][k];
			}
		}
	}

*/
}
void EpsofK::SofK(Phi & x_0, Phi & x_1){ // x_0=Phi_0 x_1 = Phi_Sol
	Array3<double> ro_r0;
	Array3<double> ro_r1;
	Array3<Complex> ro_k0;
	Array3<Complex> ro_k1;
	unsigned int nzp=nnz/2+1;
	try {
		if(x_0.getnnx() != nnx || x_0.getnny() != nny || x_0.getnnz() != nnz ) throw "EpsofK: Dimensions of Rho do not match";
	}
	catch(const char * s){
		cout << s << endl;
		exit(1);
	}
	ro_r0.Allocate(nnx,nny,nnz);
	ro_r1.Allocate(nnx,nny,nnz);
	ro_k0.Allocate(nnx,nny,nzp);
	ro_k1.Allocate(nnx,nny,nzp);

	ro_r0=x_0;
	ro_r1=x_1;
	rcfft3d Forward3D(nnz,ro_r0,ro_k0);

	Forward3D.fft(ro_r0,ro_k0);
	Forward3D.fft(ro_r1,ro_k1);

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);


	Complex zero(0.0,0.0);
	Complex v_0,v_1,v10;


	for(int i=0;i<nx0;i++){
		for(int j=0;j<ny0;j++){
			for(int k=0;k<nz0/2+1;k++){
				v_0=ro_k0[i][j][k];
				v_1=ro_k1[i][j][k];
				if(v_0 != zero) v10=v_1/v_0;
				else v10=zero;
				(*this)[0][i][j][k]=v10.real();
			}
		}
	}
}
