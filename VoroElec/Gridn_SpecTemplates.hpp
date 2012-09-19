/*
 * Gridn_SpecTemplates.hpp
 *
 *  Created on: Jun 22, 2011
 *      Author: marchi
 */

#ifndef Gridn_SPECTEMPLATES_HPP_
#define Gridn_SPECTEMPLATES_HPP_
#include "DiffCoeffs_Spec.hpp"
#define ORDER 4
#include "simpson.hpp"


template<>
vector<double> Gridn<1>::Rdf(FILE * fp,const double x[DIM], const string comm, const double & cut, const double & dx){
	double DV=Volume/static_cast<double>(nnx*nny*nnz);
	double xa[DIM];
	int nrdf=(int) (cut/dx)+1;

	double * rdf=new double [nrdf];
	int * irdf=new int [nrdf];
	for(int i=0;i<nrdf;i++) {
		rdf[i]=0.0;
		irdf[i]=0;
	}
	xa[XX] = oc[XX][XX] * x[XX] + oc[XX][YY] * x[YY] + oc[XX][ZZ] * x[ZZ];
	xa[YY] = oc[YY][XX] * x[XX] + oc[YY][YY] * x[YY] + oc[YY][ZZ] * x[ZZ];
	xa[ZZ] = oc[ZZ][XX] * x[XX] + oc[ZZ][YY] * x[YY] + oc[ZZ][ZZ] * x[ZZ];
	cout << xa[XX] << " " << xa[YY] << " " << xa[ZZ] << " " << endl;
	int nx0=(int) nnx;
	int ny0=(int) nny;
	int nz0=(int) nnz;
	double rx,ry,rz;
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
					rdf[h]+=(dist<4.0)?DV*(*this)[i][j][k][0]:0.0;
					irdf[h]++;
				}

			}
		}
	}
	vector<double> gg;
	std::ostringstream oss;
	oss<<"@  s"<<Gridn<DIM>::SetNo << " legend " << "\"" << comm << "\"" << std::endl;
	fprintf(fp,"%s",oss.str().c_str());
	Gridn<DIM>::SetNo++;
	gg.push_back(dx);

	for(int i=0;i<nrdf;i++){
		double bx=(double) (i)*dx;
		if(irdf[i]) {
			double fx=rdf[i]/(double) (irdf[i]);
			fprintf(fp," %f  %e \n",bx,fx);
			gg.push_back(fx);
		}
	}
	fprintf(fp,"& \n");

/*
	printf("&&& \n");
	double dkx=0.5;
	int nkdf=int (256.0/dkx);

	double * fu=new double[nrdf];
	double ff;
	for(int i=0;i<nkdf;i++){
		double kx=(double) (2*M_PI/co[0][0])*int (i)*dkx;
		double sum=0.0;
		for(int j=0;j<nrdf;j++){
			double bx=double (j)*dx;
			double fact;
			if(i==0 || j==0) fact=1.0;
			else fact=sin(bx*kx)/(kx*bx);
			ff=0;
			if(irdf[j]) {
				double fx=rdf[j]/(double) (irdf[j]);
				ff=4.0*M_PI*fx*fact*bx*bx;
			}
			sum+=ff*dx;
		}
		printf("%10.3f  %12.6e \n",0.1*kx,-sum);
	}
	printf("&&& \n");
*/
	return gg;
}
template<>
void Gridn<1>::RdfK(FILE * fp, const string comm, const double & cut, const double & dx){
	printf("dx = %f \n",dx);
	int nrdf=(int) (cut/dx)+1;
	double * rdf=new double [nrdf];
	int * irdf=new int [nrdf];
	for(int i=0;i<nrdf;i++) {
		rdf[i]=0.0;
		irdf[i]=0;
	}
	int nx0=(int) nnx;
	int ny0=(int) nny;
	int nz0=(int) nnz;
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;

	double mw1,mw2,mw3,dist;
	int ia,ja,ka;

	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*0.1*mw1;
				mw2=2.0*M_PI*0.1*mw2;
				mw3=2.0*M_PI*0.1*mw3;
				dist=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
				if(dist < cut) {
					int h=(int) (dist/dx);
					rdf[h]+=(*this)[i][j][k][0];
					irdf[h]++;
				}
			}
		}
	}

	std::ostringstream oss;
	oss<<"@  s"<<Gridn<DIM>::SetNo << " legend " << "\"" << comm << "\"" << std::endl;
	fprintf(fp,"%s",oss.str().c_str());
	Gridn<DIM>::SetNo++;

	for(int i=0;i<nrdf;i++){
		double bx=(double) (i)*dx;
		if(irdf[i]) {
			double fx=rdf[i]/(double) (irdf[i]);
			fprintf(fp," %f  %e\n",bx,fx);
		}
	}
	fprintf(fp,"& \n");
}
template<>
vector<double> Gridn<DIM>::Rdf(FILE * fp,const double x[DIM], const string comm, const double & cut, const double & dx){
	double xa[DIM];
	int nrdf=(int) (cut/dx)+1;

	double * rdf=new double [nrdf];
	int * irdf=new int [nrdf];
	for(int i=0;i<nrdf;i++) {
		rdf[i]=0.0;
		irdf[i]=0;
	}
	xa[XX] = oc[XX][XX] * x[XX] + oc[XX][YY] * x[YY] + oc[XX][ZZ] * x[ZZ];
	xa[YY] = oc[YY][XX] * x[XX] + oc[YY][YY] * x[YY] + oc[YY][ZZ] * x[ZZ];
	xa[ZZ] = oc[ZZ][XX] * x[XX] + oc[ZZ][YY] * x[YY] + oc[ZZ][ZZ] * x[ZZ];
	int nx0=(int) nnx;
	int ny0=(int) nny;
	int nz0=(int) nnz;
	double rx,ry,rz;
	double rnx=1.0/ (double)(nx0);
	double rny=1.0/ (double)(ny0);
	double rnz=1.0/ (double)(nz0);
	double xd[DIM],xc[DIM],dist;
	myvector<double> vr;
	double op;
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
				if(dist != 0.0) {
					vr=Versor(xc);
					if(dist < cut  ){
						int h=(int) (dist/dx);
						op=(*this)[i][j][k][XX]*vr[XX]+(*this)[i][j][k][YY]*vr[YY]+(*this)[i][j][k][ZZ]*vr[ZZ];
						rdf[h]+=op;
						irdf[h]++;
					}
				}

			}
		}
	}
	std::ostringstream oss;
	oss<<"@  s"<<SetNo << " legend " << "\"" << comm << "\"" << std::endl;
	fprintf(fp,"%s",oss.str().c_str());
	SetNo++;
	vector<double> gg;
	gg.push_back(dx);
	for(int i=0;i<nrdf;i++){
		double bx=(double) (i)*dx;
		if(irdf[i]) {
			double fx=rdf[i]/(double) (irdf[i]);
			fprintf(fp," %f  %e %5d \n",bx,fx,i);
			gg.push_back(fx);
		}
		gg.push_back(0.0);
	}
	fprintf(fp,"& \n");
	return gg;
}


template <>
void Gridn<1>::Xplor(FILE * fp, double convf){
	int nx_s,ny_s,nz_s;
	int nx_e,ny_e,nz_e;
	double convf2=convf*convf;
	double avg=0.0, avg2=0.0;
	double a=0.0,b=0.0,c=0.0,alpha=0.0, beta=0, gamma=0;
	Gridn<1> & MyGridn=*this;

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	for(int i=0; i< nx0; i++){
		for(int j=0; j< ny0; j++){
			for(int k=0; k< nz0; k++){

				avg+=MyGridn[i][j][k][0]*convf;
				avg2+=MyGridn[i][j][k][0]*MyGridn[i][j][k][0]*convf2;
			}

		}
	}
	avg/=static_cast<double>(nx0*ny0*nz0);
	avg2/=static_cast<double>(nx0*ny0*nz0);

	nx_s=0;
	ny_s=0;
	nz_s=0;
	nx_e=nx0-1+nx_s;
	ny_e=ny0-1+ny_s;
	nz_e=nz0-1+nz_s;


	alpha=Alpha()*RAD2DEG;
	beta=Beta()*RAD2DEG;
	gamma=Gamma()*RAD2DEG;

	a=norm(co[XX]);
	b=norm(co[YY]);
	c=norm(co[ZZ]);
	fprintf(fp,"\n       2 !NTITLE\n") ;
	fprintf(fp," REMARKS Electrostatic Potential from GROMACS\n") ;
	fprintf(fp," REMARKS DATE: 2010-04-21\n") ;
	fprintf(fp,"%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",nx,nx_s,nx_e,ny,ny_s,ny_e,nz,nz_s,nz_e);
	fprintf(fp,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",10*a,10*b,10*c,alpha,beta,gamma);
	fprintf(fp, "ZYX\n") ;
	int ja;
	for(int k=0;k<nz0;k++){
		fprintf(fp,"%8d\n",k);
		ja=0;
		for(int j=0;j<ny0;j++){
			for(int i=0;i<nx0;i++){
				ja++;
				fprintf(fp,"%12.5E",MyGridn[i][j][k][0]*convf);
				if(!(ja%6)) fprintf(fp,"\n");
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"%8d\n",-9999);

	double stddev=sqrt(avg2-avg*avg);
	fprintf(fp,"%12.4E %12.4E\n",avg,stddev);

}

template<>
Gridn<DIM>  Gridn<1>::Der1(){
	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);

	array3<double> fi(nnx,nny,nnz,align);
	array3<Complex> gi(nnx,nny,nzp,align);
	array3<Complex> dgi[DIM];
	for(int i=0;i<DIM;i++) dgi[i].Allocate(nnx,nny,nzp,align);

	array4<double> dfi0(DIM,nnx,nny,nnz);
	array4<double> dfi(nnx,nny,nnz,DIM);

	rcfft3d Forward3(nnz,fi,gi);
	crfft3d Backward3(nnz,gi,fi);

	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;


	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++){
				fi[i][j][k]=(*this)[i][j][k][0];
			}


	Forward3.fft(fi,gi);

	int ia,ja,ka;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	Complex imag(0.0,1.0), fact, zero(0,0);
	double mysinx, mysiny,mysinz;
	DiffCoeffs<ORDER> test(nx0,ny0,nz0,co[0][0],co[1][1],co[2][2]);

	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		mysinx=test.coeffs(XX,ia);
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			mysiny=test.coeffs(YY,ja);
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mysinz=test.coeffs(ZZ,ka);
				fact=imag*gi[i][j][k];
				dgi[XX][i][j][k]=mysinx*fact;
				dgi[YY][i][j][k]=mysiny*fact;
				dgi[ZZ][i][j][k]=mysinz*fact;
			}
		}
	}
	Backward3.fftNormalized(dgi[XX],dfi0[XX]);
	Backward3.fftNormalized(dgi[YY],dfi0[YY]);
	Backward3.fftNormalized(dgi[ZZ],dfi0[ZZ]);
	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++){
				dfi[i][j][k][XX] = dfi0[XX][i][j][k];
				dfi[i][j][k][YY] = dfi0[YY][i][j][k];
				dfi[i][j][k][ZZ] = dfi0[ZZ][i][j][k];
			}

 	return dfi;
}
template<>
Gridn<DIM>  Gridn<1>::Der2(){

	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);

	array3<double> fi(nnx,nny,nnz,align);
	array3<Complex> gi(nnx,nny,nzp,align);
	array3<Complex> gi0(nnx,nny,nnz,align);
	array3<Complex> dgi[DIM];
	for(int i=0;i<DIM;i++) dgi[i].Allocate(nnx,nny,nzp,align);

	array4<double> dfi0(DIM,nnx,nny,nnz);
	array4<double> dfi(nnx,nny,nnz,DIM);

	rcfft3d Forward3(nnz,fi,gi);
	crfft3d Backward3(nnz,gi,fi);

	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;


	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++){
				fi[i][j][k]=(*this)[0][i][j][k];
			}


	Forward3.fft(fi,gi);

	int ia,ja,ka;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	double mw1,mw2,mw3;
	DiffCoeffs<2> test(nx0,ny0,nz0,co[0][0],co[1][1],co[2][2]);

	Complex imag(0.0,1.0), fact, zero(0,0);

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
				fact=imag*gi[i][j][k];
				dgi[XX][i][j][k]=mw1*fact;
				dgi[YY][i][j][k]=mw2*fact;
				dgi[ZZ][i][j][k]=mw3*fact;
			}
		}
	}
/*
	for(int i=0;i<nx0;i++){
		for(int j=0;j<ny0;j++){
			int k=32;
			dgi[ZZ][i][j][k]=zero;
		}
	}

	{
		int i=32;int k=32;
		for(int j=0;j<ny0;j++){
			dgi[XX][i][j][k]=zero;
		}
		i=32;k=0;
		for(int j=0;j<ny0;j++){
			dgi[XX][i][j][k]=zero;
		}
	}
	{
		int j=32;int k=32;
		for(int i=0;i<nx0;i++){
			dgi[YY][i][j][k]=zero;
		}
		j=32;k=0;
		for(int i=0;i<nx0;i++){
			dgi[YY][i][j][k]=zero;
		}
	}
*/
	Backward3.fftNormalized(dgi[XX],dfi0[XX]);
	Backward3.fftNormalized(dgi[YY],dfi0[YY]);
	Backward3.fftNormalized(dgi[ZZ],dfi0[ZZ]);
	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++){
				dfi[i][j][k][XX] = dfi0[XX][i][j][k];
				dfi[i][j][k][YY] = dfi0[YY][i][j][k];
				dfi[i][j][k][ZZ] = dfi0[ZZ][i][j][k];
			}

 	return dfi;
}

template<>
Gridn<1>  Gridn<DIM>::Div(){
	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);

	array3<Complex> rhok(nnx,nny,nzp,align);
	array3<Complex> dipg[DIM];
	array3<double> rhor(nnx,nny,nnz,align);
	array3<double> dip[DIM];
	for(int i=0;i<DIM;i++) {
		dip[i].Allocate(nnx,nny,nnz,align);
		dipg[i].Allocate(nnx,nny,nzp,align);
	}
	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++){
				dip[XX][i][j][k]=(*this)[XX][i][j][k];
				dip[YY][i][j][k]=(*this)[YY][i][j][k];
				dip[ZZ][i][j][k]=(*this)[ZZ][i][j][k];
			}

	rcfft3d Forward3(nnz,dip[XX],dipg[XX]);
	crfft3d Backward3(nnz,rhok,rhor);


	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;


	Forward3.fft(dip[XX],dipg[XX]);
	Forward3.fft(dip[YY],dipg[YY]);
	Forward3.fft(dip[ZZ],dipg[ZZ]);

	int ia,ja,ka;
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	Complex imag(0.0,1.0), fact, zero(0,0), px,py,pz;
	DiffCoeffs<ORDER> test(nx0,ny0,nz0,co[0][0],co[1][1],co[2][2]);
	double mysinx, mysiny,mysinz;
	for(int i=0;i<nx0;i++){
		ia=(i<nfx)?i : i-nx0;
		mysinx=test.coeffs(XX,ia);
		for(int j=0;j<ny0;j++){
			ja=(j<nfy)?j : j-ny0;
			mysiny=test.coeffs(YY,ja);
			for(int k=0;k<nz0/2+1;k++){
				ka=(k<nfz)?k : k-nz0;
				mysinz=test.coeffs(ZZ,ka);
				px=mysinx*imag*dipg[XX][i][j][k];
				py=mysiny*imag*dipg[YY][i][j][k];
				pz=mysinz*imag*dipg[ZZ][i][j][k];
				rhok[i][j][k]=px+py+pz;
			}
		}
	}
	Backward3.fftNormalized(rhok,rhor);
	Gridn<1> temp;
	temp.Allocate();
	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++)
				temp[i][j][k][0]=rhor[i][j][k];
 	return temp;
}
template<>
	Gridn<1> Gridn<1>::Derivative(unsigned int Dir){
	Gridn<1> temp;
	return temp;
}


#endif /* Gridn_SPECTEMPLATES_HPP_ */
