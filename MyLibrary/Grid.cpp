/*
 * Grid.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: marchi
 */
#ifdef GRID_H_
#include "Grid.h"

template <unsigned int Ndim> matrix Grid<Ndim>::co={{0,0,0},{0,0,0},{0,0,0}};
template <unsigned int Ndim> matrix Grid<Ndim>::oc={{0,0,0},{0,0,0},{0,0,0}};
template <unsigned int Ndim> double Grid<Ndim>::Volume=0.0;
template <unsigned int Ndim> unsigned int Grid<Ndim>::nnx=0;
template <unsigned int Ndim> unsigned int Grid<Ndim>::nny=0;
template <unsigned int Ndim> unsigned int Grid<Ndim>::nnz=0;
template <unsigned int Ndim> const rvec Grid<Ndim>::cube[VERTEX]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
template <unsigned int Ndim> array3<Complex> * Grid<Ndim>::filter=NULL;
template <unsigned int Ndim> array3<Complex> * Grid<Ndim>::filterrc=NULL;
template <unsigned int Ndim> int Grid<Ndim>::SetNo=0;
static int Fcounter=0;

template <class T>
T sqr(T a){return a*a;}

template <unsigned int Ndim>
Grid<Ndim>::Grid()
{
		try{
			if(nnx+nny+nnz==0) throw "Warning: Metrics initialization is postponed ";
		}
		catch(const char * s){
			std::cout << s << std::endl;
		}
		Allocate();
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(Grid<Ndim> & y) {
	Allocate();
	*this=y;
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(const Grid<Ndim> & y){
	Allocate();
	*this=y;
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(const double & x){
	this->Allocate();
	*this=x;
}

template <unsigned int Ndim>
Grid<Ndim> & Grid<Ndim>::operator=(const double y){
	double * pt=&(*this)[0][0][0][0];
	for(int i=0;i<static_cast<int>(Ndim*nnx*nny*nnz);i++) pt[i]=y;
	return *this;
}

template <unsigned int Ndim>
Grid<Ndim> & Grid<Ndim>::operator+=(const Grid<Ndim> & y){
	double * pt=&(*this)[0][0][0][0];
	double * py=&y[0][0][0][0];
	for(int i=0;i<static_cast<int>(Ndim*nnx*nny*nnz);i++) pt[i]+=py[i];
	return *this;
}

template <unsigned int Ndim>
Grid<Ndim> Grid<Ndim>::operator-(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	temp-=a;
	return temp;
}

template <unsigned int Ndim>
Grid<Ndim> Grid<Ndim>::operator+(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	temp+=a;
	return temp;
}

template <unsigned int Ndim>
Grid<Ndim> Grid<Ndim>::operator*(const double & a){
	Grid<Ndim> temp=*this;
	temp*=a;
	return temp;
}

template <unsigned int Ndim>
Grid<Ndim> Grid<Ndim>::operator/(const double & a){
	Grid<Ndim> temp=*this;
	temp/=a;
	return temp;
}
template <unsigned int Ndim>
Grid<Ndim> Grid<Ndim>::operator/(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					if(a[l][i][j][k] == 0) temp[l][i][j][k]=0.0;
					else temp[l][i][j][k]/=a[l][i][j][k];
				}
	return temp;
}

template <unsigned int Ndim>
void Grid<Ndim>::setMetric(const matrix & co_in)
{
  /* Save some time by assuming lower right part is zero */
  for(int i=0;i<DIM;i++)
	  for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];

  float tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
  oc[XX][XX]=co[YY][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][XX]=0;
  oc[ZZ][XX]=0;
  oc[XX][YY]=-co[XX][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][YY]=co[XX][XX]*co[ZZ][ZZ]*tmp;
  oc[ZZ][YY]=0;
  oc[XX][ZZ]=(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ])*tmp;
  oc[YY][ZZ]=-co[YY][ZZ]*co[XX][XX]*tmp;
  oc[ZZ][ZZ]=co[XX][XX]*co[YY][YY]*tmp;
  Volume=(co[XX][XX]*(co[YY][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[YY][ZZ])
		  -co[YY][XX]*(co[XX][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[XX][ZZ])
		  +co[ZZ][XX]*(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ]));
}


template <unsigned int Ndim>
double Grid<Ndim>::Alpha(){
	double alpha;
	if (norm2(co[YY])*norm2(co[ZZ])!=0)
		alpha = acos(cos_angle(co[YY],co[ZZ]));
	else
		alpha = 0.5*M_PI;

	return alpha;
}

template <unsigned int Ndim>
double Grid<Ndim>::Beta(){
	double beta;
	if (norm2(co[XX])*norm2(co[ZZ])!=0)
		beta  = acos(cos_angle(co[XX],co[ZZ]));
	else
		beta  = 0.5*M_PI;
	return beta;
}

template <unsigned int Ndim>
double Grid<Ndim>::Gamma(){
	double gamma;
	if (norm2(co[XX])*norm2(co[YY])!=0)
		gamma = acos(cos_angle(co[XX],co[YY]));
	else
		gamma = 0.5*M_PI;
	return gamma;
}

template <unsigned int Ndim>
Grid<Ndim>::~Grid() {
	// TODO Auto-generated destructor stub
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(array4<Complex> & x){
	Allocate();
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					(*this)[l][i][j][k]=x[l][i][j][k].real();
				}
}
template <unsigned int Ndim>
Grid<Ndim> & Grid<Ndim>::operator=(const array4<Complex> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<nnx;i++)
			for(unsigned j=0;j<nny;j++)
				for(unsigned k=0;k<nnz;k++)
					(*this)[l][i][j][k]=a[l][i][j][k].real();
				}
	return *this;

}
template <unsigned int Ndim>
Grid<Ndim>::Grid(array4<double> & x){
	Allocate();
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					(*this)[l][i][j][k]=x[l][i][j][k];
				}
}
template <unsigned int Ndim>
Grid<Ndim> & Grid<Ndim>::operator=(const array4<double> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<nnx;i++)
			for(unsigned j=0;j<nny;j++)
				for(unsigned k=0;k<nnz;k++)
					(*this)[l][i][j][k]=a[l][i][j][k];
				}
	return *this;
}

template <unsigned int Ndim>
void Grid<Ndim>::TriLinear(Real gx, Real gy, Real gz, Real q[Ndim]){
	int n;
	Real dd;
	for(n=0;n<VERTEX;n++){
		dd=TriDiff(gx,cube[n][XX])*TriDiff(gy,cube[n][YY])*TriDiff(gz,cube[n][ZZ]);
		for(unsigned int o=0;o<Ndim;o++)
			qq[n][o]=dd*q[o];
	}
}

template <unsigned int Ndim>
Grid<Ndim>  operator-(const Grid<Ndim> & y){
	Grid<Ndim> temp;
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<y.nnx;i++)
			for(unsigned j=0;j<y.nny;j++)
				for(unsigned k=0;k<y.nnz;k++)
					temp[l][i][j][k]=-y[l][i][j][k];
				}
	return temp;
}
template <class T, unsigned int Ndim>
	Grid<Ndim> operator*(const T & a,const Grid<Ndim> & y){
	Grid<Ndim> temp;
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<y.nnx;i++)
			for(unsigned j=0;j<y.nny;j++)
				for(unsigned k=0;k<y.nnz;k++)
					temp[l][i][j][k]=static_cast<double> (a)*y[l][i][j][k];
				}
	return temp;
}

template<unsigned int Ndim>
void Grid<Ndim>::Filter(){
	try{
	if(!filter) throw "Warning filter must be initialised ";
	}
	catch(const char * a){
		if(!Fcounter) std::cout << a << std::endl;
		Fcounter=1;
		return;
	}
	array3<Complex> & Mfilter=*filter;

	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	array3<Complex> ro[Ndim];

	for(unsigned no=0;no<Ndim;no++){
		ro[no].Allocate(nnx,nny,nnz);
		for(int i=0;i<nx0;i++)
			for(int j=0;j<ny0;j++)
				for(int k=0;k<nz0;k++){
					ro[no][i][j][k].real()=(*this)[no][i][j][k];
					ro[no][i][j][k].imag()=0.0;
				}
	}

		fft3d Forward3(nnx,nny,nnz,-1);
		fft3d Backward3(nnx,nny,nnz,1);


		for(unsigned no=0;no<Ndim;no++){
			Forward3.fft(ro[no]);
			for(int i=0;i<nx0;i++)
				for(int j=0;j<ny0;j++)
					for(int k=0;k<nz0;k++){
						Complex a;
						ro[no][i][j][k]=ro[no][i][j][k]*Mfilter[i][j][k];
					}
		}

		for(unsigned no=0;no<Ndim;no++){
			Backward3.fftNormalized(ro[no]);
			for(int i=0;i<nx0;i++)
				for(int j=0;j<ny0;j++)
					for(int k=0;k<nz0;k++)
						(*this)[no][i][j][k]=ro[no][i][j][k].real();
		}
}


template<unsigned int Ndim>
void Grid<Ndim>::MakeFilter(int mx,int my, int mz){
	if(!filter) {
		filter=new array3<Complex>;
		filter->Allocate(nnx,nny,nnz);
	}
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;
	array3<Complex> & Mfilter=*filter;

	fft3d Forward3(nnx,nny,nnz,-1);
	for(unsigned int i=0;i<nnx;i++){
		for(unsigned int j=0;j<nny;j++){
			for(unsigned int k=0;k<nnz;k++){
				Mfilter(i,j,k).real()=0.0;
				Mfilter(i,j,k).imag()=0.0;
			}
		}
	}

	double fact=1.0/(double) ((2*mx+1)*(2*my+1)*(2*mz+1));
	for(int i=-mx;i<=mx;i++){
		int ia=(i<0)?i+nnx:i;
		for(int j=-my;j<=my;j++){
			int ja=(j<0)?j+nny:j;
			for(int k=-mz;k<=mz;k++){
				int ka=(k<0)?k+nnz:k;
				Mfilter[ia][ja][ka].real()=fact;
				Mfilter(ia,ja,ka).imag()=0;
			}
		}
	}

	Forward3.fft(Mfilter);
}
template<unsigned int Ndim>
void Grid<Ndim>::MakeFilter(const double sigma){
	unsigned int nzp=nnz/2+1;
	size_t align=sizeof(Complex);
	if(!filter && !filterrc) {
		filter=new array3<Complex>;
		filter->Allocate(nnx,nny,nnz);
		filterrc=new array3<Complex>;
		filterrc->Allocate(nnx,nny,nzp,align);
	}
	int nfx,nfy,nfz;
	nfx=(nnx % 2 == 0)? nnx/2: nnx/2+1;
	nfy=(nny % 2 == 0)? nny/2: nny/2+1;
	nfz=(nnz % 2 == 0)? nnz/2: nnz/2+1;
	array3<Complex> & Mfilter=*filter;
	array3<Complex> & Mfilterrc=*filterrc;

	fft3d Forward3(nnx,nny,nnz,-1);
	fft3d Backward3(nnx,nny,nnz,1);
	Mfilter=Complex(0.0,0.0);

	double dx=co[0][0]/static_cast<double> (nnx);
	double dy=co[1][1]/static_cast<double> (nny);
	double dz=co[2][2]/static_cast<double> (nnz);
	int mx=static_cast<int> (3.0*sigma/dx);
	int my=static_cast<int> (3.0*sigma/dy);
	int mz=static_cast<int> (3.0*sigma/dz);

	double xa,ya,za;
	double sum=0.0;
	for(int i=-mx;i<=mx;i++){
		int ia=(i<0)?i+nnx:i;
		double xd=static_cast<double> (i)/static_cast<double> (nnx);
		for(int j=-my;j<=my;j++){
			int ja=(j<0)?j+nny:j;
			double yd=static_cast<double> (j)/static_cast<double> (nny);
			for(int k=-mz;k<=mz;k++){
				int ka=(k<0)?k+nnz:k;
				double zd=static_cast<double> (k)/static_cast<double> (nnz);
				xa=co[XX][XX]*xd+co[XX][YY]*yd+co[XX][ZZ]*zd;
				ya=co[YY][XX]*xd+co[YY][YY]*yd+co[YY][ZZ]*zd;
				za=co[ZZ][XX]*xd+co[ZZ][YY]*yd+co[ZZ][ZZ]*zd;
				double rs=xa*xa+ya*ya+za*za;
				double fact1=exp(-rs/(2.0*sigma*sigma));
				Mfilter[ia][ja][ka].real()=fact1;
				Mfilter(ia,ja,ka).imag()=0;

				sum+=fact1;

			}
		}
	}

	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++){
				Mfilter(i,j,k).real()=Mfilter(i,j,k).real()/sum;
			}


	array3<double> ro(nnx,nny,nnz,align);
	array3<Complex> rok(nnx,nny,nzp,align);

	rcfft3d Forward(nnz,ro,Mfilterrc);
	crfft3d Backward(nnz,Mfilterrc,ro);

	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);

	for(int i=0;i<nx0;i++)
		for(int j=0;j<ny0;j++)
			for(int k=0;k<nz0;k++){
				ro[i][j][k]=Mfilter[i][j][k].real();
			}
	Forward.fft(ro,Mfilterrc);
	Forward3.fft(Mfilter);



}

//	real sum=0;
//			sum+=qq[n];



template <>
inline Grid<1> & Grid<1>::operator=(const array3<double> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned int i=0;i<nnx;i++)
		(*this)[0][i]=a[i];
	return *this;
}
template <>
inline Grid<1>::Grid(array3<Complex> & x){
	Allocate();
	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++)
				(*this)[0][i][j][k]=x[i][j][k].real();
}

template <>
inline Grid<1>::Grid(array3<double> & x){
	Allocate();
	(*this)[0]=x;
}


template <>
inline Grid<1> & Grid<1>::operator=(const array3<Complex> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}

	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++){
				(*this)[0][i][j][k]=a[i][j][k].real();
			}
	return *this;
}

#endif
