/*
 * MyUtilClass.h
 *
 *  Created on: Jan 20, 2012
 *      Author: marchi
 */

#ifndef MYUTILCLASS_H_
#define MYUTILCLASS_H_

#include <iostream>
#include <iomanip>
#include "linalg.h"
#include <cmath>
#include "typedefs.h"
using namespace alglib;
const double eps=1.0e-4;
const bool DEBUG=false;
namespace MATRIX{
	struct Matrix;
}
namespace DVECT{
	struct Dvect{
		double x[DIM];
		Dvect(){};
		Dvect(const double x0, const double y0, const double z0){
			x[XX]=x0;
			x[YY]=y0;
			x[ZZ]=z0;
		}
		Dvect(const double y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
		}
		Dvect(const rvec y){
			for(int o=0;o<DIM;o++) x[o]=static_cast<double> (y[o]);
		}
		Dvect(const double y){
			for(int o=0;o<DIM;o++) x[o]=y;
		}
		Dvect(const int y){
			for(int o=0;o<DIM;o++) x[o]=static_cast<double> (y);
		}
		double operator*(const Dvect & y){
			double temp=0.0;
			for(int o=0;o<DIM;o++) temp+=x[o]*y.x[o];
			return temp;
		};
		Dvect  operator^(const Dvect & y){
			Dvect temp;
			temp[XX]=x[YY]*y.x[ZZ]-x[ZZ]*y.x[YY];
			temp[YY]=x[ZZ]*y.x[XX]-x[XX]*y.x[ZZ];
			temp[ZZ]=x[XX]*y.x[YY]-x[YY]*y.x[XX];
			return temp;
		};
		Dvect  operator+(const Dvect & y){
			Dvect temp;
			for(int o=0;o<DIM;o++) temp[o]=x[o]+y.x[o];
			return temp;
		};
		Dvect  operator-(const Dvect & y){
			Dvect temp;
			for(int o=0;o<DIM;o++) temp[o]=x[o]-y.x[o];
			return temp;
		};
		bool operator==(const double y) const {
			bool is_y=true;
			for(int i=0;i<DIM;i++)
				if(x[i] != y) {is_y=false;break;}
			return is_y;
		}
		bool operator!=(const double y) const {
			bool is_y=false;
			for(int i=0;i<DIM;i++)
				if(x[i] != y) {is_y=true;break;}
			return is_y;
		}
		Dvect & operator=(const Dvect & y){
			for(int o=0;o<DIM;o++) x[o]=y.x[o];
			return *this;
		};
		Dvect & operator=(const double y){
			for(int o=0;o<DIM;o++) x[o]=y;
			return *this;
		};
		Dvect & operator=(const double y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		Dvect & operator()(const double y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		Dvect & operator=(const rvec y){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		Dvect & operator()(const rvec y){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		Dvect & operator()(const double x0, const double y0, const double z0){
			x[XX]=x0;
			x[YY]=y0;
			x[ZZ]=z0;
			return *this;
		};
		Dvect & operator+=(const Dvect & y){
			for(int o=0;o<DIM;o++) x[o]+=y.x[o];
			return *this;
		};
		Dvect & operator/=(const double y){
			for(int o=0;o<DIM;o++) x[o]/=y;
			return *this;
		}
		Dvect & operator*=(const double y){
			for(int o=0;o<DIM;o++) x[o]*=y;
			return *this;
		}
		Dvect & operator()(const double y){
			for(int o=0;o<DIM;o++) x[o]=y;
			return *this;
		}
		Dvect operator/(const double y){
			Dvect temp(*this);
			for(int o=0;o<DIM;o++) temp.x[o]/=y;
			return temp;
		}
		double & operator[](const int i){
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return x[i];
		};
		const double & operator[](const int i) const {
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return x[i];
		};
		void normalize(){
			double n=1.0/sqrt(x[XX]*x[XX]+x[YY]*x[YY]+x[ZZ]*x[ZZ]);
			*this*=n;
		}
		Dvect normal(){
			Dvect tmp=*this;
			double n=1.0/sqrt(x[XX]*x[XX]+x[YY]*x[YY]+x[ZZ]*x[ZZ]);
			tmp*=n;
			return tmp;
		}
		double Norm(){return sqrt(sqr(x[XX])+sqr(x[YY])+sqr(x[ZZ]));};
		friend std::ostream & operator<<(std::ostream &, const Dvect &);
	};
	const MATRIX::Matrix operator%(const Dvect & x,const Dvect & y);
	const MATRIX::Matrix operator%(Dvect & x,Dvect & y);

	struct Dvect2rvec{
		rvec * r;
		int nr;
		Dvect2rvec(){};
		Dvect2rvec(int n):r(new rvec[n]),nr(n) {
		}

		Dvect2rvec(std::vector<Dvect> & x):r(new rvec[x.size()]),nr(x.size()){
			copyDvect(x);
		}
		void copyDvect(std::vector<Dvect> & x){
			for(size_t n=0;n<x.size();n++){
				float tmp[DIM];
				for(int o=0;o<DIM;o++) r[n][o]=static_cast<float> (x[n][o]);
			}
		}
		Dvect2rvec & operator()(std::vector<Dvect> & x){
			if(!r) r=new rvec[x.size()];
			copyDvect(x);
			return *this;
		}
		rvec & operator[](int n){return r[n];};
		rvec * pt(){return r;}
		~Dvect2rvec(){
			if(r) delete [] r;
		}
	};
}
namespace MATRIX{
	using namespace DVECT;
	struct Matrix{
		double m[DIM][DIM];
		Matrix(){};
		Matrix(const double y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=y;
		}
		Matrix(const int y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=static_cast<double> (y);
		}
		Matrix(const double y[DIM][DIM]){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=y[o][p];
		};
		Matrix(const float y[DIM][DIM]){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=y[o][p];
		};
		double Trace() const {
			double tmp=0.0;
			for(int o=0;o<DIM;o++) tmp+=(*this)[o][o];
			return tmp;
		}
		Matrix Transpose(){
			Matrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=(*this)[p][o];
			return temp;
		}
		Matrix Inversion(){
			Matrix temp;
			real_2d_array tmp;
			tmp.setlength(DIM,DIM);
			tmp.setcontent(DIM,DIM,&m[0][0]);
			double det=rmatrixdet(tmp);
			try{if(!fabs(det)) throw " While inverting G matrix, the matrix is singular. Abort!";}
			catch(const char * s){std::cout << s << " Det = " << det << std::endl;}

			if(DEBUG){
				real_1d_array w;
				real_2d_array u, vt;
				rmatrixsvd(tmp,DIM, DIM, 2, 2, 2, w, u, vt);
				for(int o=0;o<DIM;o++) std::cout << " W[" <<o<<"] " << w[o] << std::endl;
			}


			ae_int_t info;
			matinvreport rep;
			rmatrixinverse(tmp, DIM, info, rep);

			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=tmp[o][p];
			return temp;
		}
		const Matrix & Unit(){
			for(int o=0;o<DIM;o++)for(int p=0;p<DIM;p++) m[o][p]=(o==p)?1:0;return *this;
		};
		const double * operator[](const int i) const{
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return m[i];
		};
		double * operator[](const int i){
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return m[i];
		};
		Matrix operator*(const Matrix & y) const{
			Matrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++){
					temp[o][p]=0.0;
					for(int r=0;r<DIM;r++)
						temp[o][p]+=(*this)[o][r]*y[r][p];
				}
			return temp;
		};
		Dvect operator*(const Dvect & y){
			Dvect temp;
			for(int o=0;o<DIM;o++){
				temp[o]=0.0;
				for(int r=0;r<DIM;r++)
					temp[o]+=(*this)[o][r]*y[r];
			}
			return temp;
		};
		Matrix RotTensor(const Matrix & y){
			Matrix & Rot=*this;
			Matrix Rot_t=Rot.Transpose();
			Matrix temp=y*Rot_t;
			Matrix temp1=Rot*temp;
			return temp1;
		}
		const Matrix operator+(const Matrix & y){
			Matrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=m[o][p]+y.m[o][p];
			return temp;
		};
		const Matrix operator-(const Matrix & y){
			Matrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=m[o][p]-y.m[o][p];
			return temp;
		};
		bool operator==(const double y) const{
			bool is_y=true;
			for(int i=0;i<DIM;i++)
				for(int j=0;j<DIM;j++)
					if(m[i][j] != y) {is_y=false;break;}
			return is_y;
		}
		Matrix & operator=(const Matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator=(const matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator=(matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator=(const double y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y;
			return *this;
		};
		Matrix & operator=(const double  y[DIM][DIM]){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator()(const Matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator()(const matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator()(matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		Matrix & operator+=(const Matrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]+=y[o][p];
			return *this;
		};
		Matrix operator/(const double y){
			Matrix temp(*this);
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp.m[o][p]/=y;
			return temp;
		}
		Matrix & operator/=(const double y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]/=y;
			return *this;
		};
		Matrix & operator*=(const double y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]*=y;
			return *this;
		};
		Matrix operator*(const double y){
			Matrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=(*this)[o][p]*y;

			return temp;
		};
		friend Matrix operator*(const double y, const Matrix & z);
		friend std::ostream & operator<<(std::ostream &, const Matrix &);
	};
}



#endif /* MYUTILCLASS_H_ */
