/*
 * MyUtilClass.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: marchi
 */

#include "MyUtilClass.h"
namespace DVECT{
	using MATRIX::Matrix;
	std::ostream & operator<<(std::ostream & out, const Dvect & z){
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[XX] << " "
				<< z[YY] << " " << z[ZZ] << std::endl;
		return out;
	}
	const Matrix operator%(const Dvect & x,const Dvect & y){
		Matrix temp;
		for(int o=0;o<DIM;o++)
			for(int p=0;p<DIM;p++)
				temp[o][p]=x[o]*y[p];
		return temp;
	};
	const Matrix operator%(Dvect & x,Dvect & y){
		Matrix temp;
		for(int o=0;o<DIM;o++)
			for(int p=0;p<DIM;p++)
				temp[o][p]=x[o]*y[p];
		return temp;
	};

}
namespace MATRIX{

	Matrix operator*(const double y, const Matrix & z){
		Matrix temp;
		for(int o=0;o<DIM;o++)
			for(int p=0;p<DIM;p++)
				temp[o][p]=z[o][p]*y;
		return temp;
	};
	std::ostream & operator<<(std::ostream & out, const Matrix & z){
		out << "      Matrix        " << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[XX][XX] << " " << z[XX][YY] << " " << z[XX][ZZ] << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed<< z[YY][XX] << " " << z[YY][YY] << " " << z[YY][ZZ] << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed<< z[ZZ][XX] << " " << z[ZZ][YY] << " " << z[ZZ][ZZ] << std::endl;
		return out;
	}


}

