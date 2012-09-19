/*
 * GridEps.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: marchi
 */

#include "GridEps.h"
#include "MyMPI_Spec.hpp"

namespace GridEpsNS {
	const bool TEST=false;
	unsigned int GridEps::nnx=0,GridEps::nny=0,GridEps::nnz=0;

	void GridEps::setMetric(Matrix & co_in){
/* Save some time by assuming lower right part is zero */
		is_COset=true;
		for(int i=0;i<DIM;i++)
			for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];

		double tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
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

	bool GridEps::operator()(Polarization & P, Field & Ef, bool test){
		try{if(!nnx || !nny || !nnz) throw "Grid not initialized for GridEps ";}
		catch(const char * s){std::cout << s << std::endl;return false;}
		Dvect P_t;
		P_t=P.Integrate();
		Dvect pp,ee;
		if(test){
			for(unsigned int i=0;i<nnx;i++){
				for(unsigned int j=0;j<nny;j++){
					for(unsigned int k=0;k<nnz;k++){
						pp=P[i][j][k];
						ee=Ef[i][j][k];
						D[i][j][k]+=pp % pp;
						G[i][j][k]+=ee % pp;
						for(int o=0;o<DIM;o++) M[i][j][k][o]+=P[i][j][k][o];
						for(int o=0;o<DIM;o++) E[i][j][k][o]+=Ef[i][j][k][o];
					}
				}
			}

		} else {
			for(unsigned int i=0;i<nnx;i++){
				for(unsigned int j=0;j<nny;j++){
					for(unsigned int k=0;k<nnz;k++){
						pp=P[i][j][k];
						ee=Ef[i][j][k];
						D[i][j][k]+=pp % P_t;
						G[i][j][k]+=ee % P_t;
						for(int o=0;o<DIM;o++) M[i][j][k][o]+=P[i][j][k][o];
						for(int o=0;o<DIM;o++) E[i][j][k][o]+=Ef[i][j][k][o];
					}
				}
			}
		}
		M_avg+=P_t;
		TotalCount++;
		return true;
	}

	void GridEps::Avg(){
		Matrix D_fl,G_fl,G_fli;
		Dvect M_avg_t;
		double fact0=1.0/static_cast<double> (TotalCount);
		M_avg_t=M_avg*fact0;
		for(unsigned int i=0;i<nnx;i++){
			for(unsigned int j=0;j<nny;j++){
				for(unsigned int k=0;k<nnz;k++){
					Matrix D_l=D[i][j][k]*fact0;
					Matrix G_l=G[i][j][k]*fact0;
					Dvect M_l=M[i][j][k]*fact0;
					Dvect E_l=E[i][j][k]*fact0;
					D_fl=D_l-M_l % M_avg_t;
					G_fl=G_l-E_l % M_avg_t;
					G_fli=G_fl.Inversion();
					eps[i][j][k]=4.0*M_PI*D_fl*G_fli+Matrix().Unit();
				}
			}
		}
	}
	void GridEps::Avg2(){
		Matrix D_fl,G_fl,G_fli,D_fli;
		Dvect M_avg_t;
		double fact0=1.0/static_cast<double> (TotalCount);
		M_avg_t=M_avg*fact0;
		for(unsigned int i=0;i<nnx;i++){
			for(unsigned int j=0;j<nny;j++){
				for(unsigned int k=0;k<nnz;k++){
					Matrix D_l=D[i][j][k]*fact0;
					Matrix G_l=G[i][j][k]*fact0;
					Dvect M_l=M[i][j][k]*fact0;
					Dvect E_l=E[i][j][k]*fact0;
					D_fl=D_l-M_l % M_avg_t;
					eps[i][j][k]=D_fl;
				}
			}
		}
	}
std::ofstream & operator<<(std::ofstream & fout, const GridEps & y){


	if(!Parallel::comm->Get_Rank()){
		fout.write((char *) &y.nnx, sizeof y.nnx);
		fout.write((char *) &y.nny, sizeof y.nny);
		fout.write((char *) &y.nnz, sizeof y.nnz);
		fout.write((char *) &y.co, sizeof y.co);
		fout.write((char *) &y.co, sizeof y.oc);
		fout.write((char *) &y.xRef, sizeof y.xRef);
	}
	Parallel::comm->Barrier();

	array3<Matrix> D;
	array3<Matrix> G;
	array3<Dvect>  M;
	array3<Dvect>  E;
	D.Allocate(y.nnx,y.nny,y.nnz);
	G.Allocate(y.nnx,y.nny,y.nnz);
	M.Allocate(y.nnx,y.nny,y.nnz);
	E.Allocate(y.nnx,y.nny,y.nnz);
	D=y.D;
	G=y.G;
	M=y.M;
	E=y.E;
	Dvect M_avg=y.M_avg;
	int ntot=y.TotalCount*Parallel::comm->Get_Size();
	double fact0=1.0/static_cast<double> (ntot);
	M_avg*=fact0;
	for(unsigned int i=0;i<y.nnx;i++){
			for(unsigned int j=0;j<y.nny;j++){
				for(unsigned int k=0;k<y.nnz;k++){
					D[i][j][k]*=fact0;
					G[i][j][k]*=fact0;
					M[i][j][k]*=fact0;
					E[i][j][k]*=fact0;
				}
			}
		}
	int TotalCount=y.TotalCount;
	Parallel::comm->ReduceSum(&TotalCount,1);
	if(!Parallel::comm->Get_Rank())	fout.write((char *) &TotalCount, sizeof TotalCount);
	Parallel::comm->ReduceSum(&D[0][0][0][0][0],DIM*DIM*y.nnx*y.nny*y.nnz);
	Parallel::comm->ReduceSum(&G[0][0][0][0][0],DIM*DIM*y.nnx*y.nny*y.nnz);
	Parallel::comm->ReduceSum(&M[0][0][0][0],DIM*y.nnx*y.nny*y.nnz);
	Parallel::comm->ReduceSum(&E[0][0][0][0],DIM*y.nnx*y.nny*y.nnz);

	if(!Parallel::comm->Get_Rank())	{
		fout.write((char *) &D[0][0][0][0][0], (sizeof D[0][0][0][0][0])*DIM*DIM*y.nnx*y.nny*y.nnz);
		fout.write((char *) &G[0][0][0][0][0], (sizeof G[0][0][0][0][0])*DIM*DIM*y.nnx*y.nny*y.nnz);
		fout.write((char *) &M[0][0][0][0], (sizeof M[0][0][0][0])*DIM*y.nnx*y.nny*y.nnz);
		fout.write((char *) &E[0][0][0][0], (sizeof E[0][0][0][0])*DIM*y.nnx*y.nny*y.nnz);
	}
	Parallel::comm->Barrier();
	return fout;
}
std::ifstream & operator>>(std::ifstream & fin, GridEps & y){
	Matrix co,oc;
	static unsigned int nnx,nny,nnz;
	int TotalCount;
	fin.read((char *) &nnx, sizeof nnx);
	fin.read((char *) &nny, sizeof nny);
	fin.read((char *) &nnz, sizeof nnz);
	fin.read((char *) &co, sizeof co);
	fin.read((char *) &oc, sizeof oc);
	fin.read((char *) &y.xRef, sizeof y.xRef);
	fin.read((char *) &y.TotalCount, sizeof y.TotalCount);
	y.Deallocate();
	y.Setup(nnx,nny,nnz);
	y.setMetric(co);
	y.is_xRefset=true;
	y.Allocate();
	fin.read((char *) &y.D[0][0][0][0][0], (sizeof y.D[0][0][0][0][0])*DIM*DIM*y.nnx*y.nny*y.nnz);
	fin.read((char *) &y.G[0][0][0][0][0], (sizeof y.G[0][0][0][0][0])*DIM*DIM*y.nnx*y.nny*y.nnz);
	fin.read((char *) &y.M[0][0][0][0], (sizeof y.M[0][0][0][0])*DIM*y.nnx*y.nny*y.nnz);
	fin.read((char *) &y.E[0][0][0][0], (sizeof y.E[0][0][0][0])*DIM*y.nnx*y.nny*y.nnz);
	return fin;
}

} /* namespace GridEpsNS */
