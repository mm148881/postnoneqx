/*
 * EpsilonGrid.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#include "EpsilonGrid.h"

namespace EpsilonNS {

	Matrix EpsilonGrid::zero=Matrix().Unit();
	EpsilonGrid & EpsilonGrid::operator=(const Epsilon & z){
		const EpsilonGrid & z1=dynamic_cast<const EpsilonGrid &> (z);
		Deallocate();
		nnx=z1.nnx;
		nny=z1.nny;
		nnz=z1.nnz;
		Matrix co0=z1.co;
		is_xRefset=true;
		is_COset=true;
		setMetric(co0);
		Allocate();
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++)
					for(int o=0;o<DIM;o++) {
						M[i][j][k][o]=z1.M[i][j][k][o];
						E[i][j][k][o]=z1.E[i][j][k][o];
					}
		return *this;
	}

	EpsilonGrid & EpsilonGrid::operator-=(const Epsilon & z) {
		const EpsilonGrid & z1=dynamic_cast<const EpsilonGrid &> (z);
		cout << M[0][0][0][YY] << endl;
		cout << E[0][0][0][YY] << endl;
		for(unsigned int i=0;i<nnx;i++){
			for(unsigned int j=0;j<nny;j++){
				for(unsigned int k=0;k<nnz;k++){
					for(int o=0;o<DIM;o++) {
						M[i][j][k][o]-=z1.M[i][j][k][o];
						E[i][j][k][o]-=z1.E[i][j][k][o];
					}
				}
			}
		}
		cout << M[0][0][0][YY] << endl;
		cout << E[0][0][0][YY] << endl;
		return *this;
	}

	bool EpsilonGrid::operator()(Polarization & P, Field & Ef, Matrix & Rot,  Voronoi * v){
		try{if(!nnx || !nny || !nnz) throw "Grid not initialized for EpsilonGrid ";}
		catch(const char * s){std::cout << s << std::endl;return false;}
		Dvect P_t;
		P_t=P.Integrate();
		Dvect pp,ee;
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
		M_avg+=P_t;
		TotalCount++;
		return true;
	}


	EpsilonGrid::EpsilonGrid() {

		try{
			if(nnx && nny && nnz) Allocate();
			else throw " Warning Allocation is differred !!!";
		}
		catch(const char * s) {
			std::cout << s << std::endl;
		}
	};

	EpsilonGrid::~EpsilonGrid() {Deallocate();}
	void EpsilonGrid::WriteIt(std::ofstream & fout) const {
		if(!Parallel::comm->Get_Rank()){
			fout.write((char *) &nnx, sizeof nnx);
			fout.write((char *) &nny, sizeof nny);
			fout.write((char *) &nnz, sizeof nnz);
			fout.write((char *) &co, sizeof co);
			fout.write((char *) &co, sizeof oc);
			fout.write((char *) &xRef, sizeof xRef);
		}
		Parallel::comm->Barrier();

		array3<Matrix> D1;
		array3<Matrix> G1;
		array3<Dvect>  M1;
		array3<Dvect>  E1;
		D1.Allocate(nnx,nny,nnz);
		G1.Allocate(nnx,nny,nnz);
		M1.Allocate(nnx,nny,nnz);
		E1.Allocate(nnx,nny,nnz);
		D1=D;
		G1=G;
		M1=M;
		E1=E;
		Dvect M1_avg=M_avg;
		int ntot=TotalCount*Parallel::comm->Get_Size();
		double fact0=1.0/static_cast<double> (ntot);
		M1_avg*=fact0;
		for(unsigned int i=0;i<nnx;i++){
			for(unsigned int j=0;j<nny;j++){
				for(unsigned int k=0;k<nnz;k++){
					D1[i][j][k]*=fact0;
					G1[i][j][k]*=fact0;
					M1[i][j][k]*=fact0;
					E1[i][j][k]*=fact0;
				}
			}
		}
		int Tcount=TotalCount;
		Parallel::comm->ReduceSum(&Tcount,1);
		if(!Parallel::comm->Get_Rank())	fout.write((char *) &Tcount, sizeof Tcount);
		Parallel::comm->ReduceSum(&D1[0][0][0][0][0],DIM*DIM*nnx*nny*nnz);
		Parallel::comm->ReduceSum(&G1[0][0][0][0][0],DIM*DIM*nnx*nny*nnz);
		Parallel::comm->ReduceSum(&M1[0][0][0][0],DIM*nnx*nny*nnz);
		Parallel::comm->ReduceSum(&E1[0][0][0][0],DIM*nnx*nny*nnz);

		if(!Parallel::comm->Get_Rank())	{
			fout.write((char *) &D1[0][0][0][0][0], (sizeof D1[0][0][0][0][0])*DIM*DIM*nnx*nny*nnz);
			fout.write((char *) &G1[0][0][0][0][0], (sizeof G1[0][0][0][0][0])*DIM*DIM*nnx*nny*nnz);
			fout.write((char *) &M1[0][0][0][0], (sizeof M1[0][0][0][0])*DIM*nnx*nny*nnz);
			fout.write((char *) &E1[0][0][0][0], (sizeof E1[0][0][0][0])*DIM*nnx*nny*nnz);
		}
		Parallel::comm->Barrier();
	}
	void EpsilonGrid::ReadIt(std::ifstream & fin){
		Matrix co0,oc0;
		fin.read((char *) &nnx, sizeof nnx);
		fin.read((char *) &nny, sizeof nny);
		fin.read((char *) &nnz, sizeof nnz);
		fin.read((char *) &co0, sizeof co0);
		fin.read((char *) &oc0, sizeof oc0);
		fin.read((char *) &xRef, sizeof xRef);
		fin.read((char *) &TotalCount, sizeof TotalCount);
		Deallocate();
		is_xRefset=true;
		is_COset=true;
		setMetric(co0);
		Allocate();
		cout << G.Nx() << " " << M.Nx() << " " << E.Nx()<< endl;
		fin.read((char *) &D[0][0][0][0][0], (sizeof D[0][0][0][0][0])*DIM*DIM*nnx*nny*nnz);
		fin.read((char *) &G[0][0][0][0][0], (sizeof G[0][0][0][0][0])*DIM*DIM*nnx*nny*nnz);
		fin.read((char *) &M[0][0][0][0], (sizeof M[0][0][0][0])*DIM*nnx*nny*nnz);
		fin.read((char *) &E[0][0][0][0], (sizeof E[0][0][0][0])*DIM*nnx*nny*nnz);
	}


} /* namespace EpsilonNS */
