/*
 * EpsilonVor.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: marchi
 */

#include "EpsilonVor.h"
#include <boost/algorithm/string.hpp>
#include "Parameters.h"

namespace EpsilonNS {
	static vector<Dvect> E1;
	static vector<Dvect> M1;
	static vector<Dvect> Vzero;
	static vector<int> Mycount;
	static vector<int> izero;
	int EpsilonVor::nresid=0;

	EpsilonVor & EpsilonVor::operator-=(const Epsilon & z) {
		const EpsilonVor & z1=dynamic_cast<const EpsilonVor &> (z);
		for(unsigned int q=0;q<nresid;q++){
			E[q]=E[q]-z1.E[q];
			M[q]=M[q]-z1.M[q];
		}
		return *this;
	}
	EpsilonVor & EpsilonVor::operator=(const Epsilon & z){
		const EpsilonVor & z1=dynamic_cast<const EpsilonVor &> (z);
		Matrix co0=z1.co;
		nresid=z1.nresid;
		Deallocate();
		is_xRefset=true;
		is_COset=true;
		setMetric(co0);
		Allocate();
		ResLabel=z1.ResLabel;
		TotalCount=z1.TotalCount;
		M_avg=z1.M_avg;
		VoroVol=z1.VoroVol;
		D=z1.D;
		G=z1.G;
		M=z1.M;
		E=z1.E;
		return *this;
	}

	bool EpsilonVor::operator()(Polarization & P, Field & Ef, Matrix & Rot, Voronoi * v){
		try{
			if(!v) throw " EpsilonVor: Voronoi class not initialized ";
		}
		catch(const char * s) {
			std::cout << s << std::endl;
			exit(-1);
		}
		if(!TotalCount) {
			M1=vector<Dvect>(nresid);
			E1=vector<Dvect>(nresid);
			Mycount=vector<int>(nresid);
			Vzero=vector<Dvect>(nresid,0.0);
			izero=vector<int>(nresid,0.0);
		}

		nnx=P.getnnx();nny=P.getnny();nnz=P.getnnz();
		Dvect P_t;
		P_t=Rot*P.GetTotDip();
		Dvect pp,ee,ppo,eeo;

		M1=Vzero;E1=Vzero;Mycount=izero;
		double dx=1.0/static_cast<double> (nnx);
		double dy=1.0/static_cast<double> (nny);
		double dz=1.0/static_cast<double> (nnz);
		double GridVol=getDV();
		double xd,yd,zd,xc,yc,zc;
		for(int i=0;i<nnx;i++){
			xd=static_cast<double> (i)*dx;
			for(int j=0;j<nny;j++){
				yd=static_cast<double> (j)*dy;
				for(int k=0;k<nnz;k++){
					zd=static_cast<double> (k)*dz;
					xc=co[XX][XX]*xd+co[XX][YY]*yd+co[XX][ZZ]*zd;
					yc=co[YY][XX]*xd+co[YY][YY]*yd+co[YY][ZZ]*zd;
					zc=co[ZZ][XX]*xd+co[ZZ][YY]*yd+co[ZZ][ZZ]*zd;
					double xu,yu,zu;
					int iD=v->VoronoiCellID(xc,yc,zc);
					int myres=Residue::ResiduePid(iD);
					if(myres>-1) {
						pp=P[i][j][k];
						pp*=GridVol;
						ee=Ef[i][j][k];
						Mycount[myres]++;
						for(int o=0;o<DIM;o++) M1[myres][o]+=pp[o];
						for(int o=0;o<DIM;o++) E1[myres][o]+=ee[o];
					}else{exit(1);}
				}
			}
		}
		for(unsigned int q=0;q<nresid;q++){
			double factN=1.0/static_cast<double> (Mycount[q]);
			VoroVol[q]+=v->getVolR(q);
			E1[q]*=factN;
			pp=Rot*M1[q];

			ee=Rot*E1[q];
			E0[q]=ee;
			M0[q]=pp;
			D0[q]=pp % P_t;
			G0[q]=ee % P_t;
		}
		for(unsigned int q=0;q<nresid;q++){
			E[q]+=E0[q];
			M[q]+=M0[q];
			D[q]+=D0[q];
			G[q]+=G0[q];
		}
		M_avg+=P_t;
		TotalCount++;
		return true;
	}
	void EpsilonVor::testDipole(int nres,Dvect x){
		cout << nres << " " << x << "   " << x.Norm() << "\n";
	}
	Matrix EpsilonVor::TestD(int n){
		Matrix Dd;
		Dd=D0[n];
		return Dd;
	}
	Matrix EpsilonVor::TestG(int n){
		Matrix Gg;
		Gg=G0[n];
		return Gg;
	}
	Dvect EpsilonVor::TestM(int n){
		Dvect Mm;
		Mm=M0[n];
		return Mm;
	}
	Dvect EpsilonVor::TestE(int n){
		Dvect Ee;
		Ee=E0[n];
		return Ee;
	}
	void EpsilonVor::TestConvergence(int nu){
		double fact1=(Parameters::efact/1000.0/Parameters::kT300)/Parameters::unitc/10.0;
		Matrix Dd,Gg;
		Dvect Mm,Ee,mtot;
		double Vols;
		Dd=D[nu];
		Gg=G[nu];
		Mm=M[nu];
		Ee=E[nu];
		mtot=M_avg;
		Vols=VoroVol[nu]/static_cast<double> (TotalCount);
		Gg/=static_cast<double> (TotalCount);
		Ee/=static_cast<double> (TotalCount);
		Dd/=static_cast<double> (TotalCount);
		Mm/=static_cast<double> (TotalCount);
		mtot/=static_cast<double> (TotalCount);
		double D_fl=(Dd[XX][XX]+Dd[YY][YY]+Dd[ZZ][ZZ]-Mm*mtot)/3.0/Vols;
		double G_fl=(Gg[XX][XX]+Gg[YY][YY]+Gg[ZZ][ZZ]-Ee*mtot)/3.0;
		cout <<"\n"<< D_fl*fact1 << " " << 1+G_fl*fact1  << endl;

	}
	void EpsilonVor::WriteIt(std::ofstream & fout) const {
		if(!Parallel::comm->Get_Rank()){
			fout.write((char *) &nnx, sizeof nnx);
			fout.write((char *) &nny, sizeof nny);
			fout.write((char *) &nnz, sizeof nnz);
			fout.write((char *) &co, sizeof co);
			fout.write((char *) &co, sizeof oc);
			fout.write((char *) &nresid, sizeof nresid);
			stringstream iss;
			for(int i=0;i<ResLabel.size();i++)
				i==ResLabel.size()?iss << ResLabel[i]:iss << ResLabel[i] << " ";

			string token=iss.str();
			int ntoken=token.size();
			const char * css=new char [ntoken];
			css=token.c_str();
			fout.write((char *) &ntoken, sizeof ntoken);
			fout.write(css, (sizeof css[0])*ntoken);
		}
		Parallel::comm->Barrier();
		vector<Matrix> D1(D);
		vector<Matrix> G1(G);
		vector<Dvect>  M1(M);
		vector<Dvect>  E1(E);
		vector<double> Vols(VoroVol);
		vector<string> lab(ResLabel);
		Dvect M1_avg=M_avg;

		int Tcount=TotalCount;
		Parallel::comm->ReduceSum(&Tcount,1);
		if(!Parallel::comm->Get_Rank())	fout.write((char *) &Tcount, sizeof Tcount);
		Parallel::comm->ReduceSum(&M1_avg[0],DIM);
		Parallel::comm->ReduceSum(&Vols[0],nresid);
		Parallel::comm->ReduceSum(&D1[0][0][0],DIM*DIM*nresid);
		Parallel::comm->ReduceSum(&G1[0][0][0],DIM*DIM*nresid);
		Parallel::comm->ReduceSum(&M1[0][0],DIM*nresid);
		Parallel::comm->ReduceSum(&E1[0][0],DIM*nresid);
		if(!Parallel::comm->Get_Rank())	{
			double fact0=1.0/static_cast<double> (Tcount);
			M1_avg*=fact0;
			for(int i=0;i<nresid;i++){
				double fact1=1.0/static_cast<double> (Tcount);
				Vols[i]*=fact1;
				D1[i]*=fact1;
				G1[i]*=fact1;
				M1[i]*=fact1;
				E1[i]*=fact1;
			}

			fout.write((char *) &M1_avg[0], (sizeof M1_avg[0])*DIM);
			fout.write((char *) &Vols[0], (sizeof Vols[0])*nresid);
			fout.write((char *) &D1[0][0][0], (sizeof D1[0][0][0])*DIM*DIM*nresid);
			fout.write((char *) &G1[0][0][0], (sizeof G1[0][0][0])*DIM*DIM*nresid);
			fout.write((char *) &M1[0][0], (sizeof M1[0][0])*DIM*nresid);
			fout.write((char *) &E1[0][0], (sizeof E1[0][0])*DIM*nresid);
		}
		Parallel::comm->Barrier();
	}
	void EpsilonVor::ReadIt(std::ifstream & fin){
		Matrix co0,oc0;
		fin.read((char *) &nnx, sizeof nnx);
		fin.read((char *) &nny, sizeof nny);
		fin.read((char *) &nnz, sizeof nnz);
		fin.read((char *) &co0, sizeof co0);
		fin.read((char *) &oc0, sizeof oc0);
		fin.read((char *) &nresid, sizeof nresid);
		Deallocate();
		is_xRefset=true;
		is_COset=true;
		setMetric(co0);
		Allocate();
		int ntoken;
		fin.read((char *) &ntoken, sizeof ntoken);
		char * css=new char [ntoken];
		fin.read(css, (sizeof css[0])*ntoken);
		string token(css);
		ResLabel.clear();

		boost::split(ResLabel,token, boost::is_any_of("\t "));

		fin.read((char *) &TotalCount, sizeof TotalCount);
		fin.read((char *) &M_avg[0], (sizeof M_avg[0])*DIM);
		fin.read((char *) &VoroVol[0], (sizeof VoroVol[0])*nresid);
		fin.read((char *) &D[0][0][0], (sizeof D[0][0][0])*DIM*DIM*nresid);
		fin.read((char *) &G[0][0][0], (sizeof G[0][0][0])*DIM*DIM*nresid);
		fin.read((char *) &M[0][0], (sizeof M[0][0])*DIM*nresid);
		fin.read((char *) &E[0][0], (sizeof E[0][0])*DIM*nresid);
	}



} /* namespace EpsilonNS */
