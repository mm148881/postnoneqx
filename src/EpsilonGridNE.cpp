/*
 * EpsilonGridNE.cpp
 *
 *  Created on: Jun 28, 2012
 *      Author: marchi
 */

#include "EpsilonGridNE.h"

namespace EpsilonNS {

EpsilonGridNE::EpsilonGridNE() {
	// TODO Auto-generated constructor stub

}

void EpsilonGridNE::WriteIt(std::ofstream & fout) const {
	double dx=Parameters::Input::dx;
	double cut=Parameters::Input::cut;
	const vector<IndexedMatrix> & Rdf=RdfData;
	fout << "@  s1" << " legend " << "\"" << "My comm "<< "\"" << "  " << Parameters::unitfield << std::endl;
	for(unsigned int i=0;i< Rdf.size();i++){
		Matrix one=0.0;
		double EPar=Rdf[i].e[ZZ]*Parameters::unitfield;
		double PPar=Rdf[i].p[ZZ]*Parameters::unitfield;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << Rdf[i].i*dx*10;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << EPar;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << PPar;
		fout << std::endl;
	}
	fout << "&\n";
}
void EpsilonGridNE::Rdf(const double dx, const double cut){
		vector<IndexedMatrix> Rdf;
		RdfData.clear();
		int nrdf=static_cast<int> ((cut/dx)+1);
		vector<HistData>  rdf(nrdf);

		try{
			if(!is_COset) throw " Must provide Metric ";
			if(!is_xRefset) throw " Must provide reference point ";
		}
		catch(const char * s){
			std::cerr << s << std::endl;
			exit(1);
		}
		Dvect M_avg_t;
		Dvect P_t=0.0;
		for(int i=0;i< nnx;i++)
			for(int j=0;j<nny;j++)
				for(int k=0;k<nnz;k++)
					for(int o=0;o<DIM;o++)
						P_t[o]+=M[i][j][k][o];
		for(int o=0;o<DIM;o++) P_t[o]*=this->getDV();

		M_avg=P_t;
		double fact1=Parameters::efact_nm/Parameters::kT300/Parameters::unitc/1000.0;
		cout << fact1<< " " << Parameters::kT300 << " " << Parameters::efact_nm/1000.0 << " " << Parameters::unitc << endl;

		for(unsigned int i=0;i<nnx;i++){
			for(unsigned int j=0;j<nny;j++){
				for(unsigned int k=0;k<nnz;k++){
					Dvect M_l=M[i][j][k];
					Dvect E_l=E[i][j][k];
				}
			}
		}


		Dvect xa;
		const Dvect & x=xRef;
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
		double dist;
		Dvect xc,xd;
		for(int i=0;i<nx0;i++){
			rx=static_cast<double> ((i)*rnx);
			for(int j=0;j<ny0;j++){
				ry=static_cast<double> ((j)*rny);
				for(int k=0;k<nz0;k++){
					rz=static_cast<double> ((k)*rnz);
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


						Euler Eul(xc);
						Matrix one=Matrix(0.0).Unit();

						Dvect p=M[i][j][k];
						Dvect e=E[i][j][k];
						bool bRot=false;
						Dvect LocP,LocE;
						if(bRot){
							LocP=Eul.Rotate(p);
							LocE=Eul.Rotate(e);
						}else{
							LocP=p;
							LocE=e;
						}
						int h=static_cast<int> (dist/dx);

						rdf[h].e+=LocE;
						rdf[h].p+=LocP;
						rdf[h].idx++;
					}

				}
			}
		}
		for(unsigned int i=0;i<nrdf;++i){
			IndexedMatrix tmp;
			if(rdf[i].idx){
				tmp.p=rdf[i].p/static_cast<double> (rdf[i].idx);
				tmp.e=rdf[i].e/static_cast<double> (rdf[i].idx);
				tmp.i=i;
				Rdf.push_back(tmp);

			}
		}
		RdfData=Rdf;
	}

EpsilonGridNE::~EpsilonGridNE() {
	// TODO Auto-generated destructor stub
}

} /* namespace EpsilonNS */
