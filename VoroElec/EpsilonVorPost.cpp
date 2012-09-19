/*
 * EpsilonVorPost.cpp
 *
 *  Created on: Mar 30, 2012
 *      Author: marchi
 */

#include "EpsilonVorPost.h"

namespace EpsilonNS {
EpsilonVorPost::EpsilonVorPost() {
	// TODO Auto-generated constructor stub

}
void EpsilonVorPost::Rdf(double dummy, double dummy2){
	cout << ResLabel[0] << endl;
	Data=vector<IndexedMatrix>(nresid);
	vector<Matrix> DM(nresid),GM(nresid);
	Matrix D_fl,G_fl,G_fli,D_fli;
	double fact1=Parameters::efact_nm/Parameters::kT300/Parameters::unitc/1000.0;
	cout << fact1<< " " << Parameters::kT300 << " "
			<< Parameters::efact_nm/1000.0 << " " << Parameters::unitc << endl;
	for(unsigned int i=0;i<nresid;i++){
		Matrix D_l=D[i];
		Matrix G_l=G[i];
		Dvect M_l=M[i];
		Dvect E_l=E[i];
		D_fl=D_l-M_l % M_avg;
		G_fl=G_l-E_l % M_avg;
		GM[i]=G_l/Parameters::unitc/10.0;
		DM[i]=D_l/Parameters::unitc/10.0/VoroVol[i];
	}

	Matrix one=Matrix(0.0).Unit();
	for(int i=0;i<nresid;i++){
		IndexedMatrix rdf;
		Matrix LocDM0=DM[i]*(Parameters::efact/1000.0/Parameters::kT300);
		Matrix LocGM0=GM[i]*(Parameters::efact/1000.0/Parameters::kT300);
		Dvect p=M[i];
		Dvect e=E[i];
		Matrix LocGM,LocDM;
		Dvect LocP,LocE;
		LocGM0+=one;
		LocGM=LocGM0;
		LocDM=LocDM0;
		LocP=p;
		LocE=e;
		Data[i].gm=LocGM;
		Data[i].dm=LocDM;
		Data[i].e=LocE;
		Data[i].p=LocP;
	}
}
struct duo{
	double m,g;
};
void EpsilonVorPost::Statistics(){
	const vector<IndexedMatrix> & Rdf=Data;
	Matrix one=Matrix(0.0).Unit();
	vector<duo> Dm;
	for(unsigned int i=0;i< nresid;i++){
		Matrix Locdm=Rdf[i].dm;
		Matrix Locgm=Rdf[i].gm;
		double dmt=(Locdm[XX][XX]+Locdm[YY][YY]+Locdm[ZZ][ZZ])/3.0;
		double dgt=(Locgm[XX][XX]+Locgm[YY][YY]+Locgm[ZZ][ZZ])/3.0;
		double Epsuni=(4.0*Parameters::pi)*dmt/dgt+1.0;
		duo pp={dmt,dgt};
		Dm.push_back(pp);
	}
	const int ncycle=200,ntot=200;
	int n=0;
	vector<vector<duo>::iterator > poolD;
	cout << "-- data start ---" << endl;
	do{
		vector<duo>::iterator ita=Dm.begin();
		for(;ita!=Dm.end();ita++) poolD.push_back(ita);
		do{
			int nend=(poolD.size() > ntot)?ntot:poolD.size();
			double ddm=0.0,ggm=0.0;
			for(int o=0;o<nend;o++){
				int p=roll_dice(0,poolD.size()-1);
				ddm+=poolD[p]->m;
				ggm+=poolD[p]->g;
				poolD.erase(poolD.begin()+p);
			}
			double fact0=1.0/static_cast<double> (nend);

			cout << ddm*fact0 <<  " " << ggm*fact0 << "\n";
		} while(poolD.size());
		n++;
	} while(n<ncycle);

}
void EpsilonVorPost::WriteIt(std::ofstream & fout) const {
	const vector<IndexedMatrix> & Rdf=Data;
	fout << "@  s1" << " legend " << "\"" << "My comm "<< "\"" << std::endl;
	Matrix one=Matrix(0.0).Unit();
	for(unsigned int i=0;i< nresid;i++){
		Matrix Locdm=Rdf[i].dm;
		Matrix Locgm=Rdf[i].gm;
		double dmt=(Locdm[XX][XX]+Locdm[YY][YY]+Locdm[ZZ][ZZ])/3.0;
		double dgt=(Locgm[XX][XX]+Locgm[YY][YY]+Locgm[ZZ][ZZ])/3.0;
		double Epsuni=(4.0*Parameters::pi)*dmt/dgt+1.0;

		double EPar=Rdf[i].e[ZZ];
		double PPar=Rdf[i].p[ZZ];
		fout << std::setprecision(6) << std::setw(13) << std::fixed << ResLabel[i];
		fout << std::setprecision(6) << std::setw(13) << std::fixed << i;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << Epsuni;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << EPar;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << PPar;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << dmt;
		fout << std::setprecision(6) << std::setw(13) << std::fixed << dgt;
		fout << std::endl;
	}
	fout << "&\n";
}

EpsilonVorPost::~EpsilonVorPost() {
	// TODO Auto-generated destructor stub
	Data.clear();
}

} /* namespace EpsilonNS */
