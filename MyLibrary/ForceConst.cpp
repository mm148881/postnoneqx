/*
 * ForceConst.cpp
 *
 *  Created on: Sep 22, 2011
 *      Author: marchi
 */

#include "ForceConst.h"

Array2D<double> ForceConst::transpose(Array2D<double> y){
	int n=y.dim1();
	Array2D<double> Tr(n,n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Tr[i][j]=y[j][i];
	return Tr;
}

ForceConst::ForceConst(): Percent(1.0) {
	// TODO Auto-generated constructor stub

}
ForceConst::ForceConst(ResidueAvg & z): Percent(1.0) {
	Array4<double> y=z.getRms();
	Array2<double> x=z.getCoord();


	int n=y.Nx();
	Dist=Array2D<double>(n,n);

	Rms=Array2D<double>(n,n);
	for(int o=0;o<n;o++){
		for(int p=o;p<n;p++){
			double sum=0.0;
			for(int s=0;s<DIM;s++) sum+=y[o][p][s][s];
			Rms[o][p]=sum;
			Rms[p][o]=sum;
		}
	}
	for(int o=0;o<n;o++)
		for(int p=o+1;p<n;p++){
			Dist[o][p]=0.0;
			for(int s=0;s<DIM;s++) Dist[o][p]+=(x[o][s]-x[p][s])*(x[o][s]-x[p][s]);
			Dist[o][p]=sqrt(Dist[o][p]);
			}
}
void ForceConst::Svd(){
	int n=Rms.dim1();
	int StartLoop=static_cast<int> (Percent*n);
	Array2D<double> U, V, S;
	Array1D<double> s;
/*
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++) {
			Rms[i][j]=Rms[i][j]*100.0;
			if(Dist[i][j] > 0.6) {
				cout << i << " " << j << "  " << Dist[i][j] << endl;
				Rms[i][j]=0.0;
			}
		}
*/
	SVD<double> G(Rms);
	G.getU(U);
	G.getV(V);
	G.getS(S);
	G.getSingularValues(s);
	Array2D<double> Sm1(n,n), Ks(n,n), Rms_b(n,n), V1(n,n), U1(n,n), Id(n,n);
	for(int i=StartLoop;i<n;i++) S[i][i]=0.0;
	V1=transpose(V);
	U1=transpose(U);
	Rms_b=matmult(U,matmult(S,V1));

	Sm1=S;
	for(int i=0;i<n;i++)
		Sm1[i][i]=(Sm1[i][i] == 0.0)?0.0:1.0/Sm1[i][i];
	Ks=matmult(V,matmult(Sm1,U1));
	Id=matmult(transpose(Ks),Rms_b);

	for(int i=0;i<n;i++) printf(" %5d   %12.5e %12.5e \n",i,Rms[i][i],Rms_b[i][i]);
//	for(int i=0;i<n;i++)
//		for(int j=i;j<n;j++)
//			printf(" %5d %5d %12.4f %12.5e %12.5e \n",i,j,Dist[i][j],Ks[i][j], Rms[i][j]);
}
ForceConst::~ForceConst() {
	// TODO Auto-generated destructor stub
}
std::ofstream & operator<<(std::ofstream & fout, ForceConst & y){
	int n=y.Rms.dim1();
	int m=y.Rms.dim2();
	fout.write((char *) &n, (sizeof n));
	fout.write((char *) &m, (sizeof m));
	fout.write((char *) &y.Rms[0][0], (sizeof y.Rms[0][0])*y.Rms.dim1()*y.Rms.dim2());
	fout.write((char *) &y.Dist[0][0], (sizeof y.Dist[0][0])*y.Dist.dim1()*y.Dist.dim2());
	return fout;
}
std::ifstream & operator>>(std::ifstream & fin, ForceConst & y){
	int n,m;
	fin.read((char *) &n, (sizeof n));
	fin.read((char *) &m, (sizeof m));
	y.Rms=Array2D<double>(n,m);
	y.Dist=Array2D<double>(n,m);
	fin.read((char *) &y.Rms[0][0], (sizeof y.Rms[0][0])*y.Rms.dim1()*y.Rms.dim2());
	fin.read((char *) &y.Dist[0][0], (sizeof y.Dist[0][0])*y.Dist.dim1()*y.Dist.dim2());
	return fin;
}


