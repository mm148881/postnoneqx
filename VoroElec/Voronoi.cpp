/*
 * Voronoi.cpp
 *
 *  Created on: Jun 20, 2011
 *      Author: marchi
 */

#include "Voronoi.h"
int Voronoi::nresid=0;
int Voronoi::nr=0;
int Voronoi::nc=0;
float Voronoi::time=0.0;
string * Voronoi::label=NULL;
int * Voronoi::types=NULL;


Voronoi::Voronoi() {
	// TODO Auto-generated constructor stub

}
Voronoi::Voronoi(int n, char *** s, bool bH): Mycon(NULL),porder(NULL),	TypesName(NULL),
		Vol(NULL),Neighs(NULL),Surface(NULL){
	nr=n;
	Vol=new double [nr];
	Neighs=new vector<int> [nr];
	Surface=new vector<double> [nr];
	label=new string [nr];
	types=new int [nr];

	for(int i=0;i<nr;i++){
		label[i]=*s[i];
		if(bH || label[i].compare(0,1,h,0,1)) cindex.push_back(i);
	}
}
inline bool compare(Atoms::plane x,Atoms::plane y){
	return x.xc[3] < y.xc[3];
}
void Voronoi::Start(float frame, Atoms & atm){
	time=frame;
	const int nx=NNN,ny=NNN,nz=NNN;
	const matrix & co=atm.getMt().getCO();
	double bx=co[0][0], bxy=co[0][1], by=co[1][1],bxz=co[0][2],byz=co[1][2],bz=co[2][2];
	if(Mycon) delete Mycon;
	if(porder) delete porder;
	Mycon=new container_periodic(bx,bxy,by,bxz,byz,bz,nx,ny,nz,8);
	porder=new particle_order;
	for(unsigned int o=0;o<cindex.size();o++){
		double x=atm[cindex[o]][XX];
		double y=atm[cindex[o]][YY];
		double z=atm[cindex[o]][ZZ];
		Mycon->put(*porder,o,x,y,z);
	}
	c_loop_order_periodic vl(*Mycon,*porder);
	voronoicell_neighbor c;

	double vol0;
	int ia=0;
	if(vl.start())
		do {
			if(Mycon->compute_cell(c,vl)) {
				Vol[cindex[ia]]=c.volume();
				ia++;
			}
		} while(vl.inc());

	for(int o=0;o<nresid;o++) {
		vector<int> & cindex0=ResidueCM::getind(o);
		Vols[o]=0.0;
		for(unsigned int ia=0;ia<cindex0.size();ia++){
			int i=cindex0[ia];
			Vols[o]+=Vol[i];
		}
	}
}
void Voronoi::gather(vector<int> & it){
	vector<int> tmp=it;
	for(unsigned int o=0;o<it.size();o++) {
		tmp[o]=cindex[it[o]];
	}
	it=tmp;
	tmp.clear();
}
void Voronoi::setTypes(int nc0, const AtomIndex cidx[], const string * tname){
	nc=nc0;
	nresid=ResidueCM::Size();
	area.Allocate(nresid,nc);
	Vols.Allocate(nresid);
	TypesName=new string [nc];
	for(int o=0;o<nc;o++){
		TypesName[o]=tname[o];
		for(int ia=0;ia<cidx[o].getN();ia++)
			types[cidx[o][ia]]=o;
	}
	for(int o=0;o<ResidueCM::getiDSize();o++)
		RealResidue.push_back(ResidueCM::getiD(o));
}
void Voronoi::getData(){
	c_loop_order_periodic vl(*Mycon,*porder);
	voronoicell_neighbor c;

	double vol0;
	vector<int> nei;
	vector<double> area0;
	int ia=0;
	if(vl.start())
		do {
			if(Mycon->compute_cell(c,vl)) {
				vol0=c.volume();
				c.neighbors(nei);
				c.face_areas(area0);
				gather(nei);
				Vol[cindex[ia]]=vol0;
				Neighs[cindex[ia]]=nei;
				Surface[cindex[ia]]=area0;
				ia++;
			}
		} while(vl.inc());
	nei.clear();
	area0.clear();

	area=0.0;
	for(int o=0;o<nresid;o++) {
		vector<int> & cindex=ResidueCM::getind(o);
		double sum_v=0.0;
		for(unsigned int ia=0;ia<cindex.size();ia++){
			int i=cindex[ia];
			sum_v+=Vol[i];
			for(unsigned int p=0;p<Neighs[i].size();p++)
				area[o][getTypes(Neighs[i][p])]+=Surface[i][p];
		}
		Vols[o]=sum_v;
	}
	time++;
}
std::ofstream & operator<<(std::ofstream & fout, Voronoi & y){
	VoronoiPrint g;

	float time=y.time;
	fout << "######>> At step No. " << setw(10) << setprecision(2) << fixed<< time << endl;
	for(int o=0;o<y.nresid;o++) {
		if(!g.bPrintVols)continue;
		if(y.getTypesRes(o) != g.pGroup && g.pGroup != -1) continue;
		double a=y.Vols[o]*1000.0;
		string l=ResidueCM::getlabels(o);
		int rres=y.RealResidue[o]+1;
		(o%5)?fout << setw(10) << setprecision(4) << fixed<<  a << ' ' << setw(4) << l << ' ' << setw(5) << rres << ' ':
		  fout << endl << "%$VolRes " << setw(10) << setprecision(4) << fixed << a << ' ' << setw(4) << l  <<' ' << setw(5) << rres << ' ';
	}
	fout << endl;

	array2<double> interface;
	interface.Allocate(y.nc,y.nc);
	interface=0.0;
	array1<double> VolSel;
	VolSel.Allocate(y.nc);
	VolSel=0.0;
	for(int o=0;o<y.nresid ;o++){
		int o_type=y.getTypesRes(o);
		VolSel[o_type]+=y.Vols[o]*1000.0;
		for(int p=0;p<y.nc;p++) {
			double a=y.area[o][p]*100.0;
			interface[o_type][p]+=a;
		}
	}

	if(g.bPrintAreas ){
		fout << "# Area format:        ";
		for(int o=0;o<y.nc;o++) fout << setw(1) << fixed << " areatype(" << o << ") ";
		fout << endl;
		fout  << "%$AreaRes " ;
		for(int o=0;o<y.nresid ;o++){
			string l=ResidueCM::getlabels(o);
			int o_type=y.getTypesRes(o);
			if(g.pGroup != -1 && o_type != g.pGroup) continue;
			fout  << right << setw(5)<< l << " " << setw(5) << fixed << y.RealResidue[o]+1 << ' ' << setw(3) << fixed << o_type << ' ' ;
			for(int p=0;p<y.nc;p++) {
				double a=y.area[o][p]*100.0;
				fout << setw(10) << setprecision(5) << fixed << a << ' ';
			}
			if(o+1==y.nresid) fout << endl;
			else if(!((o+1)%2)) fout << endl << "%$AreaRes " ;
		}
	}

	fout << endl;
	fout << "# Volume of selection " << endl <<"#      ";
	for(int o=0;o<y.nc;o++) fout << setw(1) << fixed << " Volumetype(" << o << ")  ";
	fout << endl;
	fout << "%$TotVol ";
	for(int o=0;o<y.nc;o++){
			fout << setw(15) << setprecision(4) << fixed << VolSel[o];
	}
	fout << endl;

	fout << "# Interface area of selection " << endl << "#             ";
	for(int o=0;o<y.nc;o++) fout << setw(1) << fixed << "  areatype(" << o << ")";
	fout << endl;
	for(int o=0;o<y.nc;o++){
		fout << "%$AreaTot " << setw(3) << o ;
		for(int p=0;p<o;p++)
			fout <<setw(13) << setprecision(3) << fixed << ' ';

		for(int p=o;p<y.nc;p++)
			fout << setw(13) << interface[o][p] ;
		fout << endl;
	}

	return fout;
}
