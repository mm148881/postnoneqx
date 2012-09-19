/*
 * PostGridEps.cpp
 *
 *  Created on: Mar 9, 2012
 *      Author: marchi
 */
#include "Gridn.h"
#include "Gridn_SpecTemplates.hpp"

#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options.hpp>
#include "EpsilonGridNE.h"
#include "EpsilonVorNE.h"
#include "MyMPI_Spec.hpp"
namespace po = boost::program_options;

using namespace std;

#include <fstream>
#include <iterator>
#include <string>
#include "Residue.h"
using std::string;
using namespace std;

#include "Parameters.h"

int main(int ac, char* av[])
{
    string fileI1,fileO,fileI2;
    bool verbose=false;
    bool Grid=false;
    bool Voronoi=false;
    bool force=false;
    vector<Residue> test;
    EpsilonNS::Epsilon * MyEps1;
    EpsilonNS::Epsilon * MyEps2;
    EpsilonNS::Epsilon * MyEps;
    try {
        po::options_description desc("Allowed options:");
        po::variables_map vm;
        desc.add_options()
            ("help,h", "produce help message")
            ("verbose,v", "verbose mode")
            ("Voronoi,V", "Voronoi mode")
            ("Grid,G", "Grid mode")
            ("input1,1", po::value<string>(&fileI1),"Input file to process")
            ("input2,2", po::value<string>(&fileI2),"Input file to process")
            ("out,o", po::value<string>(&fileO),"output file to generate")
            ("Dir,D", po::value<int>(&Parameters::Input::Dir),"Perturbation Direction ")
            ("dx,d", po::value<float>(&Parameters::Input::dx),"Histogram bin parameter ")
            ("cut,c", po::value<float>(&Parameters::Input::cut),"Historgram cutoff parameter ")
        ;

        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }
        if(!vm.count("input1") || !vm.count("input2")  ){
            cout << desc << "\n";
            return 1;
        }
        if(vm.count("verbose")) verbose=true;
        if(vm.count("Voronoi")) Voronoi=true;
        if(vm.count("Grid")) Grid=true;
        if(Voronoi && Grid) throw "You must choose between --grid or --voronoi. "
        		"Both arguments are present !!";
        if(!Voronoi && !Grid) throw "You must choose between --grid or --voronoi. "
        		"Neither of the arguments are present !!";
        if(Voronoi){
        	MyEps1=new EpsilonNS::EpsilonVorNE;
        	MyEps2=new EpsilonNS::EpsilonVorNE;
        	MyEps=new EpsilonNS::EpsilonVorNE;
        } else if(Grid){
        	MyEps1=new EpsilonNS::EpsilonGridNE;
        	MyEps2=new EpsilonNS::EpsilonGridNE;
        	MyEps=new EpsilonNS::EpsilonGridNE;
        }
        MyEps1->setdir(Parameters::Input::Dir);
        MyEps2->setdir(Parameters::Input::Dir);
        MyEps->setdir(Parameters::Input::Dir);
        if(verbose) {
        	cout << "Input1 filename was set to " << fileI1 << ".\n";
        	cout << "Input2 filename was set to " << fileI2 << ".\n";
        	cout << "Output filename was set to " <<  fileO << ".\n";
        }
    }
    catch(const char * s){
    	cerr << "error: " << s << endl;
    	return 1;
    }
    catch(std::exception & e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }

    ifstream ifs1(fileI1.c_str(),ifstream::in);
    ifstream ifs2(fileI2.c_str(),ifstream::in);
    ofstream ofs(fileO.c_str(),ofstream::out);
    do{
    	if(dynamic_cast<EpsilonNS::EpsilonGridNE *> (MyEps1))
    		cout << Parameters::Input::dx << " " << Parameters::Input::cut << endl;

    	if(!ifs1){cout << "Could not open input file: " + fileI1 + " \n";exit(1);}
    	if(!ifs2){cout << "Could not open input file: " + fileI2 + " \n";exit(1);}
    	ifs1>> *MyEps1;
    	ifs2>> *MyEps2;
    	cout << " Illo " << "\n";
    	*MyEps=*MyEps2;
    	*MyEps-=*MyEps1;
    	unsigned int nx0=MyEps1->GetNx();
    	unsigned int ny0=MyEps1->GetNy();
    	unsigned int nz0=MyEps1->GetNz();
    	Matrix CO=MyEps1->getCO();
    	cout << nx0 << " " << ny0 << " " << nz0 << endl;
    	cout << "It has been read !" << "\n";
    	MyEps->Rdf(Parameters::Input::dx,Parameters::Input::cut);
    	//MyEps->Statistics();
    	ofs << *MyEps;
    } while(ifs1.peek() != -1);

    return 0;
}
