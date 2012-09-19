/*
 * MyMPI.h
 *
 *  Created on: Mar 16, 2012
 *      Author: marchi
 */

#ifndef MYMPI_H_
#define MYMPI_H_
#ifdef __MPI
#include <mpi.h>
#else
#include "sys/time.h"
#endif
#include <iostream>
#include <vector>
using std::vector;

namespace Parallel {

class MyMPI {
	bool is_parallel;
#ifdef __MPI
	static MPI::Intracomm comm;
#else
	static void * comm;
#endif
public:
#ifdef __MPI
	MyMPI(int & argc,char ** & argv):is_parallel(true){
		MPI::Init(argc,argv);
		comm=MPI::COMM_WORLD;
		if(comm.Get_rank()){
			std::cout.setstate(std::ios_base::badbit);
			fclose(stdout);
			fclose(stderr);
		}
		std::cout << "\n" << std::endl;
		for(int i=0;i<comm.Get_size();i++) std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
	};
	MyMPI(): is_parallel(true) {
		MPI::Init();
		comm=MPI::COMM_WORLD;
		if(comm.Get_rank()){
			std::cout.setstate(std::ios_base::badbit);
			fclose(stdout);
			fclose(stderr);
		}
		std::cout << "\n" << std::endl;
		for(int i=0;i<comm.Get_size();i++) std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
	};
#else
	MyMPI(int & argc, char ** & argv):is_parallel(false){};
	MyMPI():is_parallel(false){};
#endif

	int Get_Rank(){
#ifdef __MPI
		return comm.Get_rank();
#else
		return 0;
#endif
	}
	int Get_Size(){
#ifdef __MPI
		return comm.Get_size();
#else
		return 1;
#endif
	}
	bool AmI_Parallel(){return is_parallel;};

	template <class T>
	void Broadcast(T *, const int);

	template <class T>
	void ReduceSum(T *, const int);

	template <class T>
	void Gather(vector<T> & , vector<T> &);

#ifdef __MPI
	void Barrier(){comm.Barrier();}
#else
	void Barrier(){};
#endif

#ifdef __MPI
	double Time(){return MPI::Wtime();};
#else
	double Time(){
        timeval tim;
        double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
		return temp;
	}
#endif
#ifdef __MPI
	virtual ~MyMPI(){MPI::Finalize();};
#else
	virtual ~MyMPI(){};
#endif
};


} /* namespace Parallel */

#endif /* MYMPI_H_ */
