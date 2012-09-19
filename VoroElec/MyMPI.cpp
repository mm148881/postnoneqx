/*
 * MyMPI.cpp
 *
 *  Created on: Mar 16, 2012
 *      Author: marchi
 */
#include "MyMPI.h"
namespace Parallel {
#ifdef __MPI
	MPI::Intracomm MyMPI::comm=MPI::COMM_WORLD;
#else
	void * comm=NULL;
#endif
}

