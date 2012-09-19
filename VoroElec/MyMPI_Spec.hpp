/*
 * MyMPI_Spec.hpp
 *
 *  Created on: Mar 16, 2012
 *      Author: marchi
 */

#ifndef MYMPI_SPEC_HPP_
#define MYMPI_SPEC_HPP_

namespace Parallel{
template<>
	void MyMPI::Gather<double>(vector<double> & bufferIn, vector<double> & bufferOut){
	int nIn=static_cast<int> (bufferIn.size());
	int nbufferIn,nbufferOut;
	nbufferIn=nIn;
#ifdef __MPI
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(&nIn,NULL,1,MPI::INT,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,&nIn,1,MPI::INT,MPI::SUM,0);
	comm.Bcast(&nIn,1,MPI::INT,0);
	nbufferOut=nIn;
	bufferOut.clear();
	bufferOut=vector<double>(nbufferOut,0.0);
	comm.Gather(&bufferIn[0],nbufferIn,MPI::DOUBLE,&bufferOut[0],nbufferIn,MPI::DOUBLE,0);
#else
	if(&bufferIn[0] == &bufferOut[0]) return;
	bufferOut.clear();
	bufferOut=bufferIn;
#endif
}

template<>
	void MyMPI::Broadcast<double>(double * buffer,const int nbuffer){
#ifdef __MPI
		comm.Bcast(buffer,nbuffer,MPI::DOUBLE,0);
#endif
	};

	template<>
	void MyMPI::Broadcast<int>(int * buffer,const int nbuffer){
#ifdef __MPI
		comm.Bcast(buffer,nbuffer,MPI::INT,0);
#endif
	};
	template<>
	void MyMPI::Broadcast<float>(float * buffer,const int nbuffer){
#ifdef __MPI
		comm.Bcast(buffer,nbuffer,MPI::FLOAT,0);
#endif
	};
	template<>
	void MyMPI::Broadcast<char>(char * buffer,const int nbuffer){
#ifdef __MPI
		comm.Bcast(buffer,nbuffer,MPI::CHAR,0);
#endif
	};

template<>
	void MyMPI::ReduceSum<float>(float * buffer,const int nbuffer){
#ifdef __MPI
		comm.Barrier();
		if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::FLOAT,MPI::SUM,0);
		else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::FLOAT,MPI::SUM,0);
#endif
	};
template<>
	void MyMPI::ReduceSum<double>(double * buffer,const int nbuffer){
#ifdef __MPI
		comm.Barrier();
		if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::DOUBLE,MPI::SUM,0);
		else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::DOUBLE,MPI::SUM,0);
#endif
	};
template<>
	void MyMPI::ReduceSum<int>(int * buffer,const int nbuffer){
#ifdef __MPI
		comm.Barrier();
		if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::INT,MPI::SUM,0);
		else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::INT,MPI::SUM,0);
#endif
	};


}



#endif /* MYMPI_SPEC_HPP_ */
