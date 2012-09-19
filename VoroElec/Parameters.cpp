/*
 * Parameters.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: marchi
 */


#include "Parameters.h"
namespace Parameters{
	float Parameters::Input::dx=0.01;
	float Parameters::Input::cut=1.1;
	int Parameters::Input::Dir=2;
	std::string Parameters::Input::filename=" ";
	bool    Parameters::Input::verbose=false;
	bool    Parameters::Input::PrintRoundIon=true;
	float    Parameters::Input::freq=-1.0;
	bool 	Parameters::Input::bRot=true;
	bool 	Parameters::Input::bTiming=false;
	bool 	Parameters::Input::bVoronoi=false;
	bool 	Parameters::Input::bTest=false;
	bool 	Parameters::Input::bTestConv=false;
	int 	Parameters::Input::nAverage=200;
}

