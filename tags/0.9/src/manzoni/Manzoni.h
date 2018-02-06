/**
* Manzoni
*
* Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
* 
* European Organization for Nuclear Research
*
* Copyright (c) 2010+ CERN. All rights reserved.
*
**/

#ifndef MANZONI_H
#define MANZONI_H

#include <flows/Flow.h>
#include <flows/HenonFlow.h>
#include <flows/SymplecticFlow.h>
#include <utils/Logger.h>

class Manzoni
{
	private:
		std::string sim_name;
		unsigned short int flow_id;
		Flow* flow;
		
		void createFlow();

	public:
		Manzoni();
		~Manzoni();
};

#endif
