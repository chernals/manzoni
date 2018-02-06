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

#ifndef SYMPLECTIC_FLOW_H
#define SYMPLECTIC_FLOW_H

#include <vector>
#include <string>
                    
#include <definitions.h>
#include <Integrator.h>
#include <flows/Flow.h>

class SymplecticFlow : public Flow
{
	private:
		Integrator* integrator;

	public:
		SymplecticFlow();
		~SymplecticFlow();
		void run();
};

#endif
