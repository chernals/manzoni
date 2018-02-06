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

#ifndef HENON_FLOW_H
#define HENON_FLOW_H

#include <vector>
#include <string>
                    
#include <definitions.h>
#include <flows/Flow.h>

class HenonFlow : public Flow
{
	private:

	public:
		HenonFlow();
		~HenonFlow();
		void run();
};

#endif
