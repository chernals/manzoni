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

#ifndef PROC_MMA_H
#define PROC_MMA_H

#include <string>
#include <vector>
#include <iostream>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>
#include <utils/MathLink.h>

class ProcMMA
{
	private:
    std::string file;
    std::string xpath;
    ParticlesDistribution* set;
    // Methods
    void readParameters();
    std::string constructData();

	public:
		ProcMMA(ParticlesDistribution*, std::string);
		~ProcMMA();

};

#endif

