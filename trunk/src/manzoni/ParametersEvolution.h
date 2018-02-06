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

#ifndef PARAMETERS_EVOLUTION_H
#define PARAMETERS_EVOLUTION_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <fstream>
#include <string>
#include <iomanip>

#include <distributions/ParticlesSet.h>
#include <types/CurveParameters.h>
#include <utils/Logger.h>

class ParametersEvolution
{
	private:      
	    int turns;
			std::string variable_path;
	    
	protected:

	public:
        ParametersEvolution(std::string, std::string);
        ~ParametersEvolution();
        char readEvolvingVariable(int, std::string, CurveParameters[]);
        char readEvolvingVariable(std::string p, CurveParameters&);
        char setToZero(int, CurveParameters[]);
	    	void computeDerivedParametersPolynomialCurve(CurveParameters&);
	    	void computeDerivedParametersChirpCurve(CurveParameters&);
        void computeDerivedParametersFileCurve(CurveParameters&);
        
};

#endif
