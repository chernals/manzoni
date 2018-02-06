#ifndef PARAMETERS_EVOLUTION_H
#define PARAMETERS_EVOLUTION_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <fstream>

#include <ParticlesSet.h>
#include <types/CurveParameters.h>
#include <utils/Logger.h>

class ParametersEvolution
{
	private:      
	    //
	    // Data
	    //
	    int turns;
	    
	protected:

	public:
        ParametersEvolution();
        ~ParametersEvolution();
        char readEvolvingVariable(int, std::string, CurveParameters[]);
        char readEvolvingVariable(std::string p, CurveParameters&);
        char setToZero(int, CurveParameters[]);
	    void computeDerivedParametersPolynomialCurve(CurveParameters&);
	    void computeDerivedParametersChirpCurve(CurveParameters&);
        void computeDerivedParametersFileCurve(CurveParameters&);
        
};

#endif
