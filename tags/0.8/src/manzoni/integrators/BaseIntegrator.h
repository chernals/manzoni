#ifndef BASE_INTEGRATOR_H
#define BASE_INTEGRATOR_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <string>
#include <vector>
                               
#include <definitions.h>
#include <distributions/ParticlesDistribution.h>

class BaseIntegrator
{
	private:
	    
	protected:
		double timestep;

	public:
		BaseIntegrator(double);
		~BaseIntegrator();
		virtual void integrate(double&, double&,std::vector<double>) = 0;
		virtual double getEnergy(double, double, std::vector<double>) = 0;
        
};

#endif
