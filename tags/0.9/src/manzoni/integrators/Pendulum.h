#ifndef PENDULUM_H
#define PENDULUM_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <string>
#include <vector>
                               
#include <definitions.h>
#include <utils/XMLInputParser.h>
#include <distributions/ParticlesDistribution.h>
#include <integrators/BaseIntegrator.h>

class Pendulum : public BaseIntegrator
{
	private:
		int order;
		static double coef_a[];
		static double coef_b[];
		int variation_type;
	    
	protected:

	public:
		Pendulum(double,int);
		~Pendulum();
		virtual void integrate(double&, double&,std::vector<double>, int, int,int);
    void integrate1(double&, double&,std::vector<double>, int, int);
    void integrate2(double&, double&,std::vector<double>, int, int);
    void integrate3(double&, double&,std::vector<double>, int, int);
		void integrate4(double&, double&,std::vector<double>, int, int);
		virtual double getEnergy(double, double, std::vector<double>, int, int, int);
};

#endif
