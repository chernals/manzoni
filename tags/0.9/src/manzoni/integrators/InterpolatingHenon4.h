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

#ifndef INTERPOLATING_HENON_4_H
#define INTERPOLATING_HENON_4_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <string>
#include <vector>
                               
#include <definitions.h>
#include <distributions/ParticlesDistribution.h>
#include <integrators/BaseIntegrator.h>

class InterpolatingHenon4 : public BaseIntegrator
{
	private:
		int order;
		static double coefa;
		static double coefb;
		
		double computeA(std::vector<double>, int, int, int);
		double computeOmega2(std::vector<double>, int, int, int);
		double computeTune(std::vector<double>, int, int, int);
	    
	protected:

	public:
		InterpolatingHenon4(double,int);
		~InterpolatingHenon4();
		virtual void integrate(double&, double&,std::vector<double>, int, int, int);
    void integrate1(double&, double&,std::vector<double>, int, int, int);
    void integrate2(double&, double&,std::vector<double>, int, int, int);
    void integrate4(double&, double&,std::vector<double>, int, int, int);
		virtual double getEnergy(double, double, std::vector<double>, int, int, int);
};

#endif
