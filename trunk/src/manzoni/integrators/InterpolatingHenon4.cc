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

#include <integrators/InterpolatingHenon4.h>

// For A. Chin method
//double InterpolatingHenon4::coef_a[] = {1.81712,-1.81712};
//double InterpolatingHenon4::coef_b[] = {-0.302853,0.302853,-0.151427};

double InterpolatingHenon4::coefa = (1.0/(2.0-pow(2.0,1.0/3.0)));//1.351207191959658;
double InterpolatingHenon4::coefb = -(pow(2.0,1.0/3.0)/(2.0-pow(2.0,1.0/3.0)));//-1.7024143839193155;


InterpolatingHenon4::InterpolatingHenon4(double ts, int o) : BaseIntegrator(ts), order(o)
{
}

InterpolatingHenon4::~InterpolatingHenon4()
{

}

void InterpolatingHenon4::integrate(double& PHI, double& J, std::vector<double> params, int t, int iter, int turns)
{
	switch(order)
	{
		case 1:
			integrate1(PHI, J, params, t, iter, turns);
		case 2:		
			integrate2(PHI, J, params, t, iter, turns);
			break;
		case 4:
			integrate4(PHI, J, params, t, iter, turns);
			break;
		default:
			throw("Invalid iterator order (should be 1, 2 or 4) !!!");
	}
}

double InterpolatingHenon4::getEnergy(double PHI, double J, std::vector<double> params, int t, int iter, int turns)
{
	return (2.0 * M_PI * (params.at(2) + (params.at(0)/turns) * (t + timestep * iter ))) * J 
	       + 0.5 * computeOmega2(params,t,iter,turns) * J *J 
	       + computeA(params,t,iter,turns) * J *J * cos(4*PHI) ;
}

void InterpolatingHenon4::integrate1(double& PHI, double& J, std::vector<double> params, int t, int iter, int turns)
{	
	double A = computeA(params, t, iter, turns);
	double OMEGA2 = computeOmega2(params, t, iter, turns);
	double epst = params.at(2) + (params.at(0)/turns) * (t + timestep * iter );
	
	//
	// H = H0 + H1 = (\epsilon \rho + (\Omega_2+A) \rho^2) + (-2A) p^2 q^2
	//
	// Exp(\Delta t / 2 D_{H0})
	PHI += timestep * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
	
	// Exp(\Delta t / 2 D_{H1})
	// First change the variables, in these variables
	// the equations of motion are exactly solvable
	double q =  sqrt(2.0*J) * sin(PHI);
	double p =  sqrt(2.0*J) * cos(PHI);
	double q_tmp = q * exp(-timestep*4*A*p*q);
	p = p * exp(timestep * 4 * A * p * q);
	q = q_tmp;
	J = (p*p+q*q)/2;
	PHI = atan2(q,p);

}

void InterpolatingHenon4::integrate2(double& PHI, double& J, std::vector<double> params, int t, int iter, int turns)
{	
	double A = computeA(params, t, iter, turns);
	double OMEGA2 = computeOmega2(params, t, iter, turns);
	double epst = params.at(2) + (params.at(0)/turns) * (t + timestep * iter );
	
	//
	// H = H0 + H1 = (\epsilon \rho + (\Omega_2+A) \rho^2) + (-2A) p^2 q^2
	//
	// Exp(\Delta t / 2 D_{H0})
	PHI += 0.5 * timestep * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
	
	// Exp(\Delta t / 2 D_{H1})
	double q =  sqrt(2.0*J) * sin(PHI);
	double p =  sqrt(2.0*J) * cos(PHI);
	double q_tmp = q * exp(-timestep*4*A*p*q);
	p = p * exp(timestep * 4 * A * p * q);
	q = q_tmp;
	J = (p*p+q*q)/2;
	PHI = atan2(q,p);
	
	// Exp(\Delta t / 2 D_{H0})
	PHI += 0.5 * timestep * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
}

void InterpolatingHenon4::integrate4(double& PHI, double& J, std::vector<double> params, int t, int iter, int turns)
{/*
// For A. Chin method
	J   -= coef_b[2] * (PHI + (1.0/3.0) * PHI * PHI * PHI);
  PHI += coef_a[1] * (1.0 + 0.5 * J * J);
	J   -= coef_b[1] * (PHI + (1.0/3.0) * PHI * PHI * PHI);
	PHI += coef_a[0] * (1.0 + 0.5 * J * J);
	J   -= coef_b[0] * (PHI + (1.0/3.0) * PHI * PHI * PHI);
	PHI += coef_a[0] * (1.0 + 0.5 * J * J);
  J   -= coef_b[1] * (PHI + (1.0/3.0) * PHI * PHI * PHI);
	PHI += coef_a[1] * (1.0 + 0.5 * J * J);
	J   -= coef_b[2] * (PHI + (1.0/3.0) * PHI * PHI * PHI);
	*/

	double A = computeA(params, t, iter, turns);
	double OMEGA2 = computeOmega2(params, t, iter, turns);
	double epst = params.at(2) + (params.at(0)/turns) * (t + timestep * iter );

	//
	// H = H0 + H1 = (\epsilon \rho + (\Omega_2+A) \rho^2) + (-2A) p^2 q^2
	//
	{
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefa*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
		
		// Exp(\Delta t / 2 D_{H1})
		double q =  sqrt(2.0*J) * sin(PHI);
		double p =  sqrt(2.0*J) * cos(PHI);
		double q_tmp = q * exp(-coefa*timestep*4*A*p*q);
		p = p * exp(coefa*timestep * 4 * A * p * q);
		q = q_tmp;
		J = (p*p+q*q)/2;
		PHI = atan2(q,p);
		
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefa*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
	}
	
	{
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefb*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
		
		// Exp(\Delta t / 2 D_{H1})
		double q =  sqrt(2.0*J) * sin(PHI);
		double p =  sqrt(2.0*J) * cos(PHI);
		double q_tmp = q * exp(-coefb*timestep*4*A*p*q);
		p = p * exp(coefb*timestep * 4 * A * p * q);
		q = q_tmp;
		J = (p*p+q*q)/2;
		PHI = atan2(q,p);
		
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefb*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
	}
	
	{
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefa*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
		
		// Exp(\Delta t / 2 D_{H1})
		double q =  sqrt(2.0*J) * sin(PHI);
		double p =  sqrt(2.0*J) * cos(PHI);
		double q_tmp = q * exp(-coefa*timestep*4*A*p*q);
		p = p * exp(coefa*timestep * 4 * A * p * q);
		q = q_tmp;
		J = (p*p+q*q)/2;
		PHI = atan2(q,p);
		
		// Exp(\Delta t / 2 D_{H0})
		PHI += 0.5 * (coefa*timestep) * (2.0 * M_PI * epst + (OMEGA2 + 2 * A) * J);
	}

}

double InterpolatingHenon4::computeTune(std::vector<double> params, int t, int iter, int turns)
{
	return 2.0*M_PI*(0.25 + params.at(2) + (params.at(0)/(turns-1)) * (t + timestep * iter ));
}

double InterpolatingHenon4::computeOmega2(std::vector<double> params, int t, int iter, int turns)
{
	double tune = computeTune(params, t, iter, turns);
	return (-1.0/16.0) * (3*tan(M_PI_2 - 0.5 * tune) + tan(M_PI_2 - 1.5 * tune)) - (3.0/8.0) * params.at(1);
}

double InterpolatingHenon4::computeA(std::vector<double> params, int t, int iter, int turns)
{
	double tune = computeTune(params, t, iter, turns);
	return (tune-M_PI/2) * (1.0/16.0) * (tan(M_PI_2 - 0.5 * tune) - tan(M_PI_2 - 1.5 * tune) - 2.0 * params.at(1));
}
