#include <integrators/Pendulum.h>

double Pendulum::coef_a[] = {0.6756035959798289, -0.17560359597982886, -0.175603596, 0.6756035959798289};
double Pendulum::coef_b[] = {0, 1.351207191959658, -1.7024143839193155, 1.351207191959658};

Pendulum::Pendulum(double ts, int o) : BaseIntegrator(ts), order(o)
{
	// Integrator coefficients 
	coef_a[0] = coef_a[0] * timestep;
	coef_a[1] = coef_a[1] * timestep;
	coef_a[2] = coef_a[2] * timestep;
	coef_a[3] = coef_a[3] * timestep;
  
	coef_b[0] = coef_b[0] * timestep;
	coef_b[1] = coef_b[1] * timestep;
	coef_b[2] = coef_b[2] * timestep;
	coef_b[3] = coef_b[3] * timestep;
}

Pendulum::~Pendulum()
{

}

void Pendulum::integrate(double& PHI, double& J,std::vector<double> params)
{
	switch(order)
	{
		case 1:
			integrate1(PHI, J, params);
		case 2:
			integrate2(PHI, J, params);
			break;
		case 3:
			integrate3(PHI, J, params);
			break;
		case 4:
			integrate4(PHI, J, params);
			break;
		default:
			throw("Invalid iterator order !!!");
	}
}

double Pendulum::getEnergy(double PHI, double J, std::vector<double> params)
{
	return 0.5 * J * J - params.at(0) * cos(PHI);
}

void Pendulum::integrate1(double& PHI, double& J,std::vector<double> params)
{
	// params.at(0) === AK
	// params.at(1) === DELTA
	// params.at(2) === BETA
	PHI = PHI + timestep * 1 *	(J-params.at(1));
	J   =	J   - timestep * params.at(0) * (1+params.at(2)) * sin(PHI);
}
	
void Pendulum::integrate2(double& PHI, double& J,std::vector<double> params)
{
	PHI = PHI + timestep * 0.5 *	(J-params.at(1));
	J   =	J   - timestep * params.at(0) * (1+params.at(2)) * sin(PHI);
	PHI = PHI + timestep * 0.5 *	(J-params.at(1));
}

void Pendulum::integrate3(double& PHI, double& J,std::vector<double> params)
{
	PHI = PHI + timestep * (7/24) *	(J-params.at(1));
	J   =	J   - timestep * (2/3) * params.at(0) * (1+params.at(2)) * sin(PHI);
	PHI = PHI + timestep * (3/4) *	(J-params.at(1));
	J   =	J   - timestep * (-2/3) * params.at(0) * (1+params.at(2)) * sin(PHI);
	PHI = PHI + timestep * (-1/24) *	(J-params.at(1));
	J   =	J   - timestep * (1) * params.at(0) * (1+params.at(2)) * sin(PHI);
}

void Pendulum::integrate4(double& PHI, double& J,std::vector<double> params)
{
	// params.at(0) === AK
	// params.at(1) === DELTA
	// params.at(2) === BETA
	for(int order=1;order <= 4; order++)
	{
		 	J   = J    - coef_b[order-1] * params.at(0) * (1+params.at(2)) * sin(PHI);
     	PHI = PHI  + coef_a[order-1] * (J-params.at(1));
	}
}
