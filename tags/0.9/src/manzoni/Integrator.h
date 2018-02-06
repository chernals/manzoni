#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
                               
#include <definitions.h>
#include <distributions/ParticlesDistribution.h>
#include <distributions/TrackingDistribution.h>
#include <ParametersEvolution.h>
#include <types/types.h>
#include <utils/Logger.h>
#include <integrators/BaseIntegrator.h>
#include <integrators/Pendulum.h>
#include <integrators/InterpolatingHenon4.h>

class Flow;

class Integrator
{
	private:
	    //
	    // Methods
	    //
	    void readParameters();
	    inline double evaluateVariable(CurveParameters&);
	    void createContainers();
			void separatrixSorting(int, double, double, double,BaseIntegrator*);
			void adiabaticInvariant();

	    //
	    // Data
	    //
	    // Particles
	    Flow* flow;
	    store_t* store;
			store_t* tracking_store;
	    ParticlesDistribution* part_set;
			TrackingDistribution* tracking_set;
			
			// Integrator
			BaseIntegrator* my_integrator;
	
	    // Physical parameters
	    CurveParameters* A1_;
			CurveParameters* A2_;
			CurveParameters* A3_;
    
	    // Kind of a helper class
	    ParametersEvolution* param_evolution;
	    
	    // Other parameters
			int n_kicks_;
	    int t;
			int dt;
	    int turns;
			double timestep;
	    int n_lost;
	    int intermediate_turns;
	    std::string file_name;
			bool is_phase_portrait;
			bool adiabatic_invariant;
			bool invariant_all_turns;
			int invariant_passes;
	    
	protected:

	public:
        Integrator(ParticlesDistribution*, Flow*);
        ~Integrator();
        void integrate();
        int getTurns();
				TrackingDistribution* getTrackingSet();
        
};

#endif
