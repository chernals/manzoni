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

#ifndef ITERATOR_H
#define ITERATOR_H

#include <time.h>
#include <math.h>                        

#include <iostream> 
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
                               
#include <definitions.h>
#include <distributions/ParticlesDistribution.h>
#include <ParametersEvolution.h>
#include <types/types.h>
#include <utils/Logger.h>
#include <dataProcessing/procGraph.h>

class Flow;

class Iterator
{
	private:
	    //
	    // Methods
	    //
	    void readParameters();
	    inline double evaluateVariable(CurveParameters&);
	    void createContainers();
	    inline void phaseRotation(double, int); 
	    void readKicksTypes();
			void iterateParticle4(std::vector<double>&, std::string, double, double, double, double,
                                     		 	 double, double, double, double);
	    //
	    // Data
	    //
	    // Particles
	    Flow* flow;
	    store_t* store;
	    ParticlesDistribution* part_set;
	    
	    // Tracking
	    double*** tracking_ic;
	    bool is_tracked;
	    unsigned int n_tracked;
	    
	    // Coherent motion
	    double** coherent_parameters;
	    bool is_coherent;
	    int n_lost;
	    	    
	    // Kicks
	    int n_kicks_;
	    std::vector<std::string> kick_types; // Can be S, O, SO or U4 (future addition possible)
	    
	    // Physical parameters
	    std::vector<double> normalizations;
	    CurveParameters* hor_tune_;
	    CurveParameters* ver_tune_;
	    CurveParameters* long_tune;
	    CurveParameters* kappa_;
	    std::vector<double> chi;
	    CurveParameters* hor_chromaticity;
	    CurveParameters* ver_chromaticity;
	    CurveParameters* dispersion;
	    CurveParameters* dispersion_prime;
	    #ifdef DAMPER_FLAG
	    	CurveParameters* damper;
	    #endif
	    
	    // Evolution tables for the plot of the parameters
	    double** table_tune_h;
	    double** table_tune_v;
	    double** table_tune_l;
	    double** table_kappa;
	    #ifdef DAMPER_FLAG
	    	double** table_damper;
	    #endif
	    
	    // Kind of a helper class
	    ParametersEvolution* param_evolution;
	    
	    // Other parameters
	    int turns;
			int dimensions;
	    int intermediate_turns;
	    double horizontal_phase_modification;
	    double vertical_phase_modification;
	    
	    // State variables
	    int t;
	    #ifdef DAMPER_FLAG
	        bool* isDamperSet;
	    #endif
	    
	    // Misc
	    ProcGraph* graph_;
	    bool with_plots_;
	    std::string file_name;
			bool adiabatic_invariant;
			#ifdef INVARIANT_FLAG
				std::ofstream file_invariant;
			#endif
	    
	protected:

	public:
        Iterator(ParticlesDistribution*, Flow*);
        ~Iterator();
        void iterate2d();
        void iterate4d();
        void iterate6d();
        int getTurns();
        
};

#endif
