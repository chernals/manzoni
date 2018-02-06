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

#ifndef PROC_TRACKING_H
#define PROC_TRACKING_H

#include <math.h> 

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <string>
#include <list>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>
#include <root.h>

#define X_PLANE 0
#define Y_PLANE 2

class ProcTracking
{
	private:
        // Misc
        const unsigned short int spacing;
        const unsigned short int precision;
	    std::string xpath;
	    std::string file_name_particles;
	    std::string file_name_coherent;
	    int ics_count; // Number of tracked particles
	    int turns;
	    int dims;
	    int fft_window;
	    int fft_step;
	    int fft_out_step;
	    // Data
  	    store_t* store;
	    ParticlesDistribution* part_set;
	    double*** tracking_ic;
	    double** coherent_data;
	    double* coherent_tune_data;
	    double** tracking_tune_data;
	    // Flags
	    bool withParticles;
	    bool withCoherent;
	    bool withPlots;
	    bool withFFT;
	    bool withTune;
	    
	    // Methods
	    void readParameters();
	    void write2File();
	    void drawPlots(int);
	    void computeFFT(int);
	    void plotFFT(int, double*, std::string);
	    void plotComputedTune(double*, std::string);
	    void process();
	    double tuneFunctionA(double,double,double);

	public:
		ProcTracking(ParticlesDistribution*,std::string);
		~ProcTracking();

};

#endif

