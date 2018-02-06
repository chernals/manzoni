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

#ifndef PROC_DISTR_HISTO_H
#define PROC_DISTR_HISTO_H

#include <math.h>

#include <string>
#include <iomanip>
#include <list>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>
#include <utils/TypeConversions.h>
#include <root.h>

class ProcDistrHisto
{
	private:
	    //
	    // Data
	    //
	    store_t* store;
	    ParticlesDistribution* part_set;
	    std::string xpath;
	    std::string file;
	    std::string title;
        std::vector<std::string> projections;
        std::list<TH2D*> hist_list;
        double box_size;
				bool polar;
				bool negative;
        TCanvas* canvas;
        unsigned short int hist_counter;
		int histo_counter;
		std::vector<double> tunes;
        //
        // Methods
        //
        void readParameters();

	public:
		ProcDistrHisto(ParticlesDistribution*, std::string,int,std::vector<double>);
		~ProcDistrHisto();
		void createHistogram(std::string, std::string, unsigned short int, unsigned short int);
		void draw();
		void save();
};

#endif

