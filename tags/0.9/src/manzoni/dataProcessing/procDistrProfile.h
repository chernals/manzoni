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

#ifndef PROC_DISTR_PROFILE_H
#define PROC_DISTR_PROFILE_H

#include <string>
#include <list>
#include <vector>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>
#include <root.h>

class ProcDistrProfile
{
	private:
	    //
	    // Data
	    //
	    store_t* store;
	    ParticlesDistribution* _part_set;
	    const std::string xpath_;
	    std::string _file;
	    std::string _title;
	    unsigned short int _x_divs;
        unsigned short int _y_divs;
        std::vector<std::string> profiles_;
        std::list<TH1D*> _profiles_list;
        double box_size_;
        TCanvas* _canvas;
        unsigned short int _profile_counter;
        //
        // Methods
        //
        void readParameters();

	public:
		ProcDistrProfile(ParticlesDistribution*, std::string);
		~ProcDistrProfile();
		void createProfile(std::string, unsigned short int);
		void draw();
		void save();
};

#endif

