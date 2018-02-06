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

#ifndef PROC_ISLANDS_H
#define PROC_ISLANDS_H

#include <math.h> 

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>

class ProcIslands
{
	private:
	    // Data
	    store_t* store;
	    ParticlesDistribution* part_set;
	    std::string xpath;
	    int n_processings;
	    std::vector<std::string> planes;
	    std::vector<std::string> boxsize;
	    std::vector<std::string> boxes_x;
	    std::vector<std::string> boxes_y;
			double box_corner_x;
			double box_corner_xp;
			double box_size_x, box_size_y;
	    std::vector<std::vector<std::string> > isl_list;
	    std::ofstream file;
	    std::string file_name;
	    // Flags
	    bool moments;
	    bool derived;
			bool grid_flag;
	    // Methods
	    void readParameters();
			std::vector<double> computeBoxGeometry(int, int, unsigned short int, double, int, int);
			std::vector<double> computeBoxGeometry(int, int, double, double, double, double);
			
	    void processBox(std::vector<double>);
	    int computeLosses();
	    void process();
	    void printHeaders(std::string);

	public:
		ProcIslands(ParticlesDistribution*,std::string);
		~ProcIslands();

};

#endif

