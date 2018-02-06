#ifndef PROC_ISLANDS_H
#define PROC_ISLANDS_H

#include <math.h> 

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>

#include <distributions/ParticlesDistribution.h>
#include <XMLInputParser.h>

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
	    std::vector<std::vector<std::string> > isl_list;
	    std::ofstream file;
	    std::string file_name;
	    // Flags
	    bool moments;
	    bool derived;
	    // Methods
	    void readParameters();
	    void processBox(int, int, unsigned short int, double, int, int);
	    int computeLosses();
	    void process();
	    void printHeaders(std::string);

	public:
		ProcIslands(ParticlesDistribution*,std::string);
		~ProcIslands();

};

#endif

