#ifndef FLOW_H
#define FLOW_H

#include <math.h>                        

#include <iostream> 
               
#include <definitions.h>
#include <ParticlesSet.h>
#include <Iterator.h>
#include <XMLInputParser.h>
#include <dataProcessing/procDistrHisto.h>
#include <dataProcessing/procDistrProfile.h>
#include <dataProcessing/procIslands.h>
#include <dataProcessing/procTracking.h>
#include <utils/Logger.h>

class Flow
{
	private:
		int histo_counter;
		std::vector<double> tunes;
	    
	protected:              
		//
	    // Data
	    //
	    ParticlesSet* particles_set;
	    Iterator* iterator;

	public:
        Flow();
        virtual ~Flow();
        virtual void run();
        virtual void process(std::string); 
		virtual void setTune(double);
};

#endif
