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

#include <flows/SymplecticFlow.h>

SymplecticFlow::SymplecticFlow()
{
}

SymplecticFlow::~SymplecticFlow()
{
}

void
SymplecticFlow::run()
{
  // Creation of the initial distribution
	LOG(SYMPLECTIC_FLOW, Creating the initial distribution);
	particles_set = new AADistribution();
	particles_set->setInitialDistribution();
	
	// Preprocess
    LOG(SYMPLECTIC_FLOW, Preprocessing ?);
	process("./flows/symplecticFlow/dataProcessing/initial/");

   LOG(SYMPLECTIC_FLOW, Create and start the integrator);
   integrator = new Integrator(particles_set, dynamic_cast<Flow*>(this));
   switch(particles_set->getDims())
    {
        case TWO_DIMENSIONS:
            integrator->integrate();
            break;
        default:
            throw("Invalid dimension !");
            break;
    }

		// Phase portrait ?
		if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "true")
		{
			particles_set = integrator->getTrackingSet();
		}
		
		if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dump") == "true")
		{
			// Also output the whole bunch of particles into a file
			int number_particles = particles_set->getParticlesNumber();
			store_t* end_store = particles_set->getStore();
			std::ofstream file_end;
			file_end.open(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dumpfile").c_str());
			for(int m = 0; m < number_particles; m++)
			{
				file_end << std::setw(24) << std::setprecision(15) << std::scientific 
	  	        	 << (*end_store)(m, ANGLE_INDEX)
    	           << std::setw(24) << std::setprecision(15) << std::scientific 
								 << (*end_store)(m, ACTION_INDEX)
								 << std::endl;
			}
		}

    // Preprocess
    LOG(SYMPLECTIC_FLOW, Postprocessing ?);
		process("./flows/symplecticFlow/dataProcessing/final/");

    LOG(SYMPLECTIC_FLOW, End of run.);
}
