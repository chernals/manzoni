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

#include <flows/HenonFlow.h>

HenonFlow::HenonFlow()
{
}

HenonFlow::~HenonFlow()
{
}

void
HenonFlow::run()
{
	// Creation of the initial distribution
  LOG(HENON_FLOW, Creating the initial distribution);
	particles_set = new ParticlesSet();
	particles_set->setTracking();
	particles_set->setCoherent();
	particles_set->setInitialDistribution();

	// Preprocess
  LOG(HENON_FLOW, Preprocessing ?);
	process("./flows/henonFlow/dataProcessing/initial/");

    // Creationg of the iterator and start of the iterations
    LOG(HENON_FLOW, Create and start the iterator);
    iterator = new Iterator(particles_set, dynamic_cast<Flow*>(this));
    switch(particles_set->getDims())
    {
        case TWO_DIMENSIONS:
            iterator->iterate2d();
            break;
        case FOUR_DIMENSIONS:
            iterator->iterate4d();
            break;
        case SIX_DIMENSIONS:
            iterator->iterate6d();
            break;
        default:
            throw("Invalid dimension !");
            break;
    }
    
    // Postprocess
    LOG(HENON_FLOW, Postprocessing ?);
    process("./flows/henonFlow/dataProcessing/final/");
    
    LOG(HENON_FLOW, End of run.);
}
