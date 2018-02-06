#include <AlessandroFlow.h>

AlessandroFlow::AlessandroFlow()
{
    LOG(ALESSANDRO_FLOW, Instanciating Alessandro flow);
}

AlessandroFlow::~AlessandroFlow()
{
    LOG(ALESSANDRO FLOW, Destructor called);
}

void
AlessandroFlow::run()
{
    LOG(ALESSANDRO_FLOW, Starting the flow);
        
	// Creation of the initial distribution
	LOG(ALESSANDRO_FLOW, Creating the initial distribution);
	particles_set = new ParticlesSet();
	particles_set->setTracking();
	particles_set->setCoherent();
	particles_set->setInitialDistribution();

	// Preprocess
	LOG(ALESSANDRO_FLOW, Preprocessing ?);
	process("./flows/alessandroFlow/dataProcessing/initial/");

    // Creationg of the iterator and start of the iterations
    LOG(ALESSANDRO_FLOW, Create and start the iterator);
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
    LOG(ALESSANDRO_FLOW, Postprocessing ?);
    process("./flows/alessandroFlow/dataProcessing/final/");
    
    LOG(ALESSANDRO_FLOW, End of run.);
}
