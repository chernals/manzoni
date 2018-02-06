#include <SymplecticFlow.h>

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
	LOG(ALESSANDRO_FLOW, Preprocessing ?);
	process("./flows/symplecticFlow/dataProcessing/initial/");

   LOG(ALESSANDRO_FLOW, Create and start the integrator);
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
		
		// Also output the whole bunch of particles into a file
		int number_particles = particles_set->getParticlesNumber();
		store_t* end_store = particles_set->getStore();
		std::ofstream file_end;
		file_end.open(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "dump").c_str());
		for(int m = 0; m < number_particles; m++)
		{
			file_end << std::setw(SPACING) << std::setprecision(5) << std::scientific 
	          	 << (*end_store)(m, ANGLE_INDEX)
               << std::setw(SPACING) << std::setprecision(5) << std::scientific 
							 << (*end_store)(m, ACTION_INDEX)
							 << std::endl;
		}

    // Preprocess
		LOG(ALESSANDRO_FLOW, Postprocessing ?);
		process("./flows/symplecticFlow/dataProcessing/final/");

    LOG(SYMPLECTIC_FLOW, End of run.);
}
