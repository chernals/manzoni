#include <distributions/TrackingDistribution.h>

TrackingDistribution::TrackingDistribution(unsigned int number)
{
    LOG(TRACKING DISTRIBUTION, Instanciating tracking distribution);
    // Variables
    isSet = false;
    store = NULL;	
		particlesNumber = number;
    
    // Check the type of distribution
    distr_dims = static_cast<unsigned short int>(XMLInputParser::getInstance().getAttribute("./flows/symplecticFlow/initialDistribution", "dimensions")); 
    if(distr_dims != 2U)
    {
        throw("Invalid dimension for the initial distribution !");
    }
 
}

TrackingDistribution::~TrackingDistribution()
{
    LOG(TRACKING DISTRIBUTION, Destructor called);
    if(store != NULL)
    {
        delete store;
    }
}

store_t*
TrackingDistribution::getStore() const
{       
	// For performance considerations, it is better to work directly on the container    
	return store;
}

void
TrackingDistribution::setInitialDistribution()
{
    #ifdef DEBUG_FLAG
      char* string = new char[BUFFER_SIZE];
	    sprintf(string, "Creating a store of %ld particles", static_cast<long>(particlesNumber));
	    LOGV(PARTICLES_SET, string);
	    delete string;
		#endif
	
	if(!isSet)
	{
	 	// Create the store object itself
	  store = new store_t(particlesNumber,distr_dims);

		// Set the flag
		isSet = true;
		// Output some info
		#ifdef DEBUG_FLAG
	  	showStoreInfo();
	  #endif
	}
	else
	{
		throw("The distribution is already set !");
	}
}

unsigned int
TrackingDistribution::getParticlesNumber()
{
    return static_cast<unsigned int>(particlesNumber);
}

int
TrackingDistribution::getDims()
{
    if(isSet)
    {
        return static_cast<int>(distr_dims);
    }
    else
    {
        throw("Error; distribution not set.");
    }
}

bool
TrackingDistribution::isTracking()
{
    return false;
}

unsigned int
TrackingDistribution::getTrackingNumber()
{
    return 0;
}

bool
TrackingDistribution::isCoherent()
{
    return false;
}

void
TrackingDistribution::setTracking()
{
}

void
TrackingDistribution::setCoherent()
{
}

double***
TrackingDistribution::getTracking() const
{
    return NULL;
}

double**
TrackingDistribution::getCoherent() const
{
    return NULL;
}
