#include <distributions/ParticlesSet.h>

ParticlesSet::ParticlesSet()
{
    LOG(PARTICLES SET, Instanciating particle set);
    // Variables
    isSet = false;
    store = NULL;
    
    // We need the number of turns to store the tracking info
    turns = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/iterator/turns");    
    
    // Check the type of distribution
    distr_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/alessandroFlow/initialDistribution", "type");
    if(distr_type != "gaussian" && distr_type != "uniform")
    {
        throw("Invalid type for the initial distribution !");
    }
    distr_dims = static_cast<unsigned short int>(XMLInputParser::getInstance().getAttribute("./flows/alessandroFlow/initialDistribution", "dimensions")); 
    if(distr_dims != 2U && distr_dims != 4U && distr_dims != 6U)
    {
        throw("Invalid dimension for the initial distribution !");
    }
    
    // Read and compute distribution parameters
    readParameters();
    computeDerivedParameters();
}

ParticlesSet::~ParticlesSet()
{
    LOG(PARTICLES SET, Destructor called);
    if(store != NULL)
    {
        delete store;
    }
    if(coherent_parameters != NULL)
    {
        for(int i = 0; i < turns; i++)
        {
            delete [] coherent_parameters[i];
        }
        delete [] coherent_parameters;
    }
    if(tracking_ic != NULL)
    {
        for(int i = 0; i < n_tracking; i++)
        {
            for(int j = 0; j < turns; j++)
            {
                delete [] tracking_ic[i][j];
            }
            delete [] tracking_ic[i];
        }
        delete [] tracking_ic;
    }
}

store_t*
ParticlesSet::getStore() const
{       
	// For performance considerations, it is better to work directly on the container    
	return store;
}

double***
ParticlesSet::getTracking() const
{
    return tracking_ic;
}

double**
ParticlesSet::getCoherent() const
{
    return coherent_parameters;
}

void
ParticlesSet::setCoherent()
{
    LOG(PARTICLES SET, Checking for coherent calculations.);
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking") && is_tracking == true)
    {
        if(XMLInputParser::getInstance().isFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking/coherent"))
        {
            is_coherent = true;
            // Create the containers
            coherent_parameters = new double*[turns];
            for(int i = 0; i < turns; i++)
            {
                coherent_parameters[i] = new double[distr_dims];
                for(int j=0; j < distr_dims; j++)
                {
                    coherent_parameters[i][j] = 0.0;
                }
            }
        }
        else
        {
            is_coherent = false;
        }
    }
    else
    {
        is_coherent = false;
    }
}

void
ParticlesSet::setTracking()
{
    LOG(PARTICLES SET, Checking for tracking.);
    // Check if the tracking is needed
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking"))
    {
        is_tracking = true;
        n_tracking = XMLInputParser::getInstance().howManyFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle");
        
        // Coherence check
        int n_x = XMLInputParser::getInstance().howManyFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/X");
        int n_xp = XMLInputParser::getInstance().howManyFoundFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/XP");
        if(n_x != n_xp || n_x != n_tracking)
        {
            throw("Error setTracking");
        }
        
        // Create the containers
        tracking_ic = new double**[n_tracking]; // Should try to find a better container... like Blitz
        for(int i = 0; i < n_tracking; i++)
        {
            tracking_ic[i] = new double*[turns+1];
            for(int j = 0; j < turns; j++)
            {
                tracking_ic[i][j] = new double[distr_dims];
                for(int v = 0; v < distr_dims; v++)
                {
                    tracking_ic[i][j][v] = 0.0;
                }
            }
        }
        
        // Modify the number of particles, the tracked coordinates will go at the end of the store
        particlesNumber += n_tracking;
        
        // Retrieve the coordinates
        std::vector<std::string> getTextElementsFromPath(std::string);
        std::vector<std::string> x,xp;
        x = XMLInputParser::getInstance().getTextElementsFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/X");
        xp = XMLInputParser::getInstance().getTextElementsFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/XP");
        if(distr_dims == 4)
        {
            std::vector<std::string> y, yp;
            y = XMLInputParser::getInstance().getTextElementsFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/Y");
            yp = XMLInputParser::getInstance().getTextElementsFromPath("./flows/alessandroFlow/initialDistribution/tracking/particle/YP");
        }
        int u = 0;
        for(std::vector<std::string>::iterator i = x.begin(); i < x.end(); i++)
        {
            tracking_ic[u][0][X_INDEX] = s2n<double>(*i);
            u++;
        }
        u = 0;
        for(std::vector<std::string>::iterator i = xp.begin(); i < xp.end(); i++)
        {
            tracking_ic[u][0][XP_INDEX] = s2n<double>(*i);
            u++;
        }
        if(distr_dims == 4)
        {
            int w = 0;
            for(std::vector<std::string>::iterator i = x.begin(); i < x.end(); i++)
            {
                tracking_ic[w][0][Y_INDEX] = s2n<double>(*i);
                u++;
            }
            w = 0;
            for(std::vector<std::string>::iterator i = xp.begin(); i < xp.end(); i++)
            {
                tracking_ic[w][0][YP_INDEX] = s2n<double>(*i);
                u++;
            }
        }
    }
    else
    {
        is_tracking = false;
    }
}

void
ParticlesSet::addTrack2Store()
{
    int index = 0;
    LOG(PARTICLES SET, Adding the tracking IC to the store);
    for(int i = 0; i< n_tracking; i++)
    {
        index = particlesNumber-n_tracking+i;
        (*store)(index,X_INDEX) = tracking_ic[i][0][X_INDEX];
        (*store)(index, XP_INDEX) = tracking_ic[i][0][XP_INDEX];
        if(distr_dims == 4)
        {
            (*store)(index, Y_INDEX) = tracking_ic[i][0][Y_INDEX];
            (*store)(index, YP_INDEX) = tracking_ic[i][0][YP_INDEX];    
        }
    }
}

void
ParticlesSet::setInitialDistribution()
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
	    if(distr_type == "gaussian")
	    {
	        createParticlesGaussian();
	    }
	    else if(distr_type == "uniform")
	    {
	        createParticlesUniform();
	    }
	    // Then add the tracking particles to the store
	    addTrack2Store();
		// Set the flag
		isSet = true;
		// Output some info
		#ifdef DEBUG_FLAG
	        showStoreInfo();
	    #endif
	}
	else
	{
		// Handle the error:w
		
		throw("The distribution is already set ! Want to destroy everything ?");
	}
}

void 
ParticlesSet::readParameters()
{
    LOG(PARTICLES SET, Reading the input parameters...);
    if(distr_type == "gaussian")
    {
        LOG(PARTICLES SET, Reading input for a gaussian distribution...);
        if(distr_dims == 2 || distr_dims == 4 || distr_dims == 6)
        {
            _gaussian_hor_sigma_x = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/sigmas/X");
            _gaussian_hor_sigma_xp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/sigmas/XP");
            _gaussian_hor_centroid_x = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/centroid/X");
            _gaussian_hor_centroid_xp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/centroid/XP");
            _gaussian_hor_density = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/horizontal/density");
        }
        if(distr_dims == 4 || distr_dims == 6)
        {
            _gaussian_ver_sigma_y = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/sigmas/Y");
            _gaussian_ver_sigma_yp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/sigmas/YP");
            _gaussian_ver_centroid_y = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/centroid/Y");
            _gaussian_ver_centroid_yp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/centroid/YP");
            _gaussian_ver_density = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/vertical/density");
        }
        if(distr_dims == 6)
        {
            _gaussian_long_sigma_z = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/sigmas/Z");
            _gaussian_long_sigma_zp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/sigmas/ZP");
            _gaussian_long_centroid_z = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/centroid/Z");
            _gaussian_long_centroid_zp = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/centroid/ZP");
            _gaussian_long_density = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/longitudinal/density");
        }
    }
    else if(distr_type == "uniform")
    {
        LOG(PARTICLES SET, Reading input for a gaussian distribution...);
        if(distr_dims == 2 || distr_dims == 4 || distr_dims == 6)
        {
            _initial_centroid_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/centroid/X");
            _initial_centroid_p_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/centroid/XP");
            _radius_limit_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/horizontal/radius");
            _radial_density_horizontal = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/horizontal/radialDensity");
            _angular_density_horizontal = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/horizontal/angularDensity");
        }
        if(distr_dims == 4 || distr_dims == 6)
        {
            _initial_centroid_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/centroid/Y");
            _initial_centroid_p_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/centroid/YP");
            _radius_limit_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/vertical/radius");
            _radial_density_vertical = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/vertical/radialDensity");
            _angular_density_vertical = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/vertical/angularDensity");
        }
        if(distr_dims == 6)
        {
            _initial_centroid_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/centroid/Z");
            _initial_centroid_p_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/centroid/ZP");
            _radius_limit_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/initialDistribution/longitudinal/radius");
            _radial_density_longitudinal = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/longitudinal/radialDensity");
            _angular_density_longitudinal = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/initialDistribution/longitudinal/angularDensity");
        }
    }
    LOG(PARTICLES SET, Input parameters read.);
}

void 
ParticlesSet::computeDerivedParameters()
{
    LOG(PARTICLES SET, Computing the derived parameters...);
    
    if(distr_type == "gaussian")
    {
        LOG(PARTICLES SET, Derived parameters for a gaussian distribution...);
        if(distr_dims == 2)
        {
            // Total number of particles
            particlesNumber = _gaussian_hor_density * _gaussian_hor_density;
        }
        else if(distr_dims == 4)
        {
            // Total number of particles
            particlesNumber = _gaussian_hor_density * _gaussian_hor_density * _gaussian_ver_density * _gaussian_ver_density;
        }
        else if(distr_dims == 6)
        {
            // Total number of particles
            particlesNumber = _gaussian_hor_density * _gaussian_hor_density * _gaussian_ver_density * _gaussian_ver_density * _gaussian_long_density * _gaussian_long_density;
        }
    }
    else if(distr_type == "uniform")
    {
        LOG(PARTICLES SET, Derived parameters for a uniform distribution...);
        if(distr_dims == 2)
        {
             particlesNumber = (_radial_density_horizontal+1)*(_angular_density_horizontal+1);
             _radial_increment_horizontal = _radius_limit_horizontal / (_radial_density_horizontal+1);
             _angular_increment_horizontal = 2*M_PI/(_angular_density_horizontal+1);
        }
        if(distr_dims == 4)
        {
            particlesNumber = (_radial_density_horizontal+1)*(_angular_density_horizontal+1) * (_radial_density_vertical+1)*(_angular_density_vertical+1);
            _radial_increment_vertical = _radius_limit_vertical / (_radial_density_vertical+1);
            _angular_increment_vertical = 2*M_PI/(_angular_density_vertical+1);
        }
        if(distr_dims == 6)
        {
            particlesNumber = (_radial_density_horizontal+1)*(_angular_density_horizontal+1) * (_radial_density_vertical+1)*(_angular_density_vertical+1) * (_radial_density_longitudinal+1)*(_angular_density_longitudinal+1);
            _radial_increment_longitudinal = _radius_limit_longitudinal / (_radial_density_longitudinal+1);
	        _angular_increment_longitudinal = 2*M_PI/(_angular_density_longitudinal+1);
        }
	}
    LOG(PARTICLES SET, Derived parameters computed.);
}

void
ParticlesSet::createParticlesGaussian()
{
    LOG(PARTICLES SET, Creating the particles --- Gaussian distribution...);
    // Initialize the GSL Random Number generator
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    
	// Create the gaussian distribution
	(*store) = 0.0;
	int loop_index = 0;
	for(unsigned int hxi=0;hxi < _gaussian_hor_density;hxi++)
	{
	    for(unsigned int hxpi=0;hxpi < _gaussian_hor_density;hxpi++)
		{
		    if(distr_dims == 2)
		    {
		        // X
		        (*store)(loop_index,X_INDEX) = _gaussian_hor_centroid_x + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_x));
		        // X'
		        (*store)(loop_index,XP_INDEX) = _gaussian_hor_centroid_xp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_xp));
		        loop_index++;
		    }
		    else if(distr_dims == 4 || distr_dims == 6)
		    {
		        for(unsigned int vyi=0;vyi < _gaussian_ver_density;vyi++)
		        {
		            for(unsigned int vypi=0;vypi < _gaussian_ver_density;vypi++)
		            {
		                if(distr_dims == 4)
		                {
		                    // X
		                    (*store)(loop_index,X_INDEX) = _gaussian_hor_centroid_x + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_x));
		                    // X'
		                    (*store)(loop_index,XP_INDEX) = _gaussian_hor_centroid_xp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_xp));
		                    // Y
		                    (*store)(loop_index,Y_INDEX) = _gaussian_ver_centroid_y + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_y));
		                    // Y'
		                (*store)(loop_index,YP_INDEX) = _gaussian_ver_centroid_yp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_yp));
		                    loop_index++;
		                }
		                else if(distr_dims == 6)
		                {
		                    for(unsigned int vzi=0;vzi < _gaussian_long_density;vzi++)
		                    {
		                        for(unsigned int vzpi=0;vzpi < _gaussian_long_density;vzpi++)
		                        {
		                            // X
		                            (*store)(loop_index,X_INDEX) = _gaussian_hor_centroid_x + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_x));
		                            // X'
		                            (*store)(loop_index,XP_INDEX) = _gaussian_hor_centroid_xp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_xp));
		                            // Y
		                            (*store)(loop_index,Y_INDEX) = _gaussian_ver_centroid_y + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_y));
		                            // Y'
		                            (*store)(loop_index,YP_INDEX) = _gaussian_ver_centroid_yp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_yp));
		                            // Z
		                            (*store)(loop_index,Z_INDEX) = _gaussian_long_centroid_z + gsl_ran_gaussian(r,static_cast<double>(_gaussian_long_sigma_z));
		                            // Z'
		                            (*store)(loop_index,ZP_INDEX) = _gaussian_long_centroid_zp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_long_sigma_zp));
		                            loop_index++;
		                        }
		                    }
		                }
		            }
		        }
		    }
		}
    } 
    LOG(PARTICLES SET, Gaussian distribution created.);
}

void
ParticlesSet::createParticlesUniform()
{
    LOG(PARTICLES SET, Creating the particles);

	// Initialize the store
	(*store) = 0.0;
    
	// Populate the store
	int loop_index = 0;
	for(unsigned int rdhi=0;rdhi <= _radial_density_horizontal;rdhi++)
	{
		for(unsigned int adhi=0;adhi <= _angular_density_horizontal;adhi++)
		{
		    // X
		    (*store)(loop_index,X_INDEX) = _initial_centroid_horizontal + _radial_increment_horizontal*rdhi*cos(_angular_increment_horizontal*adhi);
			// X'
			(*store)(loop_index,XP_INDEX) = _initial_centroid_p_horizontal + _radial_increment_horizontal*rdhi*sin(_angular_increment_horizontal*adhi);
			if(distr_dims == 4 || distr_dims == 6)
			{
			    for(unsigned int rdvi=0;rdvi <= _radial_density_vertical;rdvi++)
			    {
				    for(unsigned int advi=0;advi <= _angular_density_vertical;advi++)
				    {
					    // Y
					    (*store)(loop_index,Y_INDEX) = _initial_centroid_p_vertical + _radial_increment_vertical*rdvi*cos(_angular_increment_vertical*advi); 
					    // Y'
					    (*store)(loop_index,YP_INDEX) = _initial_centroid_p_vertical + _radial_increment_vertical*rdvi*sin(_angular_increment_vertical*advi);
					    loop_index++;
			        }
				}
			}
			else
			{
			    loop_index++;
			}
		}
	}
    LOG(PARTICLES SET, Uniform distribution created);
}

int
ParticlesSet::getParticlesNumber()
{
    return static_cast<int>(particlesNumber);
}

int
ParticlesSet::getDims()
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
ParticlesSet::isTracking()
{
    return is_tracking;
}

int
ParticlesSet::getTrackingNumber()
{
    return n_tracking;
}

bool
ParticlesSet::isCoherent()
{
    return is_coherent;
}

