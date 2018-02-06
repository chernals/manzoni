#include <distributions/AADistribution.h>

AADistribution::AADistribution()
{
    LOG(PARTICLES SET, Instanciating particle set);
    // Variables
    isSet = false;
    store = NULL;
    
    // Check the type of distribution
    distr_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/initialDistribution", "type");
    if(distr_type != "uniform" && distr_type != "gaussian" && distr_type != "file")
    {
        throw("Invalid type for the initial distribution !");
    }
    distr_dims = static_cast<unsigned short int>(XMLInputParser::getInstance().getAttribute("./flows/symplecticFlow/initialDistribution", "dimensions")); 
    if(distr_dims != 2U)
    {
        throw("Invalid dimension for the initial distribution !");
    }
    
    // Read and compute distribution parameters
    readParameters();
}

AADistribution::~AADistribution()
{
    LOG(PARTICLES SET, Destructor called);
    if(store != NULL)
    {
        delete store;
    }
}

store_t*
AADistribution::getStore() const
{       
	// For performance considerations, it is better to work directly on the container    
	return store;
}

void
AADistribution::setInitialDistribution()
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
		if(distr_type == "uniform")
	  {
	  	createParticlesUniform();
	  }
		else if(distr_type == "gaussian")
		{
			createParticlesGaussian ();
		}
		else if(distr_type == "file")
		{
			createParticlesFile();
		}
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

void 
AADistribution::readParameters()
{
    LOG(PARTICLES SET, Reading the input parameters...);
    if(distr_type == "uniform" || distr_type == "gaussian")
    {
        LOG(PARTICLES SET, Reading input for a uniform distribution...);
        if(distr_dims == 2)
        {
            _uniform_radius = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/initialDistribution/radius");
            _uniform_radial_density = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/initialDistribution/radialDensity");
						_radial_offset = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/initialDistribution/offset");
            _uniform_angular_density = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/initialDistribution/angularDensity");
						if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "false")
						{
							_angular_extension = true;
						}
						else if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "true")
						{
							_angular_extension= false;
  						//_uniform_angular_density = 0;
						}	
						else
						{
							throw("Not a valid input: phasePortrait: true or false");
						}
						particlesNumber = (_uniform_radial_density+1)*(_uniform_angular_density+1);
			      _radial_increment_horizontal = _uniform_radius / (_uniform_radial_density+1);
			      _angular_increment_horizontal = 2*M_PI/(_uniform_angular_density+1);
				}
    }
		else if(distr_type == "file")
		{
			_distr_file = XMLInputParser::getInstance().getFirstTextElementFromPath("./flows/symplecticFlow/initialDistribution/file");
			std::ifstream f(_distr_file.c_str());
			std::string line;
	    while (std::getline(f, line))
		  	++particlesNumber;
		}
    LOG(PARTICLES SET, Input parameters read.);
}

void
AADistribution::createParticlesUniform()
{
    LOG(PARTICLES SET, Creating the particles);

	// Initialize the store
	(*store) = 0.0;
    
	// Populate the store
	int loop_index = 0; 
	for(unsigned int rdhi=0;rdhi <= _uniform_radial_density;rdhi++)
	{
		for(unsigned int adhi=0;adhi <= _uniform_angular_density;adhi++)
		{
		    (*store)(loop_index,ACTION_INDEX) = _radial_increment_horizontal  * rdhi + _radial_offset;
				if(_angular_extension == true)
				{
					(*store)(loop_index,ANGLE_INDEX)  = -M_PI + _angular_increment_horizontal * adhi;
		    }
				else
				{
					(*store)(loop_index,ANGLE_INDEX)  = _angular_increment_horizontal * adhi;
				}
				loop_index++;
			}
	}
    LOG(PARTICLES SET, Uniform distribution created);
}

void
AADistribution::createParticlesGaussian()
{
    LOG(PARTICLES SET, Creating the particles);

	// Initialize the store
	(*store) = 0.0;
    
	// Populate the store
	int loop_index = 0; 
	
	// For GSL
	gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  
	for(unsigned int rdhi=0;rdhi <= _uniform_radial_density;rdhi++)
	{
		for(unsigned int adhi=0;adhi <= _uniform_angular_density;adhi++)
		{
			 	(*store)(loop_index,ACTION_INDEX) =  gsl_ran_gaussian(r,static_cast<double>(_uniform_radius)) + _radial_offset;
				if(_angular_extension == true)
				{
					(*store)(loop_index,ANGLE_INDEX)  = -M_PI + (2*M_PI/(_uniform_angular_density+1)) * adhi;
		    }
				else
				{
					(*store)(loop_index,ANGLE_INDEX)  = 0;
				}
				loop_index++;
			}
	}
    LOG(PARTICLES SET, Uniform distribution created);
}

void
AADistribution::createParticlesFile()
{
	LOG(PARTICLES SET, Reading the particles from a file...);
	// Initialize the store
	(*store) = 0.0;
	
	std::ifstream file;
	file.open(_distr_file.c_str());
	
	if(file.is_open())
	{
		LOG(PARTICLES SET, File is open.);
		int index = 0;
		while(!file.eof())
		{
			if(distr_dims == 2)
			{
				std::string line;
				getline(file, line );
				std::stringstream ss(line);
				std::string field;
				int i = 0;
				while (getline(ss, field, ' ' ))
				{
					std::stringstream fs(field);
					double f = 0.0;
					fs >> f;
					if(i < 2)
						(*store)(index,i) = f;				
					i++;
				}
			}
			else if(distr_dims == 4)
			{
				std::string line;
				getline(file, line );
				std::stringstream ss(line);
				std::string field;
				int i = 0;
				while (getline(ss, field, ' ' ))
				{
					std::stringstream fs(field);
					double f = 0.0;
					fs >> f;
					(*store)(index,i) = f;				
					i++;
				}
			}
			index++;
		} // While !eof
		LOG(PARTICLES SET, Finished.);
		file.close();
	} // If file is open
	else
	{
		ERROR(PARTICLES SET, Error reading the file !);
		throw;
	}
}

int
AADistribution::getParticlesNumber()
{
    return static_cast<int>(particlesNumber);
}

int
AADistribution::getDims()
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
AADistribution::isTracking()
{
    return false;
}

int
AADistribution::getTrackingNumber()
{
    return 0;
}

bool
AADistribution::isCoherent()
{
    return false;
}

void
AADistribution::setTracking()
{
}

void
AADistribution::setCoherent()
{
}

double***
AADistribution::getTracking() const
{
    return NULL;
}

double**
AADistribution::getCoherent() const
{
    return NULL;
}
