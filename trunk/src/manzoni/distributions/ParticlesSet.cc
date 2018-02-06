#include <distributions/ParticlesSet.h>

ParticlesSet::ParticlesSet()
{
    LOG(PARTICLES SET, Instanciating particle set);
    // Variables
    isSet = false;
    store = NULL;
    
    // We need the number of turns to store the tracking info
    LOG(PARTICLES SET, Reading turns);
   // turns = XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/iterator/turns");
    bool isFound = false;    
    turns = s2n<int>(XMLInputParser::getInstance().getFirstTextElementFromPath("./flows/henonFlow/iterator/turns", isFound));
    
    // Check the type of distribution
    LOG(PARTICLES SET, Checking distribution type);
    distr_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/initialDistribution", "type");
    if(distr_type != "gaussian" && distr_type != "uniform" && distr_type != "file" && distr_type != "hollow")
    {
        throw("Invalid type for the initial distribution !");
    }
    distr_dims = static_cast<unsigned short int>(XMLInputParser::getInstance().getAttribute("./flows/henonFlow/initialDistribution", "dimensions")); 
    if(distr_dims != 2U && distr_dims != 4U && distr_dims != 6U)
    {
        throw("Invalid dimension for the initial distribution !");
    }
    
    // Read and compute distribution parameters
    LOG(PARTICLES SET, Reading parameters);
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
      //  for(int i = 0; i < turns; i++)
            //delete [] coherent_parameters[i];
       
        
        //delete [] coherent_parameters;
    }
    
    if(tracking_ic != NULL)
    {
        //for(unsigned int i = 0; i < n_tracking; i++)
        //{
          //  for(int j = 0; j < turns; j++)
          //      delete [] tracking_ic[i][j];
          //  delete [] tracking_ic[i];
        //}
  //      delete [] tracking_ic;
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
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/initialDistribution/tracking") && is_tracking == true)
    {
        if(XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/initialDistribution/tracking/coherent"))
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
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/initialDistribution/tracking"))
    {
        is_tracking = true;
        n_tracking = XMLInputParser::getInstance().howManyFoundFromPath("./flows/henonFlow/initialDistribution/tracking/particle");
        
        // Coherence check
        unsigned int n_x = XMLInputParser::getInstance().howManyFoundFromPath("./flows/henonFlow/initialDistribution/tracking/particle/X");
        unsigned int n_xp = XMLInputParser::getInstance().howManyFoundFromPath("./flows/henonFlow/initialDistribution/tracking/particle/XP");
        if(n_x != n_xp || n_x != n_tracking)
        {
            throw("Error setTracking");
        }
        
        // Create the containers
        tracking_ic = new double**[n_tracking]; // Should try to find a better container... like Blitz
        for(unsigned int i = 0; i < n_tracking; i++)
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
        x = XMLInputParser::getInstance().getTextElementsFromPath("./flows/henonFlow/initialDistribution/tracking/particle/X");
        xp = XMLInputParser::getInstance().getTextElementsFromPath("./flows/henonFlow/initialDistribution/tracking/particle/XP");
        if(distr_dims == 4)
        {
            std::vector<std::string> y, yp;
            y = XMLInputParser::getInstance().getTextElementsFromPath("./flows/henonFlow/initialDistribution/tracking/particle/Y");
            yp = XMLInputParser::getInstance().getTextElementsFromPath("./flows/henonFlow/initialDistribution/tracking/particle/YP");
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
    unsigned int index = 0;
    std::cout << n_tracking << std::endl;
    LOG(PARTICLES SET, Adding the tracking IC to the store);
    for(unsigned int i = 0; i< n_tracking; i++)
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
      else if(distr_type == "hollow")
      {
        createParticlesHollow();
      }
			else if(distr_type == "file")
			{
				createParticlesFile();
			}
	    // Then add the tracking particles to the store
	    if(is_tracking)
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
		// Handle the error	
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
            _gaussian_hor_sigma_x = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/sigmas/X");
            _gaussian_hor_sigma_xp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/sigmas/XP");
            _gaussian_hor_centroid_x = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/centroid/X");
            _gaussian_hor_centroid_xp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/centroid/XP");
            _gaussian_hor_density = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/horizontal/density"));
        }
        if(distr_dims == 4 || distr_dims == 6)
        {
            _gaussian_ver_sigma_y = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/sigmas/Y");
            _gaussian_ver_sigma_yp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/sigmas/YP");
            _gaussian_ver_centroid_y = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/Y");
            _gaussian_ver_centroid_yp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/YP");
            _gaussian_ver_density = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/vertical/density"));
        }
        if(distr_dims == 6)
        {
            _gaussian_long_sigma_z = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/sigmas/Z");
            _gaussian_long_sigma_zp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/sigmas/ZP");
            _gaussian_long_centroid_z = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/centroid/Z");
            _gaussian_long_centroid_zp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/centroid/ZP");
            _gaussian_long_density = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/longitudinal/density"));
        }
    }
    else if(distr_type == "hollow")
    {
      LOG(PARTICLES SET, Reading input for a hollow distribution...);
      if(distr_dims == 2 || distr_dims == 4)
      {
        _gaussian_hor_centroid_offset = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/offset");
        _gaussian_hor_sigma_x = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/sigma");
        _radial_density_horizontal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/horizontal/radialDensity"));
        _angular_density_horizontal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/horizontal/angularDensity"));
      }
      if(distr_dims == 4)
      {
        _gaussian_ver_sigma_y = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/sigmas/Y");
        _gaussian_ver_sigma_yp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/sigmas/YP");
        _gaussian_ver_centroid_y = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/Y");
        _gaussian_ver_centroid_yp = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/YP");
        _gaussian_ver_density = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/vertical/density"));
      }
    }
    else if(distr_type == "uniform")
    {
        LOG(PARTICLES SET, Reading input for a uniform distribution...);
        if(distr_dims == 2 || distr_dims == 4 || distr_dims == 6)
        {
            _initial_centroid_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/centroid/X");
            _initial_centroid_p_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/centroid/XP");
            _radius_limit_horizontal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/horizontal/radius");
            _radial_density_horizontal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/horizontal/radialDensity"));
            _angular_density_horizontal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/horizontal/angularDensity"));
        }
        if(distr_dims == 4 || distr_dims == 6)
        {
            _initial_centroid_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/Y");
            _initial_centroid_p_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/centroid/YP");
            _radius_limit_vertical = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/vertical/radius");
            _radial_density_vertical = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/vertical/radialDensity"));
            _angular_density_vertical = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/vertical/angularDensity"));
        }
        if(distr_dims == 6)
        {
            _initial_centroid_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/centroid/Z");
            _initial_centroid_p_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/centroid/ZP");
            _radius_limit_longitudinal = XMLInputParser::getInstance().getFirstTextd("./flows/henonFlow/initialDistribution/longitudinal/radius");
            _radial_density_longitudinal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/longitudinal/radialDensity"));
            _angular_density_longitudinal = static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/initialDistribution/longitudinal/angularDensity"));
        }
    }
		else if(distr_type == "file")
		{
			_distr_file = XMLInputParser::getInstance().getFirstTextElementFromPath("./flows/henonFlow/initialDistribution/file");
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
    if(distr_type == "hollow")
    {
        LOG(PARTICLES SET, Derived parameters for a hollow distribution...);
        if(distr_dims == 2)
        {
            // Total number of particles
            particlesNumber = _radial_density_horizontal * (_angular_density_horizontal+1);
            _angular_increment_horizontal = 2*M_PI/(_angular_density_horizontal+1);
        }
        else if(distr_dims == 4)
        {
            // Total number of particles
            particlesNumber = _radial_density_horizontal * (_angular_density_horizontal+1) * _gaussian_ver_density * _gaussian_ver_density;
        _angular_increment_horizontal = 2*M_PI/(_angular_density_horizontal+1);
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
            particlesNumber = (_radial_density_horizontal)*(_angular_density_horizontal) * (_radial_density_vertical)*(_angular_density_vertical);
            //_radial_increment_horizontal = _radius_limit_horizontal / (_radial_density_horizontal+1);
            //_angular_increment_horizontal = 2*M_PI/(_angular_density_horizontal+1);           
 						//_radial_increment_vertical = _radius_limit_vertical / (_radial_density_vertical+1);
            //_angular_increment_vertical = 2*M_PI/(_angular_density_vertical+1);
        }
        if(distr_dims == 6)
        {
            particlesNumber = (_radial_density_horizontal+1)*(_angular_density_horizontal+1) * (_radial_density_vertical+1)*(_angular_density_vertical+1) * (_radial_density_longitudinal+1)*(_angular_density_longitudinal+1);
            _radial_increment_longitudinal = _radius_limit_longitudinal / (_radial_density_longitudinal+1);
	        _angular_increment_longitudinal = 2*M_PI/(_angular_density_longitudinal+1);
        }
	}
	else if(distr_type == "file")
	{
		std::ifstream f(_distr_file.c_str());
		std::string line;
    while (std::getline(f, line))
	        ++particlesNumber;
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
	unsigned int loop_index = 0;
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
		                
		            }
		        }
		    }
		}
    } 
    LOG(PARTICLES SET, Gaussian distribution created.);
}

void
ParticlesSet::createParticlesHollow()
{
    LOG(PARTICLES SET, Creating the particles --- Hollow distribution...);
    // Initialize the GSL Random Number generator
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    
	// Create the hollow distribution
	(*store) = 0.0;
	unsigned int loop_index = 0;
	for(unsigned int rdhi=0;rdhi < _radial_density_horizontal;rdhi++)
	{
		for(unsigned int adhi=0;adhi <= _angular_density_horizontal;adhi++)
		{
		        for(unsigned int vyi=0;vyi < _gaussian_ver_density;vyi++)
		        {
		            for(unsigned int vypi=0;vypi < _gaussian_ver_density;vypi++)
		            {
                  
		                 (*store)(loop_index,X_INDEX) = 1.0;
                     // X'
                     (*store)(loop_index,XP_INDEX) = 1.0;
                     // Y
                     (*store)(loop_index,Y_INDEX) =1.0;
                     // Y'
                     (*store)(loop_index,YP_INDEX) = 1.0;
                     
                     while( (*store)(loop_index,X_INDEX) >= 1.0 || 
                        (*store)(loop_index,XP_INDEX) >= 1.0 ||
                           (*store)(loop_index,Y_INDEX) >= 1.0 ||
                              (*store)(loop_index,YP_INDEX) >= 1.0   )
                     {
                       double tmp = _gaussian_hor_centroid_offset + gsl_ran_gaussian(r,static_cast<double>(_gaussian_hor_sigma_x));
                    // X
		                    (*store)(loop_index,X_INDEX) = tmp*cos(_angular_increment_horizontal*adhi);
		                    // X'
		                    (*store)(loop_index,XP_INDEX) = tmp*sin(_angular_increment_horizontal*adhi);
		                    // Y
		                    (*store)(loop_index,Y_INDEX) = _gaussian_ver_centroid_y + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_y));
		                    // Y'
		                (*store)(loop_index,YP_INDEX) = _gaussian_ver_centroid_yp + gsl_ran_gaussian(r,static_cast<double>(_gaussian_ver_sigma_yp));
                  }
		                    loop_index++;

		            }
		    }
		}
    } 
    LOG(PARTICLES SET, Gaussian distribution created.);
}

//void
//ParticlesSet::createParticlesUniform()
//{
//    LOG(PARTICLES SET, Creating the particles --- Uniform distribution...);
//
//	// Initialize the store
//	(*store) = 0.0;
//    
//	// Populate the store
//	int loop_index = 0;
//	for(unsigned int rdhi=0;rdhi <= _radial_density_horizontal;rdhi++)
//	{
//		for(unsigned int adhi=0;adhi <= _angular_density_horizontal;adhi++)
//		{
//		 	if(distr_dims == 4 || distr_dims == 6)
//			{
//			    for(unsigned int rdvi=0;rdvi <= _radial_density_vertical;rdvi++)
//			    {
//				    for(unsigned int advi=0;advi <= _angular_density_vertical;advi++)
//				    {
//					   // X
//					    (*store)(loop_index,X_INDEX) = _initial_centroid_horizontal + _radial_increment_horizontal*rdhi*cos(_angular_increment_horizontal*adhi);
//						// X'
//						(*store)(loop_index,XP_INDEX) = _initial_centroid_p_horizontal + _radial_increment_horizontal*rdhi*sin(_angular_increment_horizontal*adhi);
//					
//					    // Y
//					    (*store)(loop_index,Y_INDEX) = _initial_centroid_vertical + _radial_increment_vertical*rdvi*cos(_angular_increment_vertical*advi); 
//					    // Y'
//					    (*store)(loop_index,YP_INDEX) = _initial_centroid_p_vertical + _radial_increment_vertical*rdvi*sin(_angular_increment_vertical*advi);
//					    loop_index++;
//			        }
//				}
//			}
//			else
//			{
//				   // X
//				    (*store)(loop_index,X_INDEX) = _initial_centroid_horizontal + _radial_increment_horizontal*rdhi*cos(_angular_increment_horizontal*adhi);
//					// X'
//					(*store)(loop_index,XP_INDEX) = _initial_centroid_p_horizontal + _radial_increment_horizontal*rdhi*sin(_angular_increment_horizontal*adhi);
//				
//			    loop_index++;
//			}
//		}
//	}
//    LOG(PARTICLES SET, Uniform distribution created);
//}

void
ParticlesSet::createParticlesUniform()
{
  LOG(PARTICLES SET, Creating the particles --- Uniform distribution...); 
  
  // Initialize the GSL Random Number generator
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  
  // Initialize the store
  (*store) = 0.0;
  
  // Populate the store
  unsigned int n_max = 0;
  if(distr_dims == 2)
  {
    n_max = _radial_density_horizontal * _angular_density_horizontal;
  }
  else if(distr_dims == 4)
  {
    n_max = _radial_density_horizontal * _angular_density_horizontal * _radial_density_vertical * _angular_density_vertical;    
  }
  
	for(unsigned int loop_index=0; loop_index < n_max;loop_index++)
	{
    if(distr_dims == 2 || distr_dims == 4)
    {
      double tmp_x = gsl_rng_uniform(r);
      double tmp_xp = gsl_rng_uniform(r);
      double angle = 2.0 * M_PI * tmp_xp;
      
      // X
			(*store)(loop_index,X_INDEX) = _radius_limit_horizontal * sqrt(tmp_x) * sin(angle) + _initial_centroid_horizontal;
      // X'
			(*store)(loop_index,XP_INDEX) = _radius_limit_horizontal * sqrt(tmp_x) * cos(angle) + _initial_centroid_p_horizontal;
  
    }
		
    if(distr_dims == 4)
		{
      double tmp_y = gsl_rng_uniform(r);
      double tmp_yp = gsl_rng_uniform(r);
      double angle = 2.0 * M_PI * tmp_yp;
    
			// Y
			(*store)(loop_index,Y_INDEX) = _radius_limit_vertical * sqrt(tmp_y) * sin(angle) + _initial_centroid_vertical;
			// Y'
			(*store)(loop_index,YP_INDEX) = _radius_limit_vertical * sqrt(tmp_y) * cos(angle) + _initial_centroid_p_vertical;
		}
	}
  LOG(PARTICLES SET, Uniform distribution created);
}

void
ParticlesSet::createParticlesFile()
{
	LOG(PARTICLES SET, Reading the particles from a file...);
	// Initialize the store
	(*store) = 0.0;
	
	std::ifstream file;
	file.open(_distr_file.c_str());
	
	if(file.is_open())
	{
		LOG(PARTICLES SET, File is open.);
		unsigned int index = 0;
		while(!file.eof())
		{
			if(distr_dims == 2)
			{
				file >> (*store)(index,X_INDEX) >> (*store)(index,XP_INDEX);
			}
			else if(distr_dims == 4)
			{
				std::string line;
				getline(file, line );
				std::stringstream ss(line);
				std::string field;
				unsigned int i = 0;
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
	} // If file is open
	else
	{
		ERROR(PARTICLES SET, Error reading the file !);
		throw;
	}
}

unsigned int
ParticlesSet::getParticlesNumber()
{
    return particlesNumber;
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

unsigned int
ParticlesSet::getTrackingNumber()
{
    return n_tracking;
}

bool
ParticlesSet::isCoherent()
{
    return is_coherent;
}

