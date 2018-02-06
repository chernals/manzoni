#include <Integrator.h>
#include <SymplecticFlow.h>

Integrator::Integrator(ParticlesDistribution* set, Flow* f)
{
	LOG(INTEGRATOR, Instanciating integrator);
	
	// Particle store
	this->flow = f;
	this->part_set = set;
	this->store = part_set->getStore();
	n_kicks_ = 1;
	
	// Read the input parameters in the XML file
	std::string prefix="./flows/symplecticFlow/integrator";
	param_evolution = new ParametersEvolution(prefix+"/turns",prefix+"/hamiltonian");
	readParameters();
	
	// Set state variables
	t = 0;
	 
	// Type of integrator
	std::string integrator_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator/hamiltonian0", "type");
	if(integrator_type == "pendulum")
	{
		my_integrator = new Pendulum(timestep, XMLInputParser::getInstance().getAttribute("./flows/symplecticFlow/integrator", "order"));
	}
	else
	{
		throw("Invalid integrator type");
	}
	
	// Phase portrait ?
	if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "true")
	{
		is_phase_portrait = true;
		
		// Tracking store
		tracking_set = new TrackingDistribution(static_cast<unsigned int>(XMLInputParser::getInstance().getFirstTexti("./flows/symplecticFlow/integrator/turns")*(set->getParticlesNumber())));
		tracking_set->setInitialDistribution();
		this->tracking_store = tracking_set->getStore();
	}
	else if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "false")
	{
		is_phase_portrait= false;
		tracking_set = NULL;
		tracking_store = NULL;
	}
	
	LOG(ITERATOR, Iterator instanciated.);
}

Integrator::~Integrator()
{
	LOG(INTEGRATOR, Destructor called.);
	delete [] A1_;
	delete [] A2_;
	delete [] A3_;
	delete param_evolution;
}

void
Integrator::integrate()
{
    LOG(INTEGRATOR, Integrate...);
		
    // Definition of turn by turn variables
    n_lost = 0;

    // Other parameters
    int particles_number = part_set->getParticlesNumber();
 
    // To output the parameters to a file
    std::ofstream file;
    file.open((this->file_name+".dat").c_str());

		// To output the energy to a file
		std::ofstream file_energy;
		if(is_phase_portrait)
		{
			file_energy.open("energy.dat");
		}
    
    // For performance estimation
    time_t t1, t2;
    static_cast<void>(time(&t1));
    
    // Output headers if needed
    #ifdef DEBUG_FLAG
        std::cout << std::setw(5) << "#" << std::setw(SPACING) << "A1" 
                                         << std::setw(SPACING) << "A2"
                                         << std::setw(SPACING) << "A3"
                                         << std::endl;
    #endif
    
    // Log the output in the file
    #ifndef ITERATOR_TO_FILE_FLAG
        file << "The ITERATOR_TO_FILE_FLAG is not set; will no write out the parameters" << std::endl;
    #else
        file << std::setw(5) << "#" << std::setw(SPACING) << "A1" 
                                    << std::setw(SPACING) << "A2"
                                    << std::setw(SPACING) << "A3"
                                    << std::endl;
    #endif
    
    //
    // Iterations --- Turns
    //
    double A1_t, A2_t, A3_t;
    for(t = 0; t < turns; t++)
    { 
		    // Evaluate the variables
			  A1_t = evaluateVariable(A1_[0]);
				A2_t = evaluateVariable(A2_[0]);
				A3_t = evaluateVariable(A3_[0]);
				std::vector<double> params;
				params.push_back(A1_t);
				params.push_back(A2_t);
				params.push_back(A3_t);
				
        // Check in case of intermediate processing needed
        // Can be optimized with a flag !!!
				#ifdef INTERMEDIATE_TURNS_FLAG
        if(t % intermediate_turns == 0 && t != 0)
        {
        	Log::getInstance().log("INTEGRATOR", "Intermediate processing ?");
          flow->process("./flows/symplecticFlow/dataProcessing/intermediate/");
        }            
				#endif
        // Output if needed
        #ifdef DEBUG_FLAG
        	std::cout << std::setw(5) << std::setprecision(5) << std::scientific << t
                    << std::setw(SPACING) << std::setprecision(5) << std::scientific << A1_t
                    << std::setw(SPACING) << std::setprecision(5) << std::scientific << A2_t
                    << std::setw(SPACING) << std::setprecision(5) << std::scientific << A3_t
                  	<< std::endl;
        #endif
        // Log the output in the file
        #ifdef ITERATOR_TO_FILE_FLAG
        file << std::setw(5) << std::setprecision(5) << std::scientific << t
             << std::setw(SPACING) << std::setprecision(5) << std::scientific << A1_t
             << std::setw(SPACING) << std::setprecision(5) << std::scientific << A2_t
             << std::setw(SPACING) << std::setprecision(5) << std::scientific << A3_t
        	   << std::endl;
        #endif

				for(int m=0; m < particles_number; m++)
        {       
					if(is_phase_portrait)
					{
						(*tracking_store)(t*particles_number+m,ANGLE_INDEX) = (*store)(m,ANGLE_INDEX);
						(*tracking_store)(t*particles_number+m,ACTION_INDEX) = (*store)(m,ACTION_INDEX);
						file_energy << std::setw(SPACING) << std::setprecision(5) << std::scientific 
				          			<< my_integrator->getEnergy((*store)(m,ANGLE_INDEX),(*store)(m,ACTION_INDEX),params);
					}
					if((*store)(m,X_INDEX) != LOST || (*store)(m,XP_INDEX) != LOST)
				  {
						for(int iter=0; iter < (1/timestep); iter++)
						{
							double PHI = (*store)(m,ANGLE_INDEX);
							double J = (*store)(m,ACTION_INDEX);
							my_integrator->integrate(PHI,J,params);
							double theta = PHI/(2*M_PI)+0.5;
							int theta_i = static_cast<int>(theta);
							if(theta < 0)
							{
								theta_i = theta_i - 1;
							}
							PHI=PHI - theta_i * 2 * M_PI;
							(*store)(m,ANGLE_INDEX) = PHI;
							(*store)(m,ACTION_INDEX) = J;
						} // Iterations
					} // If not lost
				} // Particles
				
				if(is_phase_portrait)
				{
					file_energy << std::endl;
		 		}
		} // Turns
		
		// Output
		separatrixSorting(particles_number, A1_t, A2_t, A3_t);
							
    // Performance estimation
    static_cast<void>(time(&t2));
    Log::getInstance().warn("INTEGRATOR", "Duration of the iteration process: ", "s", static_cast<int>(t2 - t1));
    
    // Close the output file
    file.close();
		if(is_phase_portrait)
		{
			file_energy.close();
		}
    LOG(ITERATOR, End of the iteration process);
}

void Integrator::separatrixSorting(int particles_number, double A1_t, double A2_t, double A3_t)
{
	std::ofstream file_upper, file_island, file_lower;
	file_upper.open("upper.dat");
	file_lower.open("lower.dat");
	file_island.open("island.dat");

	for(int m=0; m < particles_number; m++)
	{
		double PHI, J;
		PHI = (*store)(m,ANGLE_INDEX);
		J = (*store)(m,ACTION_INDEX);
		double energy = -(1+A3_t) * cos(PHI)+(0.5 * (J-A2_t)* (J-A2_t)); 
		double J_upper_separatrix = A2_t +2*cos(PHI/2) * sqrt(1+A3_t); 
		double J_lower_separatrix = A2_t -2*cos(PHI/2) * sqrt(1+A3_t);
		if(J > J_upper_separatrix)
		{
			file_upper << PHI << " " << J << " "<< energy <<std::endl;
		}
		else if(J < J_lower_separatrix)
		{
			file_lower << PHI << " " << J << " " << energy <<std::endl;
			
		}
		else
		{
			file_island << PHI << " " << J << " " << energy <<std::endl; 						
		}
	}
	file_upper.close();
	file_lower.close();
	file_island.close();
}

inline double
Integrator::evaluateVariable(CurveParameters& p)
{
    if(p.curve_type == 'l')
    {
      	  return p.start+p.inc*t;
    }
    else if(p.curve_type == 'f')
    {
        return p.value;
    }
    else if(p.curve_type == 'p')
    {
        if(t < p.step_t)
        {
            return p.a0 + p.a1 * t + p.aq * pow(t,p.power1);
        }
        else
        {
            return p.b0 + p.b1 * (t-p.step_t+1) + p.bq * pow(t-p.step_t+1,p.power2);
        }
    }
    else if(p.curve_type == 'c')
    {
        if(p.is_on)
        { 
            p.t_on++;
            return p.amplitude*sin(p.t_on*p.t_on*p.frequencyIncrement);
        }
        else
        {
            return 0.0;
        }
    }
	else if(p.curve_type == 's')
	{
		if(p.is_on)
		{
			p.t_on++;
			return p.amplitude*sin(p.t_on*p.frequency);
		}                                                       
		else
		{
			return 0.0;
		}
	}
	else if(p.curve_type == 'd')
	{
	    if(p.is_on)
	    {
	        p.t_on++;
	        return p.values[p.t_on];
	    }
	    else
	    {
	        return 0.0;
	    }
	}
    else
        return 0.0;
}

void 
Integrator::readParameters()
{
    char curve_info[] = "t";
    #ifdef DEBUG_FLAG
        std::string tmp = "";
    #endif
    LOG(ITERATOR, Reading the input parameters...);        

    // File name
    LOG(ITERATOR, Reading file name...);
    this->file_name = "iterations";//XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/iterator", "file");
        
    // Turns
    LOG(ITERATOR, Reading turns...);
    turns = XMLInputParser::getInstance().getFirstTexti("./flows/symplecticFlow/integrator/turns");
        
    // Intermediate turn    
		#ifdef INTERMEDIATE_TURNS_FLAG
    	LOG(ITERATOR, Reading intermediate turns...);
    	intermediate_turns = XMLInputParser::getInstance().getFirstTexti("./flows/symplecticFlow/integrator/intermediateTurns");
    #endif

    // Timestep
    LOG(ITERATOR, Reading timestep...);
    timestep = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/integrator/timestep");
 
 		// Create the containers for the other parameters
    createContainers();
            
    LOG(ITERATOR, Reading evolv for A1);
		curve_info[0] = param_evolution->readEvolvingVariable(static_cast<std::string>("A1"),A1_[0]);
		#ifdef DEBUG_FLAG
		    tmp = "A1 curve type is: ";
		    tmp += curve_info;
		    LOGV(ITERATOR, tmp);
		#endif
		
    LOG(ITERATOR, Reading evolv for A2);
		curve_info[0] = param_evolution->readEvolvingVariable(0,static_cast<std::string>("A2"),A2_);
		#ifdef DEBUG_FLAG
		    tmp = "A2 curve type is: ";
		    tmp += curve_info;
		    LOGV(ITERATOR, tmp);
		#endif
		
		LOG(ITERATOR, Reading evolv for A3);
		curve_info[0] = param_evolution->readEvolvingVariable(0,static_cast<std::string>("A3"),A3_);
		#ifdef DEBUG_FLAG
		    tmp = "A3 curve type is: ";
		    tmp += curve_info;
		    LOGV(ITERATOR, tmp);
		#endif

    LOG(ITERATOR, Input parameters read.);
}

void
Integrator::createContainers()
{
    LOG(ITERATOR, Creating containers...);
    A1_ = new CurveParameters[n_kicks_];
    A2_ = new CurveParameters[n_kicks_];
    A3_ = new CurveParameters[n_kicks_];
    LOG(ITERATOR, Creating containers: done.);
}

int
Integrator::getTurns()
{
    return turns;
}          

TrackingDistribution* Integrator::getTrackingSet()
{
	return tracking_set;
}
