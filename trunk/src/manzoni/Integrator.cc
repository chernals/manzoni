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

#include <Integrator.h>
#include <flows/SymplecticFlow.h>

bool sorting_comparison (std::vector<double> i,std::vector<double> j) { return (i.at(0)<j.at(0)); }

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
	t = 0.0;
	dt = 0.0;
	 
	// Type of integrator
	std::string integrator_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator/hamiltonian0", "type");
	if(integrator_type == "pendulum")
	{
		my_integrator = new Pendulum(timestep, XMLInputParser::getInstance().getAttribute("./flows/symplecticFlow/integrator", "order"));
	}
	else if(integrator_type == "interpolatingHenon4")
	{
		my_integrator = new InterpolatingHenon4(timestep, XMLInputParser::getInstance().getAttribute("./flows/symplecticFlow/integrator", "order"));
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
		tracking_set->setInitialDistribution(); // Empty distribution
		this->tracking_store = tracking_set->getStore();
	}
	else if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow", "phasePortrait") == "false")
	{
		is_phase_portrait= false;
		tracking_set = NULL;
		tracking_store = NULL;
	}
	
	// Compute the adiabatic invariant ?
	if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator/hamiltonian0", "invariant") == "true")
	{
		adiabatic_invariant = true;
		if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator/hamiltonian0", "invariantallturns") == "true")
		{
			invariant_all_turns = true;
		}
		else
		{
			invariant_all_turns = false;
		}
		if(XMLInputParser::getInstance().isFoundFromPath("./flows/symplecticFlow/integrator/hamiltonian0/invariantPasses"))
		{
			invariant_passes = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/integrator/hamiltonian0/invariantPasses");
		}
		else
		{
			invariant_passes = 10;
		}
	}
	else
	{
		adiabatic_invariant = false;
	}
	
	#ifndef INVARIANT_FLAG
	if(adiabatic_invariant == true)
	{
		ERROR(INTEGRATOR, Manzoni was not compiled with support for the invariant computation !);
		throw("No support for the invariant computation !");
	}
	#endif
	
	LOG(INTEGRATOR, Integrator instanciated.);
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
		std::vector<std::ofstream*> file_energy;

		if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dumpenergy") == "true")
		{
			file_energy.resize(particles_number);
			for(int i = 0; i < particles_number; i++)
			{
				file_energy.at(i) = new std::ofstream((std::string("energy")+i2s(i)+std::string(".dat")).c_str());
			}
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
    
		// Open some output files
		#ifdef INVARIANT_FLAG
		std::ofstream file_invariant;
		if(adiabatic_invariant)
		{
			file_invariant.open("invariant.dat");
		}
		#endif

		std::ofstream file_dump_all;		
		if(is_phase_portrait)
		{
			file_dump_all.open("dumpall.dat");
		}
		

    //
    // WITH INVARIANT COMPUTATION
    //
    //
		if(adiabatic_invariant)
		{
			// Iterations --- Turns
	    //   
			double A1_t, A2_t, A3_t;
	    for(t = 0.0; t < turns; t++)
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
					#ifdef INTERMEDIATE_TURNS_FLAG
	        if(t % intermediate_turns == 0 && t != 0)
	        {
	        	Log::getInstance().log("INTEGRATOR", "Intermediate processing ?");
	          flow->process("./flows/symplecticFlow/dataProcessing/intermediate/");
	
						if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "separatrix") == "true")
							separatrixSorting(particles_number, A1_t, A2_t, A3_t,my_integrator);
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

					#ifdef INVARIANT_FLAG
					std::vector<double> invariants;
					invariants.resize(particles_number);
					bool invariant_test = true;
					if(invariant_all_turns == false)
					{
						if (t == 0 || t == (turns-1))
							invariant_test = true;
						else
							invariant_test = false;
					}
					#endif

					// Particles main loop
					for(int m=0; m < particles_number; m++)
	        {       
						double PHI = (*store)(m,ANGLE_INDEX);
						double J = (*store)(m,ACTION_INDEX);
						if(is_phase_portrait)
						{
							(*tracking_store)(t*particles_number+m,ANGLE_INDEX) = PHI;
							(*tracking_store)(t*particles_number+m,ACTION_INDEX) = J;
						}

						// Timesteps iterations
						for(int iter=0; iter < (1.0/timestep); iter++)
						{
							// Used by evaluateVariable
							dt = iter * timestep;

							// If phase_portrait, then output more things
							if(is_phase_portrait)
							{
								if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dump") == "true")
								{
									file_dump_all << std::setprecision(15) << std::scientific << PHI << " " << J << std::endl;
								}		
							}
					
							// Do the integration itself
							my_integrator->integrate(PHI,J,params,t,iter,turns);

						} // Iterations

						if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dumpenergy") == "true")
						{
							(*(file_energy.at(m))) << t << " " << std::setprecision(15) << std::scientific 
							<< my_integrator->getEnergy(PHI,J,params, t, (1.0/timestep)-1,turns) << std::endl;
						}

						// Rescaling of PHI
						double theta = PHI/(2*M_PI)+0.5;
						int theta_i = static_cast<int>(theta);
						if(theta < 0)
							theta_i = theta_i - 1;
						PHI=PHI - theta_i * 2 * M_PI;
						//

						// Store the new coordinates
						(*store)(m,ANGLE_INDEX) = PHI;
						(*store)(m,ACTION_INDEX) = J;

						// Compute the adiabatic invariant
						#ifdef INVARIANT_FLAG
						if(adiabatic_invariant && invariant_test)
						{
							bool condition;
							double PHIZERO = PHI;
							my_integrator->integrate(PHI,J,params,t,0,turns);
							std::vector<std::vector<double> > orbit_left;
							std::vector<std::vector<double> > orbit_right;
							double DELTA_PHI = PHI - PHIZERO;
							for(int loop=0;loop < invariant_passes; loop++)
							{
								int inner_iterations = 0;
								do
								{
									inner_iterations++;
									double PHI_BEFORE = PHI;								
									my_integrator->integrate(PHI,J,params,t,0,turns);
									double theta = PHI/(2*M_PI)+0.5;
									int theta_i = static_cast<int>(theta);
									if(theta < 0)
										theta_i = theta_i - 1;
									double DELTA = PHI-PHI_BEFORE;
									PHI=PHI - theta_i * 2 * M_PI;
									std::vector<double> tmp;
									tmp.push_back(PHI);
									tmp.push_back(J);
									if(DELTA > 0)
									{
										orbit_left.push_back(tmp);
									}
									else
									{
										orbit_right.push_back(tmp);
									}
									if(theta_i != 0)
										PHI_BEFORE = PHI; // 'Deactivate' the condition
									// This is the STOP condition (STOP IF TRUE)
									condition = ((PHI <= PHIZERO && PHI_BEFORE >= PHIZERO && ((PHI-PHI_BEFORE)*DELTA_PHI) > 0) || 
									             (PHI >= PHIZERO && PHI_BEFORE <= PHIZERO && ((PHI-PHI_BEFORE)*DELTA_PHI) > 0));
									if(inner_iterations > 1000000)
										break;
								}
								while(!condition);
							}
							std::sort(orbit_left.begin(), orbit_left.end(), sorting_comparison);
							std::sort(orbit_right.begin(), orbit_right.end(), sorting_comparison);
							double invariant_left = 0;
							double invariant_right = 0;
							for(int i = 1; i < static_cast<int>(orbit_left.size()); i++)
							{
								invariant_left += (orbit_left.at(i).at(0)-orbit_left.at(i-1).at(0))*0.5*(orbit_left.at(i-1).at(1)+orbit_left.at(i).at(1));
							}
							for(int i = 1; i < static_cast<int>(orbit_right.size()); i++)
							{
								invariant_right += (orbit_right.at(i).at(0)-orbit_right.at(i-1).at(0))*0.5*(orbit_right.at(i-1).at(1)+orbit_right.at(i).at(1));
							}
							if(invariant_left > invariant_right)
							{
								invariants.at(m) = invariant_left - invariant_right;
							}
							else
							{
								invariants.at(m) = invariant_right - invariant_left;
							}
						} // Invariant
						#endif

					} // Particles

					// Case where the computation of the adiabatic invariant is needed
					#ifdef INVARIANT_FLAG
					if(adiabatic_invariant && invariant_test)
					{
						file_invariant << t << " ";
						for(int m=0; m < particles_number; m++)
	        	{
							file_invariant << std::setprecision(25) << invariants.at(m) << " ";
						}
						file_invariant << std::endl;
					}
					#endif

			} // Turns

			// Output
			if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "separatrix") == "true")
				separatrixSorting(particles_number, A1_t, A2_t, A3_t,my_integrator);

	    // Performance estimation
	    static_cast<void>(time(&t2));
	    Log::getInstance().warn("INTEGRATOR", "Duration of the iteration process: ", "s", static_cast<int>(t2 - t1));

	    // Close the output file
	    file.close();
			#ifdef INVARIANT_FLAG
			file_invariant.close();
			#endif
		}

    //
    // WITHOUT INVARIANT COMPUTATION
    //
    //
		else
		{
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
			#ifdef INTERMEDIATE_TURNS_FLAG
      if(t % intermediate_turns == 0 && t != 0)
      {
      	Log::getInstance().log("INTEGRATOR", "Intermediate processing ?");
        flow->process("./flows/symplecticFlow/dataProcessing/intermediate/");

				if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "separatrix") == "true")
					separatrixSorting(particles_number, A1_t, A2_t, A3_t,my_integrator);
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
			
			// Particles main loop
			for(int m=0; m < particles_number; m++)
      {       
				double PHI = (*store)(m,ANGLE_INDEX);
				double J = (*store)(m,ACTION_INDEX);
				if(is_phase_portrait)
				{
					(*tracking_store)(t*particles_number+m,ANGLE_INDEX) = PHI;
					(*tracking_store)(t*particles_number+m,ACTION_INDEX) = J;
				}
				
				// Timesteps iterations
				for(int iter=0; iter < (1.0/timestep); iter++)
				{
					// Used by evaluateVariable
					dt = iter * timestep;
					
					// If phase_portrait, then output more things
					if(is_phase_portrait)
					{
						if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dump") == "true")
						{
							file_dump_all << std::setprecision(15) << std::scientific << PHI << " " << J << std::endl;
						}
					}

					
					// Do the integration itself
					my_integrator->integrate(PHI,J,params,t,iter,turns);
					
				} // Iterations
				
				if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dumpenergy") == "true")
				{
						(*(file_energy.at(m))) << t << " " << std::setprecision(15) << std::scientific 
							<< my_integrator->getEnergy(PHI,J,params, t, (1.0/timestep)-1,turns) << std::endl;
				}
				
				// Rescaling of PHI
				double theta = PHI/(2*M_PI)+0.5;
				int theta_i = static_cast<int>(theta);
				if(theta < 0)
					theta_i = theta_i - 1;
				PHI=PHI - theta_i * 2 * M_PI;
				//
				
				// Store the new coordinates
				(*store)(m,ANGLE_INDEX) = PHI;
				(*store)(m,ACTION_INDEX) = J;
				
			} // Particles
						
			
	} // Turns
	// Output
	if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "separatrix") == "true")
		separatrixSorting(particles_number, A1_t, A2_t, A3_t,my_integrator);
						
  // Performance estimation
  static_cast<void>(time(&t2));
  Log::getInstance().warn("INTEGRATOR", "Duration of the iteration process: ", "s", static_cast<int>(t2 - t1));
  
  // Close the output file
  file.close();
	
	
	} // Big if of the invariant

		if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "dumpenergy") == "true")
		{
			for(int i = 0; i < particles_number; i++)
			{
				file_energy.at(i)->close();
			}	
		}

		file_dump_all.close();
		
    LOG(INTEGRATOR, End of the iteration process);
}

void Integrator::adiabaticInvariant()
{
	
}

void Integrator::separatrixSorting(int particles_number, double A1_t, double A2_t, double A3_t,BaseIntegrator* integrator)
{
	static int index = 0;
	std::ofstream file_upper, file_island, file_lower;
	char buff[30];
	sprintf(buff, "%04d", index);
	file_upper.open((std::string("upper")+std::string(buff)+std::string(".dat")).c_str());
	file_lower.open((std::string("lower")+std::string(buff)+std::string(".dat")).c_str());
	file_island.open((std::string("island")+std::string(buff)+std::string(".dat")).c_str());
	index++;

	std::vector<double> params;	
	params.push_back(A1_t);
	params.push_back(A2_t);
	params.push_back(A3_t);
	
	std::string integrator_type = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator/hamiltonian0", "type");

	if(integrator_type == "pendulum")
	{
		for(int m=0; m < particles_number; m++)
		{
			double PHI, J;
			PHI = (*store)(m,ANGLE_INDEX);
			J	 = (*store)(m,ACTION_INDEX);
			double energy = my_integrator->getEnergy(PHI,J,params, 0, 0, 0);
				
			double J_upper_separatrix = A2_t +2*cos(PHI/2) * sqrt(1+A3_t); 
			double J_lower_separatrix = A2_t -2*cos(PHI/2) * sqrt(1+A3_t);
		
		if(J > J_upper_separatrix)
		{
			file_upper << PHI << " " << J << " " << energy << std::endl; 
		}
		else if(J < J_lower_separatrix)
		{
			file_lower << PHI << " " << J << " " << energy << std::endl;
			
		}
		else
		{
			file_island << PHI << " " << J << " " << energy << std::endl;					
		}
	} // For
	
} // If pendulum
else if(integrator_type == "interpolatingHenon4")
{
		for(int m=0; m < particles_number; m++)
		{
			double PHI, J;
			PHI = (*store)(m,ANGLE_INDEX);
			J	 =  (*store)(m,ACTION_INDEX);
			
		// This is for the interpolating Hamiltonian
		// At 0.245
		// For kappa=-1.1
		double J_upper_separatrix = (0.0314159 + sqrt(0.0000549474 + 0.0000549474 * cos(4*PHI))) / (0.275886 - 0.016265 * cos(4*PHI));
		double J_lower_separatrix = (0.0314159 - sqrt(0.0000549474 + 0.0000549474 * cos(4*PHI))) / (0.275886 - 0.016265 * cos(4*PHI));
		
		if(J > J_upper_separatrix)
		{
			file_upper << PHI << " " << J << std::endl; 
		}
		else if(J < J_lower_separatrix)
		{
			file_lower << PHI << " " << J << std::endl;		
		}
		else
		{
			file_island << PHI << " " << J << std::endl; 				
		}
	} // For
	
	
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
      	  return p.start+p.inc*(t+dt);
    }
    else if(p.curve_type == 'f')
    {
        return p.value;
    }
    else if(p.curve_type == 'p')
    {
        if(t < p.step_t)
        {
            return p.a0 + p.a1 * (t+dt) + p.aq * pow((t+dt),p.power1);
        }
        else
        {
            return p.b0 + p.b1 * ((t+dt)-p.step_t+1) + p.bq * pow((t+dt)-p.step_t+1,p.power2);
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
		if(p.is_on && p.t_on >= p.turnsLength)
		{
			p.is_on = false;
		}
	    if(p.is_on)
	    {
		      p.t_on++;
	        return p.values[p.t_on-1];
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
    LOG(INTEGRATOR, Reading the input parameters...);        

    // File name
    LOG(INTEGRATOR, Reading file name...);
    this->file_name = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/symplecticFlow/integrator", "file");
        
    // Turns
    LOG(INTEGRATOR, Reading turns...);
    turns = XMLInputParser::getInstance().getFirstTexti("./flows/symplecticFlow/integrator/turns");
        
    // Intermediate turn    
		#ifdef INTERMEDIATE_TURNS_FLAG
    	LOG(INTEGRATOR, Reading intermediate turns...);
    	intermediate_turns = XMLInputParser::getInstance().getFirstTexti("./flows/symplecticFlow/integrator/intermediateTurns");
    #endif

    // Timestep
    LOG(INTEGRATOR, Reading timestep...);
    timestep = XMLInputParser::getInstance().getFirstTextd("./flows/symplecticFlow/integrator/timestep");
 
 		// Create the containers for the other parameters
    createContainers();
            
    LOG(INTEGRATOR, Reading evolv for A1);
		curve_info[0] = param_evolution->readEvolvingVariable(static_cast<std::string>("A1"),A1_[0]);
		A1_[0].is_on = true;
		#ifdef DEBUG_FLAG
		    tmp = "A1 curve type is: ";
		    tmp += curve_info;
		    LOGV(INTEGRATOR, tmp);
		#endif
		
    LOG(INTEGRATOR, Reading evolv for A2);
		curve_info[0] = param_evolution->readEvolvingVariable(0,static_cast<std::string>("A2"),A2_);
		A2_[0].is_on = true;
		#ifdef DEBUG_FLAG
		    tmp = "A2 curve type is: ";
		    tmp += curve_info;
		    LOGV(INTEGRATOR, tmp);
		#endif
		
		LOG(INTEGRATOR, Reading evolv for A3);
		curve_info[0] = param_evolution->readEvolvingVariable(0,static_cast<std::string>("A3"),A3_);
		A3_[0].is_on = true;
		#ifdef DEBUG_FLAG
		    tmp = "A3 curve type is: ";
		    tmp += curve_info;
		    LOGV(INTEGRATOR, tmp);
		#endif

    LOG(INTEGRATOR, Input parameters read.);
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
