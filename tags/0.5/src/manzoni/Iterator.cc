#include <Iterator.h>
#include <AlessandroFlow.h>

Iterator::Iterator(ParticlesSet* set, Flow* f)
{
	LOG(ITERATOR, Instanciating iterator);
	
	// Particle store
	this->flow = f;
	this->part_set = set;
	this->store = part_set->getStore();
	
	// Tracking
	tracking_ic = this->part_set->getTracking();
	is_tracked = this->part_set->isTracking();
	n_tracked = this->part_set->getTrackingNumber();
	
	// Coherent motion
	coherent_parameters = this->part_set->getCoherent();
	is_coherent = this->part_set->isCoherent();
	
	// Read the input parameters in the XML file
	param_evolution = new ParametersEvolution();
	readParameters();
	
	// Set state variables
	t = 0;
	
	// Prepare the graph for the evolution of the parameters
	LOG(ITERATOR, Preparing the plot of the parameters);
	if(XMLInputParser::getInstance().getAttributeTextFromPath("./flows/alessandroFlow/iterator", "plot") == "yes")
	{
		with_plots_ = true;
		graph_ = new ProcGraph(turns, "parameters", "Evolution of the parameters");
	}
	else
	{
		with_plots_ = false;
		graph_ = NULL;
	}
	
	LOG(ITERATOR, Iterator instanciated.);
}

Iterator::~Iterator()
{
	LOG(ITERATOR, Destructor called.);
	if(part_set->getDims() == TWO_DIMENSIONS || part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
	{
		for(int u = 0; u < n_kicks_; u++)
		{
			delete [] table_tune_h[u];
			delete [] table_kappa[u];
		}
		delete [] table_tune_h;
		delete [] table_kappa;
		delete [] hor_tune_;
		delete [] kappa_;
		#ifdef DAMPER_FLAG
		if(isDamperSet)
		{
		    delete [] damper;
		    for(int u = 0; u < n_kicks_; u++)
		    {
		        delete [] table_damper[u];
		    }
		    delete [] table_damper;
		}
		#endif
	}
	if(part_set->getDims() == 4 || part_set->getDims() == 6)
	{
		for(int u = 0; u < n_kicks_; u++)
		{
			delete [] table_tune_v[u];
		}
		delete [] table_tune_v;
		delete [] ver_tune_;
	}
	if(part_set->getDims() == 6)
	{
		delete [] table_tune_l[0];
		delete [] table_tune_l;
		delete [] hor_chromaticity;
	    delete [] ver_chromaticity;
	    delete [] dispersion;
	    delete [] dispersion_prime;
	}  
	if(graph_ != NULL)
	{
	    delete graph_;
	}
	delete param_evolution;
}

void
Iterator::iterate2d()
{
    LOG(ITERATOR, Computation of the iterations using the 2d iterator...);
    
    // Definition of turn by turn variables
    double centroid_x = 0.0;
    double centroid_xp = 0.0;
    n_lost = 0;
    
    // Definition of the kick by kick variables
    double x_cos, x_sin; 
    double omega_x;
    double iota;
    double omicron;
    int rotation_kick_id;
    #ifdef DAMPER_FLAG
        double damper_kick = 0.0;
    #endif
    
    // Definition of the particle by particle variables
    double x, xp;
    double kick_x;
    
    // Physical constants
    double ppi = 2 * M_PI;
    
    // Other parameters
    int particles_number = part_set->getParticlesNumber();
    
    // To output the parameters to a file
    std::ofstream file;
    file.open((this->file_name+".dat").c_str());
    
    // For performance estimation
    time_t t1, t2;
    static_cast<void>(time(&t1));
	// time(&t1); // Test that !
    
    // Adjust the position of the islands
    LOG(ITERATOR, Rotation of the coordinates (R));
    phaseRotation(horizontal_phase_modification,X_PLANE);

    // Output headers if needed
    #ifdef DEBUG_FLAG
        std::cout << std::setw(5) << "#" << std::setw(5) << "Kick" << std::setw(5) << "Type" << std::setw(SPACING) << "Tune_h" << std::setw(13) << "Kappa" << std::setw(13) << "Norm" << std::setw(13) << "Iota" << std::setw(13) << "Omicron" 
        #ifdef DAMPER_FLAG
        << std::setw(SPACING) << "Damper"
        #endif
        << std::endl;
    #endif
    
    // Log the output in the file
    #ifndef ITERATOR_TO_FILE_FLAG
        file << "The ITERATOR_TO_FILE_FLAG is not set; will no write out the parameters" << std::endl;
    #else
        file << std::setw(5) << "#" << std::setw(5) << "Kick" << std::setw(5) << "Type" << std::setw(13) << "Tune_h" << std::setw(13) << "Kappa" << std::setw(13) << "Norm" << std::setw(13) << "Iota" << std::setw(13) << "Omicron" 
        #ifdef DAMPER_FLAG
        << std::setw(13) << "Damper"
        #endif
        << std::endl;
    #endif
    
    //
    // Iterations --- Turns
    //
    for(t = 0; t < turns; t++)
    { 
        // Check in case of intermediate processing needed
        // Can be optimized with a flag !!!
		#ifdef INTERMEDIATE_TURNS_FLAG
        	if(t % intermediate_turns == 0 && t != 0)
        	{
            	Log::getInstance().log("ITERATOR", "Intermediate processing ?");
            	flow->process("./flows/alessandroFlow/dataProcessing/intermediate/");
        	}            
		#endif
        
        //
        // Iterations --- Kicks
        //
        for(int kick_id=0; kick_id < n_kicks_; kick_id++)
        {
            //
            // Damper
            //
            #ifdef DAMPER_FLAG
            if(isDamperSet[kick_id])
            {
                // This uggly condition test if the damper is on for the given turn
                if((t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) <= damper[kick_id].turnsLength && (t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) >= 0 && damper[kick_id].index < damper[kick_id].repetitions && t >= damper[kick_id].turnStart)
                {
                    damper[kick_id].is_on = true;
                }
                else // Here the damper has to be off
                {
                    // Test if we are jumping to the next damper sequence
                    if((t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) == damper[kick_id].turnsLength+1)
                    {
                        damper[kick_id].index++; // Next "action" of the damper
                        damper[kick_id].t_on = 0;
                    }
                    damper[kick_id].is_on = false;
                }
                // In any case evaluate the variable (will be zero anyway if the damper is off)
                damper_kick = evaluateVariable(damper[kick_id]);
                table_damper[kick_id][t] = damper_kick;
            }
            #endif
        
            // Compute the turn by turn variables
            // Tune horizontal
            rotation_kick_id = (kick_id+1) % n_kicks_; // Will apply the rotation transporting to the next kick
            table_tune_h[rotation_kick_id][t] = evaluateVariable(hor_tune_[rotation_kick_id]);
            omega_x = ppi * table_tune_h[rotation_kick_id][t];
            x_cos = cos(omega_x);
            x_sin = sin(omega_x);
            
            // Iota and omicron
            if(kick_types[kick_id] == "SO")
            {
                table_kappa[kick_id][t] = evaluateVariable(kappa_[kick_id]);
                iota = normalizations[kick_id]/normalizations[0];
                omicron = iota * iota * table_kappa[kick_id][t];
            }
            else if(kick_types[kick_id] == "O")
            {
                table_kappa[kick_id][t] = 0;
                iota = normalizations[kick_id]/normalizations[0]; // Hu ?!
                omicron = iota * iota; // Got confused ?
                iota = 0; // Dirty to use iota as a tmp variable isn't ?
            }
            else if(kick_types[kick_id] == "L")
            {
                table_kappa[kick_id][t] = 0;
                iota = 0;
                omicron = 0;
            }
            else // This is a "just in case"
            {
                table_kappa[kick_id][t] = 0;
                omicron = 0;
                iota = 0;
            }
                        
            // Output if needed
            #ifdef DEBUG_FLAG
            std::cout << std::setw(5) << std::setprecision(5) << std::scientific << t
                      << std::setw(5) << std::setprecision(5) << std::scientific << kick_id
                      << std::setw(5) << kick_types[kick_id]
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << table_tune_h[rotation_kick_id][t]
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << table_kappa[kick_id][t]
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << normalizations[kick_id]
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << iota
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << omicron
                      #ifdef DAMPER_FLAG
                      << std::setw(SPACING) << std::setprecision(5) << std::scientific << damper_kick
                      #endif
                      << std::endl;
            #endif
            // Log the output in the file
            #ifdef ITERATOR_TO_FILE_FLAG
            file << std::setw(5) << std::setprecision(5) << std::scientific << t
                 << std::setw(5) << std::setprecision(5) << std::scientific << kick_id
                 << std::setw(5) << kick_types[kick_id]
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << table_tune_h[rotation_kick_id][t]
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << table_kappa[kick_id][t]
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << normalizations[kick_id]
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << iota
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << omicron
                 #ifdef DAMPER_FLAG
                 << std::setw(SPACING) << std::setprecision(5) << std::scientific << damper_kick
                 #endif
                 << std::endl;
            #endif
            
            //
            // Iterations --- Particles
            //
            
            // If damper_kick non zero       
            #ifdef DAMPER_FLAG            
				// Out of the main loop to avoid testing the condition if not needed
            	if(damper_kick != 0)
            	{
                	for(int m=0; m < particles_number; m++)
                	{
                    	(*store)(m,XP_INDEX) += damper_kick;
                	}
            	}
            #endif
          
            for(int m=0; m < particles_number; m++)
            {
                // Faster to check that than to compute invalid particles...
                if((*store)(m,X_INDEX) != LOST || (*store)(m,XP_INDEX) != LOST)
                {     
                    /*
                        !!! WARNING !!!
                        This is were 'real' physics happens !
                    */      

                    // Non linear kick
                    kick_x =   (*store)(m,XP_INDEX) 
                            + iota * (*store)(m,X_INDEX)*(*store)(m,X_INDEX)
                            + omicron * (*store)(m,X_INDEX) * (*store)(m,X_INDEX) * (*store)(m,X_INDEX);   
                    
                    // Linear rotation
                    x =  x_cos*(*store)(m,X_INDEX) + x_sin*kick_x;
                    xp= -x_sin*(*store)(m,X_INDEX) + x_cos*kick_x;
            
                    // Check and assign new coordinates
                    if(x >= 1 || xp >= 1)
                    {
                        (*store)(m,X_INDEX) = LOST;
                        (*store)(m,XP_INDEX) = LOST;
                        n_lost++;
                    }
                    else
                    {
                        (*store)(m,X_INDEX) = x;
                        (*store)(m,XP_INDEX) = xp;
                    }
                }
            } // Iterations particles
        } // Iterations kicks

        //
        // Tracking and coherent motion
        //
        if(is_tracked)
        {
            for(int i = 0; i < n_tracked; i++)
            {
                tracking_ic[i][t][X_INDEX]  = (*store)(particles_number-n_tracked+i,X_INDEX);
                tracking_ic[i][t][XP_INDEX] = (*store)(particles_number-n_tracked+i,XP_INDEX);
            }
            if(is_coherent)
            {
                centroid_x = (*store)(0,X_INDEX);
                centroid_xp = (*store)(0,XP_INDEX);
                for(int m=1; m < particles_number; m++)
                {
                    centroid_x += (*store)(m,X_INDEX);
                    centroid_xp += (*store)(m,XP_INDEX);
                }
                coherent_parameters[t][X_INDEX] = (centroid_x  - n_lost * LOST) / particles_number;
                coherent_parameters[t][XP_INDEX] = (centroid_xp - n_lost * LOST) / particles_number;
            }
        }
    } // Iterations turns

    // Adjust the position of the islands
    LOG(ITERATOR, Rotation of the coordinates (R^-1));
    phaseRotation(-horizontal_phase_modification,X_INDEX); // Note the minus sign
                 
    // Performance estimation
    static_cast<void>(time(&t2));
    Log::getInstance().warn("ITERATOR", "Duration of the iteration process: ", "s", static_cast<int>(t2 - t1));
    
    // Finish by adding the evolution of the parameters to the graph
    if(with_plots_)
    {
        for(int k = 0; k < n_kicks_; k++)
        {         
	        graph_->draw(table_tune_h[k], "Horizontal tune");
            #ifdef DAMPER_FLAG
                if(isDamperSet[k])
                {
                    graph_->draw(table_damper[k], "Damper");
                }
            #endif
            if(k==n_kicks_-1)
            {
                graph_->draw(table_kappa[k], "Kappa",true);
            }
            else
            {
                graph_->draw(table_kappa[k], "Kappa",false);
            }
        }
    }
    // Close the output file
    file.close();
    LOG(ITERATOR, End of the iteration process);
}

void
Iterator::iterate4d()
{
    LOG(ITERATOR, Computation of the iterations using the 4d iterator...);
    
    // Definition of the turn by turn variables
    double centroid_x  = 0.0;
    double centroid_xp = 0.0;
    double centroid_y  = 0.0;
    double centroid_yp = 0.0;
    n_lost = 0;
    
    // Definition of the kick by kick variables
    double omega_x, omega_y;
    double x_cos, x_sin, y_cos, y_sin;
    double iota;
    double omicron;
    int rotation_kick_id;
    double chik, chik2;
    #ifdef DAMPER_FLAG
    	double damper_kick = 0.0;
    #endif
    
    // Definition of the particle by particle variables
    double x,xp,y,yp;
    double kick_x, kick_y;
    double store_x2_tmp;
    double store_y2_tmp;
    
    // Physical constants
    double ppi = 2 * M_PI;
    
    // Other parameters
    int particles_number = part_set->getParticlesNumber();
    
    // To output the parameters to a file
    #ifdef ITERATOR_TO_FILE_FLAG
    std::ofstream file;
    file.open((this->file_name+".dat").c_str());
    #endif

    // For performance estimation
    time_t t1, t2;
    static_cast<void>(time(&t1));
    
    // Adjust the position of the islands
    LOG(ITERATOR, Rotation of the coordinates (R));
    phaseRotation(horizontal_phase_modification,X_PLANE);
    phaseRotation(vertical_phase_modification,Y_PLANE);
    
    // Output header if needed
    #ifdef DEBUG_FLAG
        std::cout << std::setw(13) << "#" << std::setw(13) << "Kick" << std::setw(13) << "Tune_h" << std::setw(13) << "Tune_v" << std::setw(13) << "Kappa" << std::setw(13) << "Norm" << std::setw(13) << "Iota" << std::setw(13) << "Omicron" 
        #ifdef DAMPER_FLAG
        << std::setw(13) << "Damper"
        #endif
        << std::endl;
    #endif
    
    // Log the output in the file
    #ifdef ITERATOR_TO_FILE_FLAG
        file << std::setw(13) << "#" << std::setw(13) << "Kick" << std::setw(13) << "Tune_h" << std::setw(13) << "Tune_v" << std::setw(13) << "Kappa" << std::setw(13) << "Norm" << std::setw(13) << "Iota" << std::setw(13) << "Omicron" 
        #ifdef DAMPER_FLAG
        << std::setw(13) << "Damper"
        #endif
        << std::endl;
    #endif

    //
    // Iterations --- Turns
    //
    for(t = 0; t < turns; t++)
    {
        //
        // Iterations --- Kicks
        //
        for(int kick_id=0; kick_id < n_kicks_; kick_id++)
        {
            //
            // Damper
            //
            #ifdef DAMPER_FLAG
            if(isDamperSet)
            {
                // This uggly condition test if the damper is on for the given turn
                if((t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) <= damper[kick_id].turnsLength && (t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) >= 0 && damper[kick_id].index < damper[kick_id].repetitions && t >= damper[kick_id].turnStart)
                {
                    damper[kick_id].is_on = true;
                }
                else // Here the damper has to be off
                {
                    // Test if we are jumping to the next damper sequence
                    if((t-damper[kick_id].index*(damper[kick_id].turnsLength+damper[kick_id].turnsInterval)-damper[kick_id].turnStart) == damper[kick_id].turnsLength+1)
                    {
                        damper[kick_id].index++;
                        damper[kick_id].t_on = 0;
                    }
                    damper[kick_id].is_on = false;
                }
                // In any case evaluate the variable (will be zero anyway if the damper is off)
                damper_kick = evaluateVariable(damper[kick_id]);
                table_damper[kick_id][t] = damper_kick;
            }
            #endif
        
            // Compute the turn by turn variables
            rotation_kick_id = (kick_id + 1) %n_kicks_; // Will apply the rotation corresponding to the next kick
            
            // Tune horizontal
            table_tune_h[rotation_kick_id][t] = evaluateVariable(hor_tune_[rotation_kick_id]);
            omega_x = ppi * table_tune_h[rotation_kick_id][t];
            x_cos = cos(omega_x);
            x_sin = sin(omega_x);
            
            // Tune vertical
            table_tune_v[rotation_kick_id][t] = evaluateVariable(ver_tune_[rotation_kick_id]);
            omega_y = ppi * table_tune_v[rotation_kick_id][t];
            y_cos = cos(omega_y);
            y_sin = sin(omega_y);
            
            // Chi
            chik = chi[kick_id];
            chik2 = chik * chik;
            
            // Iota and omicron
            if(kick_types[kick_id] == "SO")
            {
                table_kappa[kick_id][t] = evaluateVariable(kappa_[kick_id]);
                iota = normalizations[kick_id]/normalizations[0];
                omicron = iota * iota * table_kappa[kick_id][t];
            }
            else if(kick_types[kick_id] == "O")
            {
                table_kappa[kick_id][t] = 0;
                iota = normalizations[kick_id]/normalizations[0];
                omicron = iota * iota;
                iota = 0;
            }
            else if(kick_types[kick_id] == "L")
            {
                table_kappa[kick_id][t] = 0;
                iota = 0;
                omicron = 0;
            }
            else
            {
                iota = 0;
                omicron = 0;
            }
           
            // Output if needed
            #ifdef DEBUG_FLAG
            std::cout << std::setw(13) << std::setprecision(5) << std::scientific << t
                      << std::setw(13) << std::setprecision(5) << std::scientific << kick_id
                      << std::setw(13) << std::setprecision(5) << std::scientific << table_tune_h[kick_id][t]
                      << std::setw(13) << std::setprecision(5) << std::scientific << table_tune_v[kick_id][t]
                      << std::setw(13) << std::setprecision(5) << std::scientific << table_kappa[kick_id][t]
                      << std::setw(13) << std::setprecision(5) << std::scientific << normalizations[kick_id]
                      << std::setw(13) << std::setprecision(5) << std::scientific << iota
                      << std::setw(13) << std::setprecision(5) << std::scientific << omicron
                      #ifdef DAMPER_FLAG
                      << std::setw(13) << std::setprecision(5) << std::scientific << damper_kick
                      #endif
                      << std::endl;
            #endif
            // Log the output in the file
            #ifdef ITERATOR_TO_FILE_FLAG
            file << std::setw(13) << std::setprecision(5) << std::scientific << t
                 << std::setw(13) << std::setprecision(5) << std::scientific << kick_id
                 << std::setw(13) << std::setprecision(5) << std::scientific << table_tune_h[kick_id][t]
                 << std::setw(13) << std::setprecision(5) << std::scientific << table_tune_v[kick_id][t]
                 << std::setw(13) << std::setprecision(5) << std::scientific << table_kappa[kick_id][t]
                 << std::setw(13) << std::setprecision(5) << std::scientific << normalizations[kick_id]
                 << std::setw(13) << std::setprecision(5) << std::scientific << iota
                 << std::setw(13) << std::setprecision(5) << std::scientific << omicron
                 #ifdef DAMPER_FLAG
                 << std::setw(13) << std::setprecision(5) << std::scientific << damper_kick
                 #endif
                 << std::endl;
            #endif
                            
            //
            // Iterations --- Particles
            //
            
            // Damper kick
            #ifdef DAMPER_FLAG
            if(damper_kick != 0)
            {
                for(int m=0; m < particles_number; m++)
                {
                    (*store)(m,XP_INDEX) += damper_kick;
                }
            }
            #endif
            
            for(int m=0; m < particles_number; m++)
            {
                // Faster to check that than to compute invalid particles...
                if((*store)(m,X_INDEX) != LOST || (*store)(m,XP_INDEX) != LOST || (*store)(m,Y_INDEX) != LOST || (*store)(m,YP_INDEX) != LOST)
                {     
                    /*
                        !!! WARNING !!!
                        This is were 'real' physics happens !
                    */
                    
                    // Non linear kicks  
                    store_x2_tmp = (*store)(m,X_INDEX) * (*store)(m,X_INDEX);
                    store_y2_tmp = (*store)(m,Y_INDEX) * (*store)(m,Y_INDEX); 
                    kick_x =   (*store)(m,XP_INDEX) 
                            + iota * (store_x2_tmp - chik * store_y2_tmp)
                            + omicron * ((*store)(m,X_INDEX) * store_x2_tmp - 3 * chik * (*store)(m,X_INDEX) * store_y2_tmp);
                         
                    kick_y =   (*store)(m,YP_INDEX) 
                            - 2 * chik * iota * (*store)(m,X_INDEX) * (*store)(m,Y_INDEX)
                            - omicron * (chik2 * store_y2_tmp * (*store)(m,Y_INDEX)-3 * chik * (*store)(m,Y_INDEX) * store_x2_tmp);
                    
                    
                    // Linear rotations
                    x =  x_cos*(*store)(m,X_INDEX) + x_sin*kick_x;
                    xp= -x_sin*(*store)(m,X_INDEX) + x_cos*kick_x;
                    y =  y_cos*(*store)(m,Y_INDEX) + y_sin*kick_y;
                    yp= -y_sin*(*store)(m,Y_INDEX) + y_cos*kick_y;
               
                    // Check and assign new coordinates
                    if(x >= 1 || xp >= 1 || y >= 1 || yp >= 1)
                    {
                        (*store)(m,X_INDEX) = LOST;
                        (*store)(m,XP_INDEX) = LOST;
                        (*store)(m,Y_INDEX) = LOST;
                        (*store)(m,YP_INDEX) = LOST;
                        n_lost++;
                    }
                    else
                    {
                        (*store)(m,X_INDEX) = x;
                        (*store)(m,XP_INDEX) = xp;
                        (*store)(m,Y_INDEX) = y;
                        (*store)(m,YP_INDEX) = yp;
                    }
                }
            } // Iterations particles
        } // Iterations kicks
        //
        // Tracking and coherent motion
        //
        if(is_tracked)
        {
            for(int i =0; i < n_tracked; i++)
            {
                tracking_ic[i][t][X_INDEX] = (*store)(particles_number-n_tracked+i,X_INDEX);
                tracking_ic[i][t][XP_INDEX] = (*store)(particles_number-n_tracked+i,XP_INDEX);
                tracking_ic[i][t][Y_INDEX] = (*store)(particles_number-n_tracked+i,Y_INDEX);
                tracking_ic[i][t][YP_INDEX] = (*store)(particles_number-n_tracked+i,YP_INDEX);
            }
            if(is_coherent)
            {
                centroid_x = (*store)(0,X_INDEX);
                centroid_xp = (*store)(0,XP_INDEX);
                centroid_y = (*store)(0,Y_INDEX);
                centroid_yp = (*store)(0,YP_INDEX);
                for(int m = 1; m < particles_number; m++)
                {
                    centroid_x += (*store)(m,X_INDEX);
                    centroid_xp += (*store)(m,XP_INDEX);
                    centroid_y += (*store)(m,Y_INDEX);
                    centroid_yp += (*store)(m,YP_INDEX);
                }
                coherent_parameters[t][X_INDEX] = (centroid_x-n_lost*LOST)/particles_number;
                coherent_parameters[t][XP_INDEX] = (centroid_xp-n_lost*LOST)/particles_number;
                coherent_parameters[t][Y_INDEX] = (centroid_y-n_lost*LOST)/particles_number;
                coherent_parameters[t][YP_INDEX] = (centroid_yp-n_lost*LOST)/particles_number;
            }
        }

// Check in case of intermediate processing needed
// Can be optimized with a flag !!!      
#ifdef INTERMEDIATE_TURNS_FLAG
	//if(t == intermediate_turns)
	if(t%intermediate_turns==0)
	{
    	LOG(ITERATOR, Intermediate processing ?);
		flow->setTune(table_tune_h[0][t]);
    	flow->process("./flows/alessandroFlow/dataProcessing/intermediate/");
	}                             
#endif
    } // Iterations turns
    
    // Adjust the position of the islands
    LOG(ITERATOR, Rotation of the coordinates (R^-1));
    phaseRotation(-horizontal_phase_modification, X_INDEX); // Note the minus sign
    phaseRotation(-vertical_phase_modification, Y_INDEX); // Note the minus sign
   
    // Performance estimation
    static_cast<void>(time(&t2));
    Log::getInstance().warn("ITERATOR", "Duration of the iteration process: ", "s", static_cast<int>(t2 - t1));  
    
    // Finish by adding the evolution of the parameters to the graph
    if(with_plots_)
    {
        for(int k = 0; k < n_kicks_; k++)
        {         
	        graph_->draw(table_tune_h[k], "Horizontal tune");
            graph_->draw(table_tune_v[k], "Vertical tune");
            #ifdef DAMPER_FLAG
            if(isDamperSet[k])
            {
                graph_->draw(table_damper[k], "Damper");
            }
            #endif                        
            if(k==n_kicks_-1)
            {
                graph_->draw(table_kappa[k], "Kappa",true);
            }
            else
            {
                graph_->draw(table_kappa[k], "Kappa",false);
            }
        }
    }
    
    // Close the output file
    #ifdef ITERATOR_TO_FILE_FLAG
    file.close();
    #endif
    LOG(ITERATOR 4D, End of the iteration process);
}

void Iterator::iterate6d() {}

inline void
Iterator::phaseRotation(double phase, int index)
{
    if(index != X_PLANE && index != Y_PLANE)
    {
        throw("Invalid index for the phase rotation !");
    }
    int particles_number = part_set->getParticlesNumber();
    double angle = phase * 2 * M_PI;
    double cos_adjust = cos(angle);
    double sin_adjust = sin(angle);
    double x, xp;
    for(int m=0; m < particles_number; m++)
    {
        // Faster to check that than to compute invalid particles...
        if((*store)(m,X_INDEX) != LOST || (*store)(m,XP_INDEX) != LOST)
        {   
            x =  cos_adjust*(*store)(m,index) + sin_adjust*(*store)(m,index+1);
            xp= -sin_adjust*(*store)(m,index) + cos_adjust*(*store)(m,index+1);
                
            (*store)(m,index) = x;
            (*store)(m,index+1) = xp;
            
            if((*store)(m,index) >= 1 || (*store)(m,index+1) >= 1)
            {
                (*store)(m,index) = LOST;
                (*store)(m,index+1) = LOST;
            }
        }
    } // Iterations particles (for phase adjustment)
}

inline double
Iterator::evaluateVariable(CurveParameters& p)
{
    if(p.curve_type == 'l')
    {
	    if(t>20000)
		    return p.start+p.inc*20000;
	    else

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
Iterator::readParameters()
{
    std::string kick_path;
    char curve_info[] = "t";
    #ifdef DEBUG_FLAG
        std::string tmp = "";
    #endif
    LOG(ITERATOR, Reading the input parameters...);        

	//
    // Options for any number of dimensions             
	//
    if(part_set->getDims() == TWO_DIMENSIONS || part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
    {
        LOG(ITERATOR, Reading the input parameters for 2D or 4D or 6D...);
        // File name
        LOG(ITERATOR, Reading file name...);
        this->file_name = "iterations";//XMLInputParser::getInstance().getAttributeTextFromPath("./flows/alessandroFlow/iterator", "file");
        
        // Turns
        LOG(ITERATOR, Reading turns...);
        turns = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/iterator/turns");
        
        // Intermediate turn    
	#ifdef INTERMEDIATE_TURNS_FLAG
        	LOG(ITERATOR, Reading intermediate turns...);
        	intermediate_turns = XMLInputParser::getInstance().getFirstTexti("./flows/alessandroFlow/iterator/intermediateTurns");
        #endif

        // Phase adjustment
        LOG(ITERATOR, Reading phase modification...);
        horizontal_phase_modification = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/iterator/phaseAdjustment/horizontal");
        
        // Kicks
        LOG(ITERATOR, Reading number of kicks...);
        n_kicks_ = static_cast<int>(XMLInputParser::getInstance().getAttribute("./flows/alessandroFlow/iterator","kicks"));
        
        // Damper state variable
	#ifdef DAMPER_FLAG
	         isDamperSet = new bool[n_kicks_];
	 #endif

        // This function read the type variable of each kick
        readKicksTypes();
        
        // With these information we can create the containers for the other parameters
        createContainers();
        
        // Read the other parameters for each kick
        LOG(ITERATOR, Reading the input parameters for 2D or 4D or 6D --- each kick...);
        for(int i=0; i < n_kicks_; i++)
        {
            kick_path = "./flows/alessandroFlow/iterator/kick" + i2s(i); 
            
            // Reading the normalization, depending on the type of the kick
            LOG(ITERATOR, Reading normalization);
            if(kick_types[i] == "SO" || kick_types[i] == "S")
            {
                normalizations[i] = XMLInputParser::getInstance().getFirstTextd(kick_path + "/lambda");
            }
            else if(kick_types[i] == "O")
            {
                normalizations[i] = XMLInputParser::getInstance().getFirstTextd(kick_path + "/tau");
            }
            else if(kick_types[i] == "L")
            {
                normalizations[i] = 0.0;
            }
            
            // Nu horizontal
            LOG(ITERATOR, Reading evolv for nu horizontal);
            curve_info[0] = param_evolution->readEvolvingVariable(i,static_cast<std::string>("nuHorizontal"), hor_tune_);
            if(curve_info[0] == 'c' || curve_info[0] == 's' || curve_info[0] == 'd')
            {
                throw("Horizontal tune cannot be a windowed curve !");
            }
            #ifdef DEBUG_FLAG
            tmp = "Nu horizontal curve type is: ";
            tmp += curve_info;
            LOGV(ITERATOR, tmp);
            #endif
            
            // Kappa
            if(kick_types[i] == "SO")
            {
                LOG(ITERATOR, Reading evolv for kappa);
                curve_info[0] = param_evolution->readEvolvingVariable(i,static_cast<std::string>("kappa"),kappa_);
                if(curve_info[0] == 'c' || curve_info[0] == 's' || curve_info[0] == 'd')
                {
                    throw("Kappa cannot be a windowed curve !");
                }
                #ifdef DEBUG_FLAG
                    tmp = "Kappa horizontal curve type is: ";
                    tmp += curve_info;
                    LOGV(ITERATOR, tmp);
                #endif
            } 
            else if(kick_types[i] == "S" || kick_types[i] == "O")
            {
                LOG(ITERATOR, No kappa for this kick !);
                param_evolution->setToZero(i, kappa_);
            }
            else if(kick_types[i] == "L")
            {
                LOG(ITERATOR, No kappa for this kick !);
                param_evolution->setToZero(i, kappa_);
            }
            else
            {
                throw("Invalid kick type!");
            }
            
            // Is there a damper ?
            #ifdef DAMPER_FLAG
            	std::string tmp_path = "./flows/alessandroFlow/iterator/kick" + i2s(i) + "/damper";
            	isDamperSet[i] = XMLInputParser::getInstance().isFoundFromPath(tmp_path);
            #endif
        
            // Damper
            #ifdef DAMPER_FLAG
            	if(isDamperSet[i])
            	{
                	LOG(ITERATOR, Reading evolv for damper);
                	curve_info[0] = param_evolution->readEvolvingVariable(i,static_cast<std::string>("damper"), damper);     
                	if(curve_info[0] == 'f' || curve_info[0] == 'l' || curve_info[0] == 'p')
                    {
                        throw("Horizontal tune must be a windowed curve !");
                    }
                	#ifdef DEBUG_FLAG
                		tmp = "Damper curve type is: ";
                   		tmp += curve_info;
                   		LOGV(ITERATOR, tmp);
                	#endif 
            }
            #endif
        }
    }
    // Only in the 4D or 6D cases
    if(part_set->getDims() == 4 || part_set->getDims() == 6)
    {
        LOG(ITERATOR, Reading the input parameters for 4D or 6D...);
        // Phase adjustment
        vertical_phase_modification = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/iterator/phaseAdjustment/vertical");

        // Read the other parameters for each kick
        LOG(ITERATOR, Reading the input parameters for 4D --- each kick...);
        for(int i=0; i < n_kicks_; i++)
        {
            kick_path = "./flows/alessandroFlow/iterator/kick" + i2s(i); 
            // Chi
            LOG(ITERATOR, Reading chi);
            if(kick_types[i] == "L")
            {
                chi[i] = 0;
            }
            else
            {
                chi[i] = XMLInputParser::getInstance().getFirstTextd(kick_path + "/chi");
            }
            
            // Nu vertical
            LOG(ITERATOR, Reading evolv for nu vertical);
            curve_info[0] = param_evolution->readEvolvingVariable(i,static_cast<std::string>("nuVertical"), ver_tune_);
            if(curve_info[0] == 'c' || curve_info[0] == 's' || curve_info[0] == 'd')
            {
                throw("Vertical tune cannot be a windowed curve !");
            }
            #ifdef DEBUG_FLAG
            tmp = "Nu vertical curve type is: ";
            tmp += curve_info;
            LOGV(ITERATOR, tmp);
            #endif
        }
    }
    // Only in the 6D cases
    if(part_set->getDims() == 6)
    {
        // Nu longitudinal
        LOG(ITERATOR, Reading evolv for nu longitudinal);
        param_evolution->readEvolvingVariable(0,static_cast<std::string>("nuLongitudinal"), long_tune);
        // Horizontal chromaticity
        LOG(ITERATOR, Reading evolv for hor_chroma);
        param_evolution->readEvolvingVariable(0,static_cast<std::string>("chromaticityHorizontal"), hor_chromaticity);
        // Vertical chromaticity
        LOG(ITERATOR, Reading evolv for ver_chroma);
        param_evolution->readEvolvingVariable(0,static_cast<std::string>("chromaticityVertical"), ver_chromaticity);
        
        // Read the other parameters for each kick
        for(int i=0; i < n_kicks_; i++)
        {
            // Dispersion
            LOG(ITERATOR, Reading evolv for dispersion);
            param_evolution->readEvolvingVariable(i,static_cast<std::string>("dispersion"), dispersion);
            // Dispersion prime
            LOG(ITERATOR, Reading evolv for dispersion prime);
            param_evolution->readEvolvingVariable(i,static_cast<std::string>("dispersionPrime"), dispersion_prime);
        }
    }
    LOG(ITERATOR, Input parameters read.);
}

void
Iterator::createContainers()
{
    LOG(ITERATOR, Creating containers.);
	// 2 or 4 or 6 dimensions
	if(part_set->getDims() == TWO_DIMENSIONS || part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
	{
	    // Physical parameters
	    hor_tune_ = new CurveParameters[n_kicks_];
	    kappa_ = new CurveParameters[n_kicks_];
	    normalizations.resize(n_kicks_);
	    #ifdef DAMPER_FLAG
	    damper = new CurveParameters[n_kicks_];
	   #endif

	    // Evolution tables
		table_tune_h = new double*[n_kicks_];
		table_kappa = new double*[n_kicks_];
		#ifdef DAMPER_FLAG
		table_damper = new double*[n_kicks_];
		#endif
		for(int u = 0; u < n_kicks_; u++)
		{
			table_tune_h[u] = new double[turns];
			table_kappa[u] = new double[turns];
			#ifdef DAMPER_FLAG
            table_damper[u] = new double[turns];
			#endif
		}
	}
	// 4 or 6 dimensions
	if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
	{
	    // Physical parameters
	    ver_tune_ = new CurveParameters[n_kicks_];
	    chi.resize(n_kicks_);
	    
	    // Evolution tables
		table_tune_v = new double*[n_kicks_];
		for(int u = 0; u < n_kicks_; u++)
		{
			table_tune_v[u] = new double[turns];
		}
	}
	// 6 dimensions
	if(part_set->getDims() == SIX_DIMENSIONS)
	{
	    // Curve parameters
	    long_tune = new CurveParameters[1];
	    hor_chromaticity = new CurveParameters[1];
	    ver_chromaticity = new CurveParameters[1];
	    dispersion = new CurveParameters[n_kicks_];
	    dispersion_prime = new CurveParameters[n_kicks_];
	    
	    // Evolution tables
		table_tune_l = new double*[1];
		table_tune_l[0] = new double[turns];
	}
}

void 
Iterator::readKicksTypes()
{
    // Create the string table
  kick_types.resize(n_kicks_);
    
    // Read each kick
    LOG(ITERATOR, Reading kicks types...);
    for(int i=0; i < n_kicks_; i++)
    {
        std::string path = "./flows/alessandroFlow/iterator/kick" + i2s(i);
        if(XMLInputParser::getInstance().isFoundFromPath(path))
        {
            // Kick types
            kick_types[i] = XMLInputParser::getInstance().getAttributeTextFromPath(path, "type");
            if(kick_types[i] != static_cast<std::string>("SO") && kick_types[i] != static_cast<std::string>("S") && kick_types[i] != static_cast<std::string>("O") && kick_types[i] != static_cast<std::string>("L"))
            {
                throw("Invalid kick type !");
            }
            if(i == 0 && kick_types[i] == static_cast<std::string>("L"))
            {
                throw("The first kick as to be non linear !");
            }
            // Log
            #ifdef DEBUG_FLAG
                std::string tmp = "Kick " + i2s(i) + " type is: " + kick_types[i];
                LOGV(ITERATOR, tmp);
            #endif
        }
        else
        {
            throw("Something strange happened while reading the kicks...");
        }
    }
}

int
Iterator::getTurns()
{
    return turns;
}          
