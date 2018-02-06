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

#include <Manzoni.h>

Manzoni::Manzoni()
{           
    LOG(MANZONI, Instanciating manzoni);
    flow = NULL;
    try // Select and create the flow                              
    {
			// Get flow from the input
			std::string flow_name = XMLInputParser::getInstance().getFlowName();
			if(flow_name == "henon")
			{
				flow_id = HENON_FLOW;
			}
			else if(flow_name == "symplectic")
			{
				flow_id = SYMPLECTIC_FLOW;
			}
			else
			{
				throw("Invalid 'flow-id' detected in the input file !");
			}
 			createFlow();
		} // Try
	catch(const char* excep)
	{
		ERRORV(MANZONI FLOW CREATION CATCH, excep);
		// Cleaning actions
		if(flow != NULL)
		{
		    delete flow;
		}
		//
		ERROR(MANZONI FLOW CREATION CATCH, Forwarding...);
		throw("Exception catched in MANZONI FLOW CREATION");
	} // Catch
	
    try // Instanciate the XMLInputParser and set the simulation name                
    {
        sim_name = XMLInputParser::getInstance().getSimulationName();
        LOG(MANZONI, Simulation name:);
        LOGV(MANZONI, sim_name);
    } // Try
    catch(const char* excep)
    {
        ERRORV(MANZONI INPUT CATCH, excep);
        ERROR(MANZONI INPUT CATCH, Forwarding...);
        throw("Exception catched in MANZONI INPUT");
    } // Catch
    
	// Run the flow
	flow->run();
	
	LOG(MANZONI, End of the run of the flow.);
}

Manzoni::~Manzoni()
{
    LOG(MANZONI, Destructor called);
    if(flow != NULL)
		{
	    delete flow;
    }
}

void
Manzoni::createFlow()
{
    #ifdef DEBUG_FLAG
	    char* string = new char[BUFFER_SIZE];
	    sprintf(string, "Creating the flow %ld", static_cast<long>(flow_id));
	    LOGV(MANZONI, string);
	    delete string;
    #endif
	
	switch(flow_id)
	{
	    case HENON_FLOW:
	        flow = new HenonFlow();
	        break;
					case SYMPLECTIC_FLOW:
				flow = new SymplecticFlow();
				break;
			default:
	        ERROR(MANZONI, Wrong flow selection !);
	        throw("Wrong flow selection !");   
	        break;
	}
}
