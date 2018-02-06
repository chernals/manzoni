#include <Manzoni.h>

Manzoni::Manzoni()
{           
    LOG(MANZONI, Instanciating manzoni);
    flow = NULL;
    try // Select and create the flow                              
    {
        // Flow selection
 	    const char* f = Arg::getInstance().get_argument('f');
 	    if(f != NULL)
 	    {
 		    flow_id = static_cast<unsigned short int>(strtol(f, NULL, 0));
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
	    case ALESSANDRO_FLOW:
	        flow = new AlessandroFlow();
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
