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

#include <main.h>

int
main(int argc, char** argv)
{
	try 
	{       		
		// Parsing of the arguments
		Arg::getInstance().init("--hdi:l:", argc, argv);    
		
		// Define the various signals
		signalsDefinition();

		// Create the logger object so it is initialized here
		Log::getInstance();
        
    // Check if help argument is defined
    if(Arg::getInstance().get_argument('h') != NULL)
    {
			helpMessage();
			return EXIT_SUCCESS;
		}                                   
		
		// All set
		startMessage();
		  
		// Launch Manzoni
		LOG(MAIN,Launching Manzoni...);
    Manzoni manzoni;
    LOG(MAIN, Ending Manzoni.);
  
		return EXIT_SUCCESS;
	}  // Try
	catch(const char* excep)
	{
		ERRORV(MAIN CATCH, excep);
		Log::getInstance().perror_();
		ERROR(MAIN CATCH, This exception terminates Manzoni !);
		return EXIT_FAILURE;
	} // Catch

}  

void
signalsDefinition()
{
    // SIGINT signal
    // The cast if to avoid SIG_ERR old style cast warning
	if(signal(SIGINT, Signals::signals_sigint) == reinterpret_cast<void(*)(int)>(-1))
    {
    	throw("SIGINT signal definition failed.");
	}
}

void
version()
{
    #ifdef VERSION
    double v = VERSION;
	std::cout << "     ---    Version " << d2s(v) << std::endl;
    #else
	std::cout << "     ---    Version not set" << std::endl; 
    #endif
}

void
startMessage()
{
	std::cout << "========================================================================" << std::endl;
    std::cout << " Manzoni" << std::endl;
    version();
    std::cout << "========================================================================" << std::endl;
	if(Arg::getInstance().get_argument('d') != NULL)
	{
       		fprintf(stdout, " Syntax : %s -i%s -l%s -d\n", Arg::getInstance().get_exec_name(), Arg::getInstance().get_argument('i'), Arg::getInstance().get_argument('l'));
	}
	else
	{
        	fprintf(stdout, " Syntax : %s -i%s -l%s \n", Arg::getInstance().get_exec_name(), Arg::getInstance().get_argument('i'), Arg::getInstance().get_argument('l'));
    }
    compilationOptions();
    fprintf(stdout, "========================================================================\n");
	Arg::getInstance().show();
}
 
void
helpMessage()
{
  	std::cout << "========================================================================" << std::endl;
	std::cout << " Manzoni" << std::endl;
        version();
		 std::cout << "========================================================================" << std::endl;
		fprintf(stdout, " Usage : %s  [-f flow_id] [-i input_file] [-l log_level] [-d] [-h]\n", Arg::getInstance().get_exec_name());
		fprintf(stdout, "   * -l log_level : 0 (errors only), 1 (warnings) or 2 (verbose log)\n");
		fprintf(stdout, "   * -d           : debug mode, log messages on screen instead of on file\n");
		fprintf(stdout, "   * -h           : help, print this message\n");
		fprintf(stdout, "   * -i input_file: the input file\n");
		fprintf(stdout, "\n");
		fprintf(stdout, "   NOTE: Verbose log is allowed only if Manzoni was compiled with DEBUG_FLAG set\n");
         std::cout << "========================================================================" << std::endl;
		fprintf(stdout, " Flows :\n");
		fprintf(stdout, "   * Hénon: iteration of the generalized Hénon map \n"); 
		fprintf(stdout, "   * Symplectic: symplectic integration of Hamiltonian models \n"); 
		compilationOptions();
		std::cout << "========================================================================" << std::endl;
}

void
compilationOptions()
{
	std::cout << "========================================================================" << std::endl;
  	fprintf(stdout, " Compilation options :\n");
  #ifdef DEBUG_FLAG
		fprintf(stdout, "   * DEBUG_FLAG is SET (verbose log allowed)\n");
	#else
	  fprintf(stdout, "   * DEBUG_FLAG is NOT SET (verbose log not allowed)\n");
	#endif    
	#ifdef ITERATOR_TO_FILE_FLAG
		fprintf(stdout, "   * ITERATOR_TO_FILE_FLAG is SET (iterator will output in a file)\n");
	#else
		fprintf(stdout, "   * ITERATOR_TO_FILE_FLAG is NOT SET (iterator will not output in a file)\n");
	#endif               
	#ifdef INTERMEDIATE_TURNS_FLAG
		fprintf(stdout, "   * INTERMEDIATE_TURNS_FLAG is SET (intermediate data processing allowed)\n");
	#else
    fprintf(stdout, "   * INTERMEDIATE_TURNS_FLAG is NOT SET (intermediate data processing not allowed)\n");
	#endif
	#ifdef DAMPER_FLAG
	  fprintf(stdout, "   * DAMPER_FLAG is SET (dampers definition and execution are allowed)\n");
	#else
	  fprintf(stdout, "   * DAMPER_FLAG is NOT SET (damper definition and execution are not allowed)\n");
	#endif
	#ifdef BZ_DEBUG
    fprintf(stdout, "   * BZ_DEBUG is SET (blitz++ library used in debug mode, EXTREMELY SLOW !)\n"); 
  #else   
  	fprintf(stdout, "   * BZ_DEBUG is NOT SET (blitz++ library not used in debug mode)\n"); 
	#endif
	#ifdef XML_PARSER_DEBUG
  	fprintf(stdout, "   * XML_PARSER_DEBUG is SET (XML parser will verbose output XML processing)\n"); 
  #else
  	fprintf(stdout, "   * XML_PARSER_DEBUG is NOT SET (XML parser will not output XML processing)\n"); 
	#endif
}               
                          
