#include <utils/Arguments.h>

void 
Arguments::init(const char* const allowed_args, int argc, char** argv)
{ 
	// Executable name
    exec_name = argv[0];

	// No automatic error message
    opterr = 0;
    
	// Parse the arguments list  
	int opt;
    while((opt = getopt(argc, argv, allowed_args)) != -1)
	{
		if(opterr == 1)
			std::cout << "Opt: " << static_cast<char>(opt) << std::endl;
		
		switch(opt)
		{
			// Unrecognized option
			case '?':
         		break;  
			case ':':
						break;
 			// Other options
			default:
         		if (optarg != NULL)
                 	map[static_cast<char>(opt)] = optarg; 
         		else
         		{
         		    const char* tmp = "1";
                 	map[static_cast<char>(opt)] = const_cast<char*>(tmp);
                }
				break;
     	}       
	}       
}

const char*
Arguments::get_argument(char a)
{
        if(map.count(a) == 0)
        {
					return NULL;
        }
        return map[a]; 
}

const char*
Arguments::get_exec_name()
{        
        return exec_name;
}

void Arguments::show()
{        
    #ifdef DEBUG_FLAG
        for(ArgsMap::const_iterator it=map.begin(); it!=map.end(); ++it)
        {       
            char* string = new char[BUFFER_SIZE];
            sprintf(string, "Argument -- %c - %s", it->first, it->second);
            LOGV(ARGUMENTS, string);
        }
    #endif       
}
