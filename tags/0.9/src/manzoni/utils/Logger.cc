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

#include <utils/Logger.h>

Logger::Logger()
{                   
	// Check the log level
	const char* lvl = Arg::getInstance().get_argument('l');
	if(lvl != NULL)               
	{
		log_level = strtol(lvl, NULL, 0);
    }
    else
    {
        log_level = 0;
    }
    if(log_level < 0 || log_level > 2)
    {
        log_level = 2;
    }        
	
	// Check the debug mode
	const char* debug = Arg::getInstance().get_argument('d');
	if(debug != NULL)
	{
		debug_mode = static_cast<bool>(strtol(debug, NULL, 0));
    }
  else
    {
        debug_mode = false;
    }
    
	// Take the executable name
	exec_name = Arg::getInstance().get_exec_name();
            
	// Redirect the std output to a file
	if(!debug_mode)
	{
		redirect_to_file();
    }

	// For the timestamp
	old_time = time(NULL);
}

void
Logger::log(const char* const prefix, const char* const message, const char* const suffix, int value)
{
	if(log_level >= LOGG)
	{
		fprintf(stdout, "%s", construct_string(prefix, message, value, "LOG", suffix));
	}
}

void
Logger::log(const char* const prefix, std::string string, const char* const suffix, int value)
{
	log(prefix, string.c_str(), suffix, value);
}

void
Logger::warn(const char* const prefix, const char* const message, const char* const suffix, int  value)
{
	if(log_level >= WARNN)
	{
		fprintf(stdout, "%s", construct_string(prefix, message, value, "WARNING", suffix));
	}
}

void
Logger::warn(const char* const prefix, std::string string, const char* const suffix, int value)
{
	warn(prefix, string.c_str(), suffix, value);
}

void
Logger::error(const char* const prefix, const char* const message, const char* const suffix, int value)
{
	if(log_level >= ERRORR)   
	{
		fprintf(stderr, "%s", construct_string(prefix, message, value, "ERROR", suffix));
	}
}

void
Logger::error(const char* const prefix, std::string string, const char* const suffix, int value) 
{
	error(prefix, string.c_str(), suffix, value);
}

char*
Logger::construct_string(const char* const prefix, const char* const message, int value, const char* const type, const char* const suffix)
{
	// Pour gerer le timestamp
	time_t t = time(NULL);
	
	// Construct the string
	char* string = new char[BUFFER_SIZE];
	if(t != old_time)
	{
		if(value == BAD_VALUE)
		{
		        sprintf(string, "=======================   %s%s:: %ld : %s : %s: %s\n", ctime(&t), exec_name, static_cast<long>(getpid()), type, prefix, message);
		}
		else
		{    
		        sprintf(string, "=======================   %s%s:: %ld : %s : %s: %s [%s] [%i]\n", ctime(&t), exec_name, static_cast<long>(getpid()), type, prefix, message, suffix, value); 
		}
	}
	else 
	{  
		if(value == BAD_VALUE)
		{
		        sprintf(string, "%s:: %ld : %s : %s: %s\n", exec_name, static_cast<long>(getpid()), type, prefix, message);
		}
		else     
		{
			sprintf(string, "%s:: %ld : %s : %s: %s [%s] [%i]\n", exec_name, static_cast<long>(getpid()), type, prefix, message, suffix, value); 
		}
	}       
	
	// Pour gerer le timestamp
	old_time = time(NULL);
	
	return string; 
}

void
Logger::redirect_to_file() const
{ 
   	mkdir("logs", 0777);
    char file_name[] = "logs/manzoni.log";    
    // See a good system C book for help 
    int fd; 
    if((fd = open(file_name, O_WRONLY | O_CREAT | O_APPEND, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == 0)
    {
        throw errno;
    }
    close(STDOUT_FILENO);
    if(dup(fd) < 0)
    {
        throw("Error while redirecting the std output to file");
    }
}

void
Logger::perror_() const
{
	char* string = new char[BUFFER_SIZE];
  sprintf(string, "%s:: %ld : %s", exec_name, static_cast<long>(getpid()), "ERROR : MAIN CATCH: Now logging perror():");
  perror(string);
}

