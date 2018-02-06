#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>

#include <string>
#include <iostream>

#include <utils/Singleton.h>
#include <utils/Arguments.h>

#define LOGG 2
#define WARNN 1
#define ERRORR 0
#define BAD_VALUE -1

#if defined DEBUG_FLAG
#define LOG(A,B) Log::getInstance().log(#A, #B)
#define LOGV(A,B) Log::getInstance().log(#A,B)
#else
#define LOG(A,B) ;
#define LOGV(A,B) ;
#endif

#define WARN(A,B) Log::getInstance().warn(#A, #B)
#define WARNV(A,B) Log::getInstance().warn(#A, B)

#define ERROR(A,B) Log::getInstance().error(#A, #B)
#define ERRORV(A,B) Log::getInstance().error(#A, B)

class Logger
{
	private:
		int log_level;
		bool debug_mode;
		const char* exec_name;
		time_t old_time;
		char* construct_string(const char* const, const char* const, int, const char* const, const char* const);
		void redirect_to_file() const;

	public:
		Logger();
		~Logger() {};
		void log(const char* const, const char* const, const char* const = "", int value = BAD_VALUE);
		void log(const char* const, std::string, const char* const = "", int value = BAD_VALUE);
		void warn(const char* const, const char* const, const char* const = "", int value = BAD_VALUE);
		void warn(const char* const, std::string, const char* const = "", int value = BAD_VALUE);
		void error(const char* const, const char* const, const char* const = "", int value = BAD_VALUE);
		void error(const char* const, std::string, const char* const = "", int value = BAD_VALUE);
		void perror_() const;

};

typedef Singleton<Logger> Log;

#endif

