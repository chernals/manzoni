#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <unistd.h>
#include <stdio.h>

#include <iostream>
#include <map>

#include <utils/Singleton.h>
#include <utils/Logger.h>

typedef std::map<char, char*> ArgsMap;

class Arguments
{
	private:
		ArgsMap map;
		char* exec_name;

	public:
		void init(const char* const, int, char**);
		const char* get_argument(char);
		const char* get_exec_name();
		void show();

};

typedef Singleton<Arguments> Arg;

#endif
