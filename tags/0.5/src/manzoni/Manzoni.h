#ifndef MANZONI_H
#define MANZONI_H

#include <Flow.h>
#include <AlessandroFlow.h>
#include <SymplecticFlow.h>
#include <utils/Logger.h>

class Manzoni
{
	private:
	    //
	    // Data
	    //
	    std::string sim_name;
		unsigned short int flow_id;
		Flow* flow;
		//
	    // Methods
	    //
	    void createFlow();

	public:
		Manzoni();
		~Manzoni();
};

#endif
