#ifndef SYMPLECTIC_FLOW_H
#define SYMPLECTIC_FLOW_H

#include <vector>
#include <string>
                    
#include <definitions.h>
#include <Integrator.h>
#include <Flow.h>

class SymplecticFlow : public Flow
{
	private:
		Integrator* integrator;

	public:
		SymplecticFlow();
		~SymplecticFlow();
		void run();
};

#endif
