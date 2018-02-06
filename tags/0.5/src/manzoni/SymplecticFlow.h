#ifndef SYMPLECTIC_FLOW_H
#define SYMPLECTIC_FLOW_H

#include <vector>
#include <string>
                    
#include <definitions.h>
#include <Flow.h>

class SymplecticFlow : public Flow
{
	private:

	public:
		SymplecticFlow();
		~SymplecticFlow();
		void run();
};

#endif
