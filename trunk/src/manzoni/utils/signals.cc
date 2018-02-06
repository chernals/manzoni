#include <utils/signals.h>

extern bool stop_flag;

void 
Signals::signals_sigint(int) 
{
	WARN("SIGNAL", "SIGINT signal catched.");
	throw("Abord on SIGINT");
}
