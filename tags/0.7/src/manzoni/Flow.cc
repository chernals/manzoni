#include <Flow.h>

Flow::Flow()
{
    LOG(FLOW, Instanciating flow);
    iterator = NULL;
    particles_set = NULL;
}

Flow::~Flow()
{
    LOG(FLOW, Destructor called);
	if(particles_set != NULL)
	{
	    delete particles_set;
    }
    if(iterator != NULL)
    {
        delete iterator;
    }
}

void
Flow::run()
{
    LOG(FLOW, Running flow.);
}

void
Flow::process(std::string path)
{
	LOG(FLOW, Checking if processing is needed);
	std::string xpath = path + "*";
	if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	{
	    LOG(FLOW, Processing needed);
        // Processing of DistrHisto
        xpath = path + "distrHist";
	    if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	    {
	        LOG(FLOW, Processing with ProcDistrHist found);
	        ProcDistrHisto histo(particles_set, xpath);
	        histo.draw();
	        histo.save();
	    }
	    
	    // Processing of DistrProfile
        xpath = path + "distrProfile";
	    if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	    {
	        LOG(FLOW, Processing with ProcDistrProfile found);
	        ProcDistrProfile pro(particles_set, xpath);
	        pro.draw();
	        pro.save();
	    }
	    
	    // Processing of Islands
        xpath = path + "islands";
	    if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	    {
	        LOG(FLOW, Processing with ProcIslands found);	         
	        ProcIslands isl(particles_set,xpath);
	    }
	    
	    // Processing of the tracking
	    xpath = path + "tracking";
	    if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	    {
	        LOG(FLOW, Processing with ProcTracking found);
	        ProcTracking track(particles_set,xpath);
	    }
	}
	else
	{
	    LOG(FLOW, No processing for this time);
	}
	LOG(FLOW, End of processing.);
}
