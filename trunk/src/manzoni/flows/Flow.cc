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

#include <flows/Flow.h>

Flow::Flow()
{
    iterator = NULL;
    particles_set = NULL;
  	histo_counter = 0;
}

Flow::~Flow()
{
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
Flow::setTune(double t)
{
	tunes.push_back(t);
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
	        ProcDistrHisto histo(particles_set, xpath,histo_counter,tunes);
	        histo.draw();
	        histo.save();
			histo_counter++;
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
      
      // Processing of MMA
      xpath = path + "mathematica";
      if(XMLInputParser::getInstance().isFoundFromPath(xpath))
      {
        LOG(FLOW, Processing with Mathematica found);
        ProcMMA mma(particles_set, xpath);
      }
	    
	    // Processing of Islands
        xpath = path + "islands";
	    if(XMLInputParser::getInstance().isFoundFromPath(xpath))
	    {
	        LOG(FLOW, Processing with ProcIslands found);	         
	        ProcIslands* tmp = new ProcIslands(particles_set,xpath);
					delete tmp;
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
