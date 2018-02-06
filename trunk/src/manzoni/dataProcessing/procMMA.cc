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

#include <dataProcessing/procMMA.h>

ProcMMA::ProcMMA(ParticlesDistribution* s, std::string path) :
  xpath(path),
  set(s)
{
  readParameters();
  LOG(PROC_MATHEMATICA, Executing.);
  ML::getInstance().execute(constructData());
  ML::getInstance().execute("SetDirectory[\"/Users/cedric/Documents/CERN/Manzoni/SVN/trunk/src/manzoni\"]");
  std::string str("Import[\""+file+"\"]");
  ML::getInstance().execute(str.c_str());
}

ProcMMA::~ProcMMA()
{
}

std::string
ProcMMA::constructData()
{
  std::string str = "data={";
  for(unsigned int i = 0; i < set->getParticlesNumber(); i++)
  {
    str+="{";
    for(unsigned int j = 0; j < static_cast<unsigned int>(set->getDims()); j++)
    {
      str+=d2s(static_cast<double>((*(set->getStore()))(i,j)));
      if(j<static_cast<unsigned int>(set->getDims())-1)
        str+=",";
    }
    str+="}";
    if(i<((set->getParticlesNumber())-1))
      str+=",";
  }
  str+="};\"Data constructed.\"";
  return str;
}
  
void
ProcMMA::readParameters()
{
    LOG(PROC_MATHEMATICA, Reading parameters...);
    file = XMLInputParser::getInstance().getAttributeTextFromPath(xpath, "file");
}

