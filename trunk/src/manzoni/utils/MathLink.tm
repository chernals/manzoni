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

#include <utils/MathLink.h>

#define STRNAME1(n) #n
#define STRNAME(n) STRNAME1(n)

MathLink::MathLink()
{           
  #ifdef MATH_ENABLED
    LOG(MATH LINK, Instanciating MathLink);
    env = MLInitialize(NULL);
    if(env == NULL) return;
    std::cout << "-linkmode launch -linkname '"STRNAME(MATH_KERNEL_PATH)" -mathlink'" << std::endl << std::flush;
    lp = MLOpenString(env,"-linkmode launch -linkname '"STRNAME(MATH_KERNEL_PATH)" -mathlink'",&pkt);
    if(lp == NULL) return;      
#endif 
}

MathLink::~MathLink()
{
#ifdef MATH_ENABLED
  LOG(MATH LINK, Destroying MathLink);
  MLClose(lp);
  MLDeinitialize(env);   
  
#endif
}

std::string
MathLink::execute(std::string str)
{
#ifdef MATH_ENABLED
  // Send for execution
  MLPutFunction(lp, "EvaluatePacket", 1);
  MLPutFunction(lp, "ToString", 1);
  MLPutFunction(lp, "ToExpression", 1);
  MLPutString(lp,str.c_str());
  MLEndPacket(lp);
  
  // Loop to retrieve the correct packet
  while ((pkt = MLNextPacket(lp),pkt) && pkt != RETURNPKT) {
    MLNewPacket(lp);
    if (MLError(lp)) fprintf(stdout, "Error!\n");
  }

  // Get the string
  const char * string;
  if(! MLGetString(lp, &string)) return "";
  std::string tmp(string);
  MLReleaseString(lp, string);
  
  LOG(MATH LINK, Result of Mathematica execution:);
  LOGV(MATH LINK, tmp);
  
  return tmp;
#endif
#ifndef MATH_ENABLED
  return str;
#endif
}
