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

#ifndef MATHLINK_H
#define MATHLINK_H   
          
#include <string>
#include <mathlink.h>
#include <utils/Singleton.h>
#include <utils/Arguments.h>

class MathLink
{
private:
    MLINK lp;
    MLEnvironment env;
    int pkt;

public:
    MathLink();   
    ~MathLink();
    std::string execute(std::string);
};

typedef Singleton<MathLink> ML;

#endif
