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

#ifndef TYPES_H
#define TYPES_H

#include <blitz.h>
#include <definitions.h>

// We typedef the storage type to allow change of implementation
// See blitz doc for more information on the blitz::Array
typedef blitz::Array<double,TWO_DIMENSIONS> store_t;

// Include other types
#include <types/CurveParameters.h>

#endif

