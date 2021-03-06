#ifndef TYPES_H
#define TYPES_H

#include <blitz/array.h>

#include <definitions.h>

// We typedef the storage type to allow change of implementation
// See blitz doc for more information on the blitz::Array
typedef blitz::Array<double,TWO_DIMENSIONS> store_t;

// Include other types
#include <types/CurveParameters.h>

#endif

