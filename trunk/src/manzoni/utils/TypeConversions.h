#ifndef TYPE_CONVERSION_H
#define TYPE_CONVERSION_H

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

// Because those are inline functions !
#include <utils/TypeConversions.cc>

// Some very short inline functions
template <class T> inline T s2n(std::string const&);
inline std::string i2s(int const&);
inline std::string d2s(double const&);

#endif
