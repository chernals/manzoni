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

#ifndef STORE_T_H
#define STORE_T_H

class store_t {
  private:
    std::vector<std::vector<double> > vec;
    double p_, d_;
  public:
  store_t(unsigned int p, int d) : p_(p), d_(d) { 
    vec.resize(static_cast<size_t>(p));
    for(unsigned int i = 0; i < p; i++)
      vec.at(i).resize(static_cast<size_t>(d));
  }
  double& operator()(unsigned int a, unsigned int b) { 
    return vec.at(static_cast<size_t>(a)).at(static_cast<size_t>(b)); 
  }
  store_t& operator=(const double& a) {
    for(unsigned int i = 0; i < p_; i++)
      for(int j = 0; j < d_; j++)
        vec.at(i).at(static_cast<size_t>(j)) = a;
    return *this;
  }
};

#endif

