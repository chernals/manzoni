#ifndef TRACKING_DISTRIBUTION_H
#define TRACKING_DISTRIBUTION_H

#include <math.h>

#include <gsl/gsl_rng.h>

#include <iostream>

#include <distributions/ParticlesDistribution.h>
#include <XMLInputParser.h>

class TrackingDistribution : public ParticlesDistribution
{
  private:
		//
		// Data
		//
		// This private variable is an array-like container containing the particles set
		store_t* store;
		int turns;
		
	  public:
		TrackingDistribution(unsigned int);
		~TrackingDistribution();
	  void setInitialDistribution();
		//
		// Getters --- Setters
		//
		store_t* getStore() const;
		void setTracking();
    void setCoherent();
		double*** getTracking() const;
   	double** getCoherent() const;
		int getParticlesNumber();
		int getDims();
		bool isTracking();
    int getTrackingNumber();
    bool isCoherent();
		
};

#endif
