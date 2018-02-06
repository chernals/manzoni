#ifndef AA_DISTRIBUTION_H
#define AA_DISTRIBUTION_H

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>

#define ACTION_INDEX 1
#define ANGLE_INDEX 0

class AADistribution : public ParticlesDistribution
{
  private:
		//
		// Data
		//
		// This private variable is an array-like container containing the particles set
		store_t* store;

		double _uniform_radius;
		bool _angular_extension;
		double _radial_offset;
		double _uniform_radial_density;
		double _uniform_angular_density;
		double _radial_increment_horizontal;
		double _angular_increment_horizontal;
		std::string _distr_file;

		//
		// Methods
		//
		void readParameters();
		void createParticlesUniform();
		void createParticlesGaussian();
		void createParticlesFile();

	  public:
		AADistribution();
		~AADistribution();
	  void setInitialDistribution();
		//
		// Getters --- Setters
		//
		store_t* getStore() const;
		void setTracking();
    void setCoherent();
		double*** getTracking() const;
   double** getCoherent() const;
		unsigned int getParticlesNumber();
		int getDims();
		bool isTracking();
    unsigned int getTrackingNumber();
    bool isCoherent();
		
};

#endif
