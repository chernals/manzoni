#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <distributions/ParticlesDistribution.h>
#include <utils/XMLInputParser.h>

#define X_INDEX 0
#define XP_INDEX 1
#define Y_INDEX 2
#define YP_INDEX 3
#define Z_INDEX 4
#define ZP_INDEX 5

class ParticlesSet : public ParticlesDistribution
{
  private:
    //
    // Data
    //
	//
	// Shared
	int turns;
	bool is_tracking;
	int n_tracking;
	double*** tracking_ic;
	double** coherent_parameters;
	bool is_coherent;
	//
	// Gaussian
	double _gaussian_hor_sigma_x;
	double _gaussian_hor_sigma_xp;
	double _gaussian_hor_centroid_x;
	double _gaussian_hor_centroid_xp;
	unsigned int _gaussian_hor_density;
	double _gaussian_ver_sigma_y;
	double _gaussian_ver_sigma_yp;
	double _gaussian_ver_centroid_y;
	double _gaussian_ver_centroid_yp;
	unsigned int _gaussian_ver_density;
    double _gaussian_long_sigma_z;
	double _gaussian_long_sigma_zp;
	double _gaussian_long_centroid_z;
	double _gaussian_long_centroid_zp;
	unsigned int _gaussian_long_density;
	//
	// Uniform
	double _initial_centroid_horizontal;
	double _initial_centroid_p_horizontal;
	double _initial_centroid_vertical;
	double _initial_centroid_p_vertical;
	double _initial_centroid_longitudinal;
	double _initial_centroid_p_longitudinal;
	double _radius_limit_horizontal;
	double _radius_limit_vertical;
	double _radius_limit_longitudinal;
	unsigned int _radial_density_horizontal;
	unsigned int _radial_density_vertical;
	unsigned int _radial_density_longitudinal;
	unsigned int _angular_density_horizontal;
	unsigned int _angular_density_vertical;
	unsigned int _angular_density_longitudinal;
	double _radial_increment_horizontal;
	double _radial_increment_vertical;
	double _radial_increment_longitudinal;
	double _angular_increment_horizontal;
	double _angular_increment_vertical;
	double _angular_increment_longitudinal;
	std::string _distr_file;
	//
	// Methods
	//
	void readParameters();
	void computeDerivedParameters();
	void createParticlesGaussian();
	void createParticlesUniform();
	void createParticlesFile();
	void addTrack2Store();

  public:
	ParticlesSet();
	virtual ~ParticlesSet();
    void setInitialDistribution();
    void setTracking();
    void setCoherent();
	//
	// Getters --- Setters
	//
	store_t* getStore() const;
	double*** getTracking() const;
	int getParticlesNumber();
	int getDims();
    bool isTracking();
    int getTrackingNumber();
    double** getCoherent() const;
    bool isCoherent();
    
};

#endif

