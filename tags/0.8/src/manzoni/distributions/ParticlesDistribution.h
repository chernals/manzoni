#ifndef PARTICLES_DISTRIBUTION_H
#define PARTICLES_DISTRIBUTION_H

#include <types/types.h>
#include <utils/Logger.h>
#include <utils/TypeConversions.h>

class ParticlesDistribution
{
  protected:
		// This variable is an array-like container containing the particles set
		store_t* store;
		// General distribution information
    std::string distr_type;
    unsigned short int distr_dims;
		// This flag is true if the distribution is defined
		bool isSet;
		unsigned int particlesNumber;

  public:
		virtual store_t* getStore() const = 0;
		virtual int getParticlesNumber() = 0;
		virtual int getDims() = 0;
    virtual bool isTracking() = 0;
    virtual int getTrackingNumber() = 0;
    virtual double** getCoherent() const = 0;
    virtual bool isCoherent() = 0;
	  virtual double*** getTracking() const = 0;
	  virtual void setInitialDistribution() = 0;
    virtual void setTracking() = 0;
    virtual void setCoherent() = 0;
		virtual ~ParticlesDistribution() {};

	  #ifdef DEBUG_FLAG
			void showStoreInfo(bool with_store=false);
		#endif
	
};

#endif
