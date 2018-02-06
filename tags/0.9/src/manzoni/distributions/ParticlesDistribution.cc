#include <distributions/ParticlesDistribution.h>

#ifdef DEBUG_FLAG
void
ParticlesDistribution::showStoreInfo(bool with_store)
{
    LOG(PARTICLES_DISTRIBUTION, Informations about the store:);
    char string[BUFFER_SIZE];
		sprintf(string, ".... Store size: %ld", static_cast<long>(particlesNumber));
    LOGV(PARTICLES_DISTRIBUTION, string);
    if(isSet)
    {
        LOG(PARTICLES_DISTRIBUTION, .... The distribution is: set);
    }
    else
    {
        LOG(PARTICLES_DISTRIBUTION, .... The distribution is: null);
    }
    std::string tmp_str = ".... Store type: " + distr_type;
    LOGV(PARTICLES_DISTRIBUTION, tmp_str);
    char string2[BUFFER_SIZE];
    sprintf(string2,".... Store dimensions: %ld", static_cast<long>(distr_dims));
    LOGV(PARTICLES_DISTRIBUTION, string2);
    if(with_store)
    {
        LOG(PARTICLES_DISTRIBUTION, .... Outputing the store:);
        std::cout << (*store) << std::endl;
    }
}
#endif