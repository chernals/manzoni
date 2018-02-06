#ifndef PROC_DISTR_PROFILE_H
#define PROC_DISTR_PROFILE_H

#include <string>
#include <list>
#include <vector>

#include <ParticlesSet.h>
#include <root.h>

class ProcDistrProfile
{
	private:
	    //
	    // Data
	    //
	    store_t* store;
	    ParticlesSet* _part_set;
	    const std::string xpath_;
	    std::string _file;
	    std::string _title;
	    unsigned short int _x_divs;
        unsigned short int _y_divs;
        std::vector<std::string> profiles_;
        std::list<TH1D*> _profiles_list;
        double box_size_;
        TCanvas* _canvas;
        unsigned short int _profile_counter;
        //
        // Methods
        //
        void readParameters();

	public:
		ProcDistrProfile(ParticlesSet*, std::string);
		~ProcDistrProfile();
		void createProfile(std::string, unsigned short int);
		void draw();
		void save();
};

#endif

