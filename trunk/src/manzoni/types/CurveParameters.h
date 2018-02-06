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

#ifndef CURVE_PARAMETERS_H
#define CURVE_PARAMETERS_H

#include <string>

// Simple container for the curve parameters
class CurveParameters
{
    public:
        //
        // Curve parameters
        //
        char curve_type;
        //
        // Fixed value curve
        //
        double value;
        //
        // Linear curve
        //
        double start; // Also for polynomial
        double end; // Also for polynomial
        double inc;
        //
        // Polynomial
        //
        double step;
        int step_t;
        double step_slope;
        int power1;
        int power2;
        // Polynomial parameters of the curve
        double a0;
        double a1;
        double aq;
        double b0;
        double b1;
        double bq;
        //
        // Windowed
        //
        int repetitions;
        int turnStart;
        int turnsLength;
        int turnsInterval;
        // File
        std::string file;
        double* values; // Should be a std::vector
				//std::vector<double> values;
        // Sin
        double frequency;
        // Chirp and sin
        double frequencyMax;
        double frequencyMin;
        double amplitude;
        double frequencyIncrement;
        int t_on;
        int index;
	    	bool is_on;
};

#endif

