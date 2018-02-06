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

#include <ParametersEvolution.h>

ParametersEvolution::ParametersEvolution(std::string s, std::string varpath)
{
	LOG(PARAMETERS EVOLUTION, Instanciating ParamatersEvolution);
	
	// Read some parameters
	turns = XMLInputParser::getInstance().getFirstTexti(s);
	variable_path = varpath;
	LOG(PARAMETER EVOLUTION, ParamatersEvolution instanciated.);
}

ParametersEvolution::~ParametersEvolution()
{
    LOG(PARAMETERS EVOLUTION, Destructor called.);
}

char
ParametersEvolution::readEvolvingVariable(std::string p, CurveParameters& param)
{
    return readEvolvingVariable(0, p, &param);
}

char
ParametersEvolution::readEvolvingVariable(int kick_id, std::string p, CurveParameters param[])
{
    std::string path = variable_path + i2s(kick_id) + "/" + p;
    // Check if the parameter is OK in the input file
    if(XMLInputParser::getInstance().isFoundFromPath(path))
    {		
        // Polynomial curve type
        if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "polyn_slope")
        {
            std::string tmp_path = path + "/start";
            param[kick_id].start = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/step";
            param[kick_id].step = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/step_t";
            param[kick_id].step_t = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/step_slope";
            param[kick_id].step_slope = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/end";
            param[kick_id].end = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/power1";
            param[kick_id].power1 = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/power2";
            param[kick_id].power2 = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            param[kick_id].curve_type = 'p';
            // In this case, also compute the derived parameters
            computeDerivedParametersPolynomialCurve(param[kick_id]);
         }
         // Linear curve type
         else if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "linear")
         {
            std::string tmp_path = path + "/start";
            param[kick_id].start = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/end";
            param[kick_id].end = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            param[kick_id].inc = (param[kick_id].end-param[kick_id].start)/turns;
            param[kick_id].curve_type = 'l';
         }
         // Fixed value
         else if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "fixed")
         {
            std::string tmp_path = path + "/value";
            param[kick_id].value = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            param[kick_id].curve_type = 'f';
         }
         // Chirp
         else if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "chirp")
         {
            std::string tmp_path = path + "/repetitions";
            param[kick_id].repetitions = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnStart";
            param[kick_id].turnStart = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsLength";
            param[kick_id].turnsLength = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsInterval";
            param[kick_id].turnsInterval = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/frequencyMax";
            param[kick_id].frequencyMax = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/frequencyMin";
            param[kick_id].frequencyMin = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/amplitude";
            param[kick_id].amplitude = XMLInputParser::getInstance().getFirstTextd(tmp_path);       
            param[kick_id].curve_type = 'c';
            // Compute the derived parameters
            computeDerivedParametersChirpCurve(param[kick_id]);
         }
         // Sin
         else if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "sin")
         {
            std::string tmp_path = path + "/repetitions";
            param[kick_id].repetitions = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnStart";
            param[kick_id].turnStart = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsLength";
            param[kick_id].turnsLength = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsInterval";
            param[kick_id].turnsInterval = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/frequency";
            param[kick_id].frequency = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            tmp_path = path + "/amplitude";
            param[kick_id].amplitude = XMLInputParser::getInstance().getFirstTextd(tmp_path);       
            param[kick_id].curve_type = 's';
         }
         // File
         else if(XMLInputParser::getInstance().getAttributeTextFromPath(path, "curveType") == "file")
         {
            std::string tmp_path = path + "/repetitions";
            param[kick_id].repetitions = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnStart";
            param[kick_id].turnStart = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsLength";
            param[kick_id].turnsLength = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/turnsInterval";
            param[kick_id].turnsInterval = XMLInputParser::getInstance().getFirstTexti(tmp_path);
            tmp_path = path + "/file";
            param[kick_id].file = XMLInputParser::getInstance().getFirstTextElementFromPath(tmp_path);    
            tmp_path = path + "/amplitude";
            param[kick_id].amplitude = XMLInputParser::getInstance().getFirstTextd(tmp_path);
            param[kick_id].curve_type = 'd';
            computeDerivedParametersFileCurve(param[kick_id]);
         }
         else
         {
            throw("Invalid curveType !");
         }
    }
    else
    {
        throw("Invalid parameter !");        
    }
    return param[kick_id].curve_type;
}

void
ParametersEvolution::computeDerivedParametersPolynomialCurve(CurveParameters& p)
{
    /*
       !!! WARNING !!!
       This method is a bit confusing
       It computes the polynomial coefficients using the inverse of the matrix of the constraints
       See the doc and the Mathematica notebook for more explanations
   */
    
    // First part of the curve
    if(p.power1 == 1)
    {
        p.a0 = p.start;
        p.a1 = (p.step-p.start)/p.step_t;
        p.aq = 0;
    }
    else if(p.power1 > 1)
    {
        p.a0 = p.start;
        double tmp = pow(p.step_t,p.power1);
        double tmp2 = -tmp+tmp*p.power1;
        p.a1 = ((-p.step_slope*tmp)-((tmp/p.step_t)*p.start*p.power1)+((tmp/p.step_t)*p.step*p.power1))/tmp2;
        p.aq = ((p.step_slope*p.step_t)/(tmp2))+((p.start)/(tmp2))-((p.step)/(tmp2));
    }
    // Second part of the curve
    if(p.power2 == 1)
    {
        p.b0 = p.step;
        p.b1 = (p.end-p.step)/(turns-p.step_t);
        p.bq = 0;
    }
    else if(p.power2 > 1)
    {
        int rturns = turns - p.step_t;
        p.b0 = p.step;
        p.b1 = p.step_slope;
        double tmp3 = pow(rturns,-p.power2);
        p.bq = tmp3*p.end - tmp3 * p.step - pow(rturns,1-p.power2) * p.step_slope;
    }
    
    #ifdef DEBUG_FLAG
        std::cout << "Computed polynomial parameters:" << std::endl;
        std::cout << std::setw(SPACING) << "a0" 
				  << std::setw(SPACING) << "a1" 
				  << std::setw(SPACING) << "aq"
				  << std::setw(SPACING) << "b0" 
				  << std::setw(SPACING) << "b1"
				  << std::setw(SPACING) << "bq" 
				  << std::setw(SPACING) << std::endl;
        std::cout << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.a0
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.a1
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.aq
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.b0
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.b1
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.bq << std::endl;
    #endif
}

void
ParametersEvolution::computeDerivedParametersChirpCurve(CurveParameters& p)
{   
    // State variables of the curve
	p.t_on = 0;
	p.is_on = false;
	p.index = 0;
	
    // Calculate the frequency increment of the chirp
    double revolution_frequency_ps_14GeV = 476.190; // That's really bad...
    p.frequencyIncrement = 2*M_PI*(p.frequencyMax/revolution_frequency_ps_14GeV)/p.turnsLength;
    
    #ifdef DEBUG_FLAG
        std::cout << "Computed chirp parameters:" << std::endl;
        std::cout << std::setw(SPACING) << "freq_rev" 
                  << std::setw(SPACING) << "freq_min"
                  << std::setw(SPACING) << "freq_max"
                  << std::setw(SPACING) << "freq_inc" << std::endl;
        std::cout << std::setw(SPACING) << std::setprecision(5) << std::scientific << revolution_frequency_ps_14GeV
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.frequencyMin
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.frequencyMax
                  << std::setw(SPACING) << std::setprecision(5) << std::scientific << p.frequencyIncrement << std::endl;
    #endif
}

void
ParametersEvolution::computeDerivedParametersFileCurve(CurveParameters& p)
{
    p.t_on = 0;
    p.is_on = false;
    p.index = 0;
    std::ifstream file(p.file.c_str(), std::ios::in);
    if(!file.is_open())
    {
        throw("File couldn't be opened !");
    }
    p.values = new double[p.turnsLength];
    std::string tmp;
    for(int i = 0; i < p.turnsLength; i++)
    {
        if(file.eof())
        {
					WARN(PARAMETERS EVOLUTION, File is too short: filling the last value.);
					p.values[i] = p.values[i-1];
        }
				else
				{        
					file >> tmp;
          p.values[i] = s2n<double>(tmp) * p.amplitude;

			}
    }
}

char
ParametersEvolution::setToZero(int kick_id, CurveParameters param[])
{
    param[kick_id].value = 0.0;
    param[kick_id].curve_type = 'f';
    return 'f';
}


