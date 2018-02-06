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

#include <dataProcessing/procIslands.h>

ProcIslands::ProcIslands(ParticlesDistribution* set, std::string xp) : 
    part_set(set)
{
    LOG(PROC_ISLANDS, Instanciating ProcIslands);
    store = part_set->getStore();
    xpath = xp;
    moments = false;
    derived = false;
    grid_flag = true;
    readParameters();
  
    // Processing the boxes
    LOG(PROC_ISLANDS, Processing the islands...);
    process();
}

ProcIslands::~ProcIslands()
{
    LOG(PROC ISLANDS, Destructor called);
}

void
ProcIslands::process()
{
    // For each couple of planes
    for(int i = 0; i < n_processings; i++)
    {
        // Check for wich planes we are working and define the file name
        int plane1, plane2;
        if(planes.at(i) == "XXP")
        {
            plane1 = X_PLANE;
            plane2 = XP_PLANE;
            file.open((file_name + "XXP.dat").c_str());
            printHeaders("X - XP");
        }
        else if(planes.at(i) == "YYP")
        {
            plane1 = Y_PLANE;
            plane2 = YP_PLANE;
            file.open((file_name + "YYP.dat").c_str());
            printHeaders("Y - YP");
        }
        else if(planes.at(i) == "XY")
        {
            plane1 = X_PLANE;
            plane2 = Y_PLANE;
            file.open((file_name + "XY.dat").c_str());
            printHeaders("X - Y");
        }
        else if(planes.at(i) == "XYP")
        {
            plane1 = X_PLANE;
            plane2 = YP_PLANE;
            file.open((file_name + "XYP.dat").c_str());
            printHeaders("X - YP");
        }
        else if(planes.at(i) == "XPY")
        {
            plane1 = XP_PLANE;
            plane2 = Y_PLANE;
            file.open((file_name + "XPY.dat").c_str());
            printHeaders("XP - Y");
        }
        else if(planes.at(i) == "XPYP")
        {
            plane1 = XP_PLANE;
            plane2 = YP_PLANE;
            file.open((file_name + "XPYP.dat").c_str());
            printHeaders("XP - YP");
        }
        else
        {
            plane1 = plane2 = 0;
            throw("Invalid planes");
        }
        
        // Iterator on each islands
				int losses;
if(grid_flag)
{
 losses = computeLosses();
        for(std::vector<std::string>::iterator it = isl_list.at(i).begin(); it != isl_list.at(i).end(); it++)
        {
           if(s2n<int>(*it) > s2n<int>(boxes_x.at(i)) * s2n<int>(boxes_y.at(i)))
            {
                throw("Invalid box id");
            }
            else
            {
							std::vector<double> tmp = computeBoxGeometry(plane1, plane2, s2n<unsigned short int>(*it), s2n<double>(boxsize.at(i)), s2n<int>(boxes_x.at(i)), s2n<int>(boxes_y.at(i)));
               processBox(tmp);
            }
        } // Islands
}
else
{
	 losses = computeLosses();
	std::vector<double> tmp = computeBoxGeometry(plane1, plane2, box_corner_x, box_corner_xp, box_size_x, box_size_y);
	processBox(tmp);
	
}

	#ifdef DEBUG_FLAG
		std::cout << std::setw(SPACING) << "Total:" 
		          << std::setw(SPACING) << part_set->getParticlesNumber() 
		          << std::endl;
	  std::cout << std::setw(SPACING) << "Losses:" << std::setw(SPACING) << losses << std::endl;
	#endif
	file << "Total:" << std::setw(SPACING) << part_set->getParticlesNumber() 
	     << std::endl;
	file << "Losses:" << std::setw(SPACING) << losses << std::endl;
	
	// Close the file
  file.close();
        
    } // Planes
}

void 
ProcIslands::readParameters()
{
    LOG(PROC_ISLANDS, Reading the input parameters...);
    // Find the number of processing
    n_processings = XMLInputParser::getInstance().howManyFoundFromPath(xpath + "/analysis");
    
		if(XMLInputParser::getInstance().getAttributeTextFromPath(xpath + "/analysis/boxes", "grid") == "false")
		{
			grid_flag = false;
			box_corner_x = s2n<double>(XMLInputParser::getInstance().getFirstTextElementFromPath(xpath + "/analysis/boxes/box/x"));
			box_corner_xp = s2n<double>(XMLInputParser::getInstance().getFirstTextElementFromPath(xpath + "/analysis/boxes/box/xp"));
			box_size_x = s2n<double>(XMLInputParser::getInstance().getFirstTextElementFromPath(xpath + "/analysis/boxes/box/size_x"));
			box_size_y = s2n<double>(XMLInputParser::getInstance().getFirstTextElementFromPath(xpath + "/analysis/boxes/box/size_y"));
			
			// Planes
		  planes = XMLInputParser::getInstance().getTextElementsFromPath(xpath + "/analysis/planes"); 
		}
		else if(XMLInputParser::getInstance().getAttributeTextFromPath(xpath + "/analysis/boxes", "grid") == "true")
		{
			grid_flag = true;
			// Boxsizes
		  boxsize = XMLInputParser::getInstance().getTextElementsFromPath(xpath + "/analysis/boxes/size");

		  // Boxes horizontal
		  boxes_x = XMLInputParser::getInstance().getTextElementsFromPath(xpath + "/analysis/boxes/horizontal");

		  // Boxes vertical
		  boxes_y = XMLInputParser::getInstance().getTextElementsFromPath(xpath + "/analysis/boxes/vertical");

		  // Planes
		  planes = XMLInputParser::getInstance().getTextElementsFromPath(xpath + "/analysis/planes");
		
		  // Islands
      std::vector<xercesc::DOMNode*> isl_list_nodes = XMLInputParser::getInstance().getNodesFromPath(xpath + "/analysis");
      for(std::vector<xercesc::DOMNode*>::iterator i = isl_list_nodes.begin(); i != isl_list_nodes.end(); i++)
      {
        isl_list.push_back(XMLInputParser::getInstance().getTextElementsFromNode((*i), "islands/i"));
      }
    
      #ifdef DEBUG_FLAG
      for(std::vector<std::vector<std::string> >::iterator i = isl_list.begin(); i != isl_list.end(); i++)
      {
        LOG(PROC ISLANDS, List:);
        for(std::vector<std::string>::iterator j = (*i).begin(); j != (*i).end(); j++)
        {
            LOGV(PROC ISLANDS, (*j));
        }
      }
      #endif
		}

    // Flags
    if(XMLInputParser::getInstance().isFoundFromPath(xpath + "/moments"))
    {
        moments = true;
    }
    if(XMLInputParser::getInstance().isFoundFromPath(xpath + "/derived"))
    {
        derived = true;
    }
        
    // Take the generic file name
    file_name = XMLInputParser::getInstance().getAttributeTextFromPath(xpath, "file");    

}

int
ProcIslands::computeLosses()
{
  int losses = 0;
  for(int j = 0; j < part_set->getParticlesNumber();++j)
    {
      if(part_set->getDims() == TWO_DIMENSIONS)
	{
	  if((*store)(j,0) == LOST || (*store)(j,1) == LOST)
	    {
	      losses++;
	    }
	}
      else if(part_set->getDims() == FOUR_DIMENSIONS)
	{
	  if((*store)(j,0) == LOST && (*store)(j,1) == LOST && (*store)(j,2) == LOST && (*store)(j,3) == LOST)
	    {
	     losses++;
	    }
	}

    }
  return losses;
}

std::vector<double>
ProcIslands::computeBoxGeometry(int plane1, int plane2, double box_corner_x, double box_corner_xp, double box_size_x, double box_size_y)
{
	  // Compute the box boundaries
		double min_y = box_corner_xp - box_size_y;
		double max_y = box_corner_xp;
		double min_x = box_corner_x;
		double max_x = box_corner_x + box_size_x;
    // Compute the center of the box
    double center_plane1 = min_x + box_size_x/2;
    double center_plane2 = min_y + box_size_y/2;
		unsigned short int box_id = 0;

		std::vector<double> result;
		result.push_back(min_y);
		result.push_back(max_y);
		result.push_back(min_x);
		result.push_back(max_x);
		result.push_back(center_plane1);
		result.push_back(center_plane2);
		result.push_back(plane1);
		result.push_back(plane2);
		result.push_back(box_id);	
		
		return result;
}

std::vector<double>
ProcIslands::computeBoxGeometry(int plane1, int plane2, unsigned short int box_id, double box_size, int b_x, int b_y)
{
	  // Compute the box boundaries
    double min_y = -box_size+(2*box_size/b_x)*static_cast<int>((box_id-1)/b_x);
    double max_y = min_y + 2*box_size/b_x;
    double min_x = -box_size+(2*box_size/b_y)*((box_id-1)%b_y);
    double max_x = min_x + 2*box_size/b_x;
    
    // Compute the center of the box
    double center_plane1 =min_x + box_size/b_x;
    double center_plane2 =min_y + box_size/b_y;

		std::vector<double> result;
		result.push_back(min_y);
		result.push_back(max_y);
		result.push_back(min_x);
		result.push_back(max_x);
		result.push_back(center_plane1);
		result.push_back(center_plane2);
		result.push_back(plane1);
		result.push_back(plane2);
		result.push_back(box_id);
		
		return result;
}

void
ProcIslands::processBox(std::vector<double> geometry)
{    
    /*
    std::cout << geometry.at(0) << std::endl;
	  std::cout << geometry.at(1) << std::endl;
	  std::cout << geometry.at(2) << std::endl;
	  std::cout << geometry.at(3) << std::endl;
  	std::cout << geometry.at(4) << std::endl;
  	std::cout << geometry.at(5) << std::endl;
  	std::cout << geometry.at(6) << std::endl;
  	std::cout << geometry.at(7) << std::endl;
  	std::cout << geometry.at(8) << std::endl;
    */
							
    // Compute the box boundaries
		double min_y = geometry.at(0);
    double max_y = geometry.at(1);
    double min_x = geometry.at(2);
    double max_x = geometry.at(3);
    
		std::cout << "Min x : " << min_x << std::endl;
		std::cout << "Max x : " << max_x << std::endl;
		std::cout << "Min y : " << min_y << std::endl;
		std::cout << "Max y : " << max_y << std::endl;

    // Compute the center of the box
    double center_plane1 = geometry.at(4);
    double center_plane2 = geometry.at(5);

		int plane1 = geometry.at(6);
		int plane2 = geometry.at(7);
		unsigned short int box_id = geometry.at(8);

    // Initialize variables
    int box_count = 0;
    double mom_x = 0.0;
    double mom_xp = 0.0;
    double mom_x2 = 0.0;
    double mom_xp2 = 0.0;
    double mom_xxp = 0.0;
    double mom_x3 = 0.0;
    double mom_x2xp = 0.0;
    double mom_xxp2 = 0.0;
    double mom_xp3 = 0.0;
    double mom_x4 = 0.0;
    double mom_x3xp = 0.0;
    double mom_x2xp2 = 0.0;
    double mom_xxp3 = 0.0;
    double mom_xp4 = 0.0;
    
    // Rescaled coordinates
    double rescaled1 = 0.0;
    double rescaled2 = 0.0;
    
    // Find the particles in the box
    for(int j=0;j<part_set->getParticlesNumber();++j)
    {
        if((*store)(j,plane1) <= max_x && (*store)(j,plane2) <= max_y && (*store)(j,plane1) >= min_x && (*store)(j,plane2) >= min_y)
	{
           box_count++;
           mom_x += (*store)(j,plane1);
	   mom_xp += (*store)(j,plane2);
	} // If
    } // For
    mom_x /= box_count;
    mom_xp /= box_count;
            if(moments)
            {
    for(int j=0;j<part_set->getParticlesNumber();++j)
	{
        if((*store)(j,plane1) < max_x && (*store)(j,plane2) < max_y && (*store)(j,plane1) >= min_x && (*store)(j,plane2) >= min_y)
        {
	      rescaled1 = (*store)(j,plane1)-mom_x;
	      rescaled2 = (*store)(j,plane2)-mom_xp;

                // For every moment
                mom_x2 += rescaled1*rescaled1;
                mom_xp2 +=  rescaled2*rescaled2;
                mom_xxp +=  rescaled1*rescaled2;
                mom_x3 += rescaled1*rescaled1*rescaled1;
                mom_x2xp += rescaled1*rescaled1*rescaled2;
                mom_xxp2 += rescaled1*rescaled2*rescaled2;
                mom_xp3 += rescaled2*rescaled2*rescaled2;
                mom_x4 += rescaled1*rescaled1*rescaled1*rescaled1;
                mom_x3xp += rescaled1*rescaled1*rescaled1*rescaled2;
                mom_x2xp2 += rescaled1*rescaled1*rescaled2*rescaled2;
                mom_xxp3 += rescaled1*rescaled2*rescaled2*rescaled2;
                mom_xp4 += rescaled2*rescaled2*rescaled2*rescaled2;
            }
        }
	} 
	// Normalize the moments
	if(moments)
    {
        mom_x2 /= box_count;
        mom_xxp /= box_count;
        mom_xp2 /= box_count;
        mom_x3 /= box_count;
        mom_x2xp /= box_count;
        mom_xxp2 /= box_count;
        mom_xp3 /= box_count;
        mom_x4 /= box_count;
        mom_x3xp /= box_count;
        mom_x2xp2 /= box_count;
        mom_xxp3 /= box_count;
        mom_xp4 /= box_count;
    }
    
    // Compute derived quantities
    double emit_x = 0.0;
    double halo = 0.0;
    if(moments && derived)
    {
        emit_x = sqrt(mom_x2*mom_xp2-mom_xxp*mom_xxp);
        halo = (mom_x4/(mom_x2*mom_x2))-2.0; // Definition of the halo parameter
    }
    
    // Output the results
    unsigned short int spacing = SPACING;
    unsigned short int precision = PRECISION;
    #ifdef DEBUG_FLAG
    std::cout << std::setw(spacing) << std::setprecision(precision) << std::scientific << box_id
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << box_count
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << emit_x
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << halo
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x4
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x3xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2xp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp4
              << std::endl;
    #endif
    file      << std::setw(spacing) << std::setprecision(precision) << std::scientific << box_id
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << box_count
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << emit_x
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << halo
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x4
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x3xp
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_x2xp2
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xxp3
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << mom_xp4
              << std::endl;
    
}

void
ProcIslands::printHeaders(std::string p)
{
    // For a nicely formated output
    unsigned short int spacing = 13;
    #ifdef DEBUG_FLAG
    std::cout << "Islands processing for the planes " + p << std::endl;
    std::cout << std::setw(spacing) << "Isl"
              << std::setw(spacing) << "Nb"
              << std::setw(spacing) << "Emit_x"
              << std::setw(spacing) << "Halo"
              << std::setw(spacing) << "<x>"
              << std::setw(spacing) << "<x'>"
              << std::setw(spacing) << "<x2>"
              << std::setw(spacing) << "<xx'>"
              << std::setw(spacing) << "<x'2>"
              << std::setw(spacing) << "<x3>"
              << std::setw(spacing) << "<x2x'>"
              << std::setw(spacing) << "<xx'2>"
              << std::setw(spacing) << "<x'3>"
              << std::setw(spacing) << "<x4>"
              << std::setw(spacing) << "<x3x'>"
              << std::setw(spacing) << "<x2x'2>"
              << std::setw(spacing) << "<xx'3>"
              << std::setw(spacing) << "<x'4>" << std::endl; 
    #endif
    file << "Islands processing for the planes " + p << std::endl;
    file << std::setw(spacing) << "Isl"
         << std::setw(spacing) << "Nb"
         << std::setw(spacing) << "Emit_x"
         << std::setw(spacing) << "Halo"
         << std::setw(spacing) << "<x>"
         << std::setw(spacing) << "<x'>"
         << std::setw(spacing) << "<x2>"
         << std::setw(spacing) << "<xx'>"
         << std::setw(spacing) << "<x'2>"
         << std::setw(spacing) << "<x3>"
         << std::setw(spacing) << "<x2x'>"
         << std::setw(spacing) << "<xx'2>"
         << std::setw(spacing) << "<x'3>"
         << std::setw(spacing) << "<x4>"
         << std::setw(spacing) << "<x3x'>"
         << std::setw(spacing) << "<x2x'2>"
         << std::setw(spacing) << "<xx'3>"
         << std::setw(spacing) << "<x'4>" << std::endl; 
}

